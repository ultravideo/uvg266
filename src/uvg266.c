/*****************************************************************************
 * This file is part of uvg266 VVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

#include "uvg266.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bitstream.h"
#include "cfg.h"
#include "checkpoint.h"
#include "encoder.h"
#include "encoder_state-bitstream.h"
#include "encoder_state-ctors_dtors.h"
#include "encoderstate.h"
#include "global.h"
#include "image.h"
#include "input_frame_buffer.h"
#include "uvg266_internal.h"
#include "strategyselector.h"
#include "threadqueue.h"
#include "videoframe.h"
#include "rate_control.h"


static void uvg266_close(uvg_encoder *encoder)
{
  if (encoder) {
    // The threadqueue must be stopped before freeing states.
    if (encoder->control) {
      uvg_threadqueue_stop(encoder->control->threadqueue);
    }

    if (encoder->states) {
      // Flush input frame buffer.
      uvg_picture *pic = NULL;
      while ((pic = uvg_encoder_feed_frame(&encoder->input_buffer,
                                           &encoder->states[0],
                                           NULL,
                                           1)) != NULL) {
        uvg_image_free(pic);
        pic = NULL;
      }

      for (unsigned i = 0; i < encoder->num_encoder_states; ++i) {
        uvg_encoder_state_finalize(&encoder->states[i]);
      }
    }
    FREE_POINTER(encoder->states);

    uvg_free_rc_data();
    // Discard const from the pointer.
    uvg_encoder_control_free((void*) encoder->control);
    encoder->control = NULL;
  }
  FREE_POINTER(encoder);
}


static uvg_encoder * uvg266_open(const uvg_config *cfg)
{
  uvg_encoder *encoder = NULL;

  //Initialize strategies
  // TODO: Make strategies non-global
  if (!uvg_strategyselector_init(cfg->cpuid, UVG_BIT_DEPTH)) {
    fprintf(stderr, "Failed to initialize strategies.\n");
    goto uvg266_open_failure;
  }

  encoder = calloc(1, sizeof(uvg_encoder));
  if (!encoder) {
    goto uvg266_open_failure;
  }

  encoder->control = uvg_encoder_control_init(cfg);
  if (!encoder->control) {
    goto uvg266_open_failure;
  }

  encoder->num_encoder_states = encoder->control->cfg.owf + 1;
  encoder->cur_state_num = 0;
  encoder->out_state_num = 0;
  encoder->frames_started = 0;
  encoder->frames_done = 0;

  // Assure that the rc data allocation was successful
  if(!uvg_get_rc_data(encoder->control)) {
    goto uvg266_open_failure;
  }

  uvg_init_input_frame_buffer(&encoder->input_buffer);

  encoder->states = calloc(encoder->num_encoder_states, sizeof(encoder_state_t));
  if (!encoder->states) {
    goto uvg266_open_failure;
  }

  for (unsigned i = 0; i < encoder->num_encoder_states; ++i) {
    encoder->states[i].encoder_control = encoder->control;
    if (!uvg_encoder_state_init(&encoder->states[i], NULL)) {
      goto uvg266_open_failure;
    }

    encoder->states[i].frame->QP = (int8_t)cfg->qp;
  }

  for (uint32_t i = 0; i < encoder->num_encoder_states; ++i) {
    if (i == 0) {
      encoder->states[i].previous_encoder_state = &encoder->states[encoder->num_encoder_states - 1];
    } else {
      encoder->states[i].previous_encoder_state = &encoder->states[(i - 1) % encoder->num_encoder_states];
    }
    uvg_encoder_state_match_children_of_previous_frame(&encoder->states[i]);
  }

  encoder->states[encoder->cur_state_num].frame->num = -1;

  return encoder;

uvg266_open_failure:
  uvg266_close(encoder);
  return NULL;
}


static void set_frame_info(uvg_frame_info *const info, const encoder_state_t *const state)
{
  info->poc = state->frame->poc,
  info->qp = state->frame->QP;
  info->nal_unit_type = state->frame->pictype;
  info->slice_type = state->frame->slicetype;

  memset(info->ref_list[0], 0, 16 * sizeof(int));
  memset(info->ref_list[1], 0, 16 * sizeof(int));

  for (size_t i = 0; i < state->frame->ref_LX_size[0]; i++) {
    info->ref_list[0][i] = state->frame->ref->pocs[state->frame->ref_LX[0][i]];
  }

  for (size_t i = 0; i < state->frame->ref_LX_size[1]; i++) {
    info->ref_list[1][i] = state->frame->ref->pocs[state->frame->ref_LX[1][i]];
  }

  info->ref_list_len[0] = state->frame->ref_LX_size[0];
  info->ref_list_len[1] = state->frame->ref_LX_size[1];
}


static int uvg266_headers(uvg_encoder *enc,
                           uvg_data_chunk **data_out,
                           uint32_t *len_out)
{
  if (data_out) *data_out = NULL;
  if (len_out) *len_out = 0;

  bitstream_t stream;
  uvg_bitstream_init(&stream);

  uvg_encoder_state_write_parameter_sets(&stream, &enc->states[enc->cur_state_num]);

  // Get stream length before taking chunks since that clears the stream.
  if (len_out) *len_out = (uint32_t)(uvg_bitstream_tell(&stream) / 8);
  if (data_out) *data_out = uvg_bitstream_take_chunks(&stream);

  uvg_bitstream_finalize(&stream);
  return 1;
}


/**
* \brief Separate a single field from a frame.
*
* \param frame_in           input frame to extract field from
* \param source_scan_type   scan type of input material (0: progressive, 1:top field first, 2:bottom field first)
* \param field parity   
* \param field_out
*
* \return              1 on success, 0 on failure
*/
static int yuv_io_extract_field(const uvg_picture *frame_in, unsigned source_scan_type, unsigned field_parity, uvg_picture *field_out)
{
  if ((source_scan_type != 1) && (source_scan_type != 2)) return 0;
  if ((field_parity != 0)     && (field_parity != 1))     return 0;

  int32_t offset = 0;
  if (source_scan_type == 1) offset = field_parity ? 1 : 0;
  else if (source_scan_type == 2) offset = field_parity ? 0 : 1;  

  //Luma
  for (int32_t i = 0; i < field_out->height; ++i){
    uvg_pixel *row_in  = frame_in->y + MIN(frame_in->height - 1, 2 * i + offset) * frame_in->stride;
    uvg_pixel *row_out = field_out->y + i * field_out->stride;
    memcpy(row_out, row_in, sizeof(uvg_pixel) * frame_in->stride);
  }

  //Chroma
  for (int32_t i = 0; i < field_out->height / 2; ++i) {
    uvg_pixel *row_in = frame_in->u + MIN(frame_in->height / 2 - 1, 2 * i + offset) * frame_in->stride / 2;
    uvg_pixel *row_out = field_out->u + i * field_out->stride / 2;
    memcpy(row_out, row_in, sizeof(uvg_pixel) * frame_in->stride / 2);
  }

  for (int32_t i = 0; i < field_out->height / 2; ++i) {
    uvg_pixel *row_in = frame_in->v + MIN(frame_in->height / 2 - 1, 2 * i + offset) * frame_in->stride / 2;
    uvg_pixel *row_out = field_out->v + i * field_out->stride / 2;
    memcpy(row_out, row_in, sizeof(uvg_pixel) * frame_in->stride / 2);
  }

  return 1;
}


static int uvg266_encode(uvg_encoder *enc,
                          uvg_picture *pic_in,
                          uvg_data_chunk **data_out,
                          uint32_t *len_out,
                          uvg_picture **pic_out,
                          uvg_picture **src_out,
                          uvg_frame_info *info_out)
{
  if (data_out) *data_out = NULL;
  if (len_out) *len_out = 0;
  if (pic_out) *pic_out = NULL;
  if (src_out) *src_out = NULL;

  encoder_state_t *state = &enc->states[enc->cur_state_num];

  if (!state->frame->prepared) {
    uvg_encoder_prepare(state);
  }

  if (pic_in != NULL) {
    // FIXME: The frame number printed here is wrong when GOP is enabled.
    CHECKPOINT_MARK("read source frame: %d", state->frame->num + enc->control->cfg.seek);
  }

  uvg_picture* frame = uvg_encoder_feed_frame(
    &enc->input_buffer, state, pic_in,
    enc->frames_done || state->encoder_control->cfg.rc_algorithm != UVG_OBA
  );
  if (frame) {
    assert(state->frame->num == enc->frames_started);
    // Start encoding.
    uvg_encode_one_frame(state, frame);
    enc->frames_started += 1;
  }

  // If we have finished encoding as many frames as we have started, we are done.
  if (enc->frames_done == enc->frames_started) {
    return 1;
  }

  if (!state->frame->done) {
    // We started encoding a frame; move to the next encoder state.
    enc->cur_state_num = (enc->cur_state_num + 1) % (enc->num_encoder_states);
  }

  encoder_state_t *output_state = &enc->states[enc->out_state_num];
  if ((!output_state->frame->done &&
       (pic_in == NULL || enc->cur_state_num == enc->out_state_num)) ||
       (state->frame->num == 0  && state->encoder_control->cfg.rc_algorithm == UVG_OBA)) {

    uvg_threadqueue_waitfor(enc->control->threadqueue, output_state->tqj_bitstream_written);
    // The job pointer must be set to NULL here since it won't be usable after
    // the next frame is done.
    uvg_threadqueue_free_job(&output_state->tqj_bitstream_written);

    // Get stream length before taking chunks since that clears the stream.
    if (len_out) *len_out = (uint32_t)(uvg_bitstream_tell(&output_state->stream) / 8);
    if (data_out) *data_out = uvg_bitstream_take_chunks(&output_state->stream);
    if (pic_out) *pic_out = uvg_image_copy_ref(output_state->tile->frame->rec);
    if (src_out) *src_out = uvg_image_copy_ref(output_state->tile->frame->source);
    if (info_out) set_frame_info(info_out, output_state);

    output_state->frame->done = 1;
    output_state->frame->prepared = 0;
    enc->frames_done += 1;

    enc->out_state_num = (enc->out_state_num + 1) % (enc->num_encoder_states);
  }

  return 1;
}


static int uvg266_field_encoding_adapter(uvg_encoder *enc,
                                          uvg_picture *pic_in,
                                          uvg_data_chunk **data_out,
                                          uint32_t *len_out,
                                          uvg_picture **pic_out,
                                          uvg_picture **src_out,
                                          uvg_frame_info *info_out)
{
  if (enc->control->cfg.source_scan_type == UVG_INTERLACING_NONE) {
    // For progressive, simply call the normal encoding function.
    return uvg266_encode(enc, pic_in, data_out, len_out, pic_out, src_out, info_out);
  }

  // For interlaced, make two fields out of the input frame and call encode on them separately.
  encoder_state_t *state = &enc->states[enc->cur_state_num];
  uvg_picture *first_field = NULL, *second_field = NULL;
  struct {
    uvg_data_chunk* data_out;
    uint32_t len_out;
  } first = { 0, 0 }, second = { 0, 0 };

  if (pic_in != NULL) {
    first_field = uvg_image_alloc(state->encoder_control->chroma_format, state->encoder_control->in.width, state->encoder_control->in.height);
    if (first_field == NULL) {
      goto uvg266_field_encoding_adapter_failure;
    }
    second_field = uvg_image_alloc(state->encoder_control->chroma_format, state->encoder_control->in.width, state->encoder_control->in.height);
    if (second_field == NULL) {
      goto uvg266_field_encoding_adapter_failure;
    }

    yuv_io_extract_field(pic_in, pic_in->interlacing, 0, first_field);
    yuv_io_extract_field(pic_in, pic_in->interlacing, 1, second_field);
    
    first_field->pts = pic_in->pts;
    first_field->dts = pic_in->dts;
    first_field->interlacing = pic_in->interlacing;

    // Should the second field have higher pts and dts? It shouldn't affect anything.
    second_field->pts = pic_in->pts;
    second_field->dts = pic_in->dts;
    second_field->interlacing = pic_in->interlacing;
  }

  if (!uvg266_encode(enc, first_field, &first.data_out, &first.len_out, pic_out, NULL, info_out)) {
    goto uvg266_field_encoding_adapter_failure;
  }
  if (!uvg266_encode(enc, second_field, &second.data_out, &second.len_out, NULL, NULL, NULL)) {
    goto uvg266_field_encoding_adapter_failure;
  }

  uvg_image_free(first_field);
  uvg_image_free(second_field);

  // Concatenate bitstreams.
  if (len_out != NULL) {
    *len_out = first.len_out + second.len_out;
  }
  if (data_out != NULL) {
    *data_out = first.data_out;
    if (first.data_out != NULL) {
      uvg_data_chunk *chunk = first.data_out;
      while (chunk->next != NULL) {
        chunk = chunk->next;
      }
      chunk->next = second.data_out;
    }
  }

  if (src_out != NULL) {
    // TODO: deinterlace the fields to one picture.
  }

  return 1;

uvg266_field_encoding_adapter_failure:
  uvg_image_free(first_field);
  uvg_image_free(second_field);
  uvg_bitstream_free_chunks(first.data_out);
  uvg_bitstream_free_chunks(second.data_out);
  return 0;
}


static const uvg_api uvg_8bit_api = {
  .config_alloc = uvg_config_alloc,
  .config_init = uvg_config_init,
  .config_destroy = uvg_config_destroy,
  .config_parse = uvg_config_parse,

  .picture_alloc = uvg_image_alloc_420,
  .picture_free = uvg_image_free,

  .chunk_free = uvg_bitstream_free_chunks,

  .encoder_open = uvg266_open,
  .encoder_close = uvg266_close,
  .encoder_headers = uvg266_headers,
  .encoder_encode = uvg266_field_encoding_adapter,

  .picture_alloc_csp = uvg_image_alloc,
};


const uvg_api * uvg_api_get(int bit_depth)
{
  return &uvg_8bit_api;
}
