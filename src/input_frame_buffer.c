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

#include "input_frame_buffer.h"

#include "encoder.h"
#include "encoderstate.h"
#include "image.h"


void uvg_init_input_frame_buffer(input_frame_buffer_t *input_buffer)
{
  FILL(input_buffer->pic_buffer, 0);
  FILL(input_buffer->pts_buffer, 0);
  input_buffer->num_in = 0;
  input_buffer->num_out = 0;
  input_buffer->delay = 0;
  input_buffer->gop_skipped = 0;
}

/**
 * \brief Pass an input frame to the encoder state.
 *
 * Returns the image that should be encoded next if there is a suitable
 * image available.
 *
 * The caller must not modify img_in after calling this function.
 *
 * \param buf         an input frame buffer
 * \param state       a main encoder state
 * \param img_in      input frame or NULL
 * \param first_done  whether the first frame has been done,
 *                    needed for the OBA rc
 * \return        pointer to the next picture, or NULL if no picture is
 *                available
 */
uvg_picture* uvg_encoder_feed_frame(input_frame_buffer_t *buf,
                                    encoder_state_t *const state,
                                    uvg_picture *const img_in, 
                                    int first_done)
{
  const encoder_control_t* const encoder = state->encoder_control;
  const uvg_config* const cfg = &encoder->cfg;

  const int gop_buf_size = 3 * cfg->gop_len;

  bool is_closed_gop = false;

  // Check for closed gop, we need an extra frame in the buffer in this case
  if (!cfg->open_gop && cfg->intra_period > 0 && cfg->gop_len > 0) is_closed_gop = true;

  if (cfg->gop_len == 0 || cfg->gop_lowdelay) {
    // No reordering of output pictures necessary.

    if (img_in == NULL) return NULL;

    img_in->dts = img_in->pts;
    state->frame->gop_offset = 0;
    if (cfg->gop_len > 0) {
      // Using a low delay GOP structure.
      uint64_t frame_num = buf->num_out;
      if (cfg->intra_period) {
        frame_num %= cfg->intra_period;
      }
      state->frame->gop_offset = (frame_num + cfg->gop_len - 1) % cfg->gop_len;
    }
    buf->num_in++;
    buf->num_out++;
    return uvg_image_copy_ref(img_in);
  }
  
  if (img_in != NULL) {
    // Index of the next input picture, in range [-1, +inf). Values
    // i and j refer to the same indices in buf->pic_buffer iff
    // i === j (mod gop_buf_size).
    int64_t idx_in = buf->num_in - 1;

    // Index in buf->pic_buffer and buf->pts_buffer.
    int buf_idx = (idx_in + gop_buf_size) % gop_buf_size;

    // Save the input image in the buffer.
    assert(buf_idx >= 0 && buf_idx < gop_buf_size);
    assert(buf->pic_buffer[buf_idx] == NULL);
    buf->pic_buffer[buf_idx] = uvg_image_copy_ref(img_in);
    buf->pts_buffer[buf_idx] = img_in->pts;
    buf->num_in++;

    if (buf->num_in < cfg->gop_len + is_closed_gop ? 1 : 0) {
      // Not enough frames to start output.
      return 0;

    } else if (buf->num_in == cfg->gop_len + is_closed_gop ? 1 : 0) {
      // Now we known the PTSs that are needed to compute the delay.
      buf->delay = buf->pts_buffer[gop_buf_size - 1] - img_in->pts;
    }
  }

  if (buf->num_out == buf->num_in) {
    // All frames returned.
    return NULL;
  }

  if (img_in == NULL && buf->num_in < cfg->gop_len + is_closed_gop ? 1 : 0) {
    // End of the sequence but we have less than a single GOP of frames. Use
    // the difference between the PTSs of the first and the last frame as the
    // delay.
    int first_pic_idx = gop_buf_size - 1;
    int last_pic_idx  = (buf->num_in - 2 + gop_buf_size) % gop_buf_size;
    buf->delay = buf->pts_buffer[first_pic_idx] - buf->pts_buffer[last_pic_idx];
  }

  // Index of the next output picture, in range [-1, +inf). Values
  // i and j refer to the same indices in buf->pic_buffer iff
  // i === j (mod gop_buf_size).
  int64_t idx_out;

  // DTS of the output picture.
  int64_t dts_out;

  // Number of the next output picture in the GOP.
  int gop_offset;

  if (buf->num_out == 0) {
    // Output the first frame.
    idx_out = -1;
    dts_out = buf->pts_buffer[gop_buf_size - 1] + buf->delay;
    gop_offset = 0; // highest quality picture

  } else if(first_done) {
    gop_offset = (buf->num_out - 1) % cfg->gop_len;
    
    // For closed gop, calculate the gop_offset again
    if (!cfg->open_gop && cfg->intra_period > 0) {
      // Offset the GOP position for each extra I-frame added to the structure
      // in closed gop case
      int32_t num_extra_frames = (int32_t)((buf->num_out - 1) / (cfg->intra_period + 1));
      gop_offset = (buf->num_out - 1 - num_extra_frames) % cfg->gop_len;
    }

    // Index of the first picture in the GOP that is being output.
    int32_t gop_start_idx = (int32_t)(buf->num_out - 1 - gop_offset);

    // Skip pictures until we find an available one.
    gop_offset += buf->gop_skipped;

    // Every closed-gop IRAP handled here
    if (is_closed_gop && (!cfg->open_gop && ((buf->num_out - 1) % (cfg->intra_period + 1)) == cfg->intra_period)) {
      idx_out = gop_start_idx;
    } else {
      for (;;) {
        assert(gop_offset < cfg->gop_len + is_closed_gop ? 1 : 0);
        idx_out = gop_start_idx + cfg->gop[gop_offset].poc_offset - 1;
        if (idx_out < (int64_t)buf->num_in - 1) {
          // An available picture found.
          break;
        }
        buf->gop_skipped++;
        gop_offset++;
      }
    }

    if (buf->num_out < cfg->gop_len - 1) {
      // This picture needs a DTS that is less than the PTS of the first
      // frame so the delay must be applied.
      int32_t dts_idx = (int32_t)(buf->num_out - 1);
      dts_out = buf->pts_buffer[dts_idx % gop_buf_size] + buf->delay;
    } else {
      int32_t dts_idx = (int32_t)(buf->num_out - (cfg->gop_len - 1));
      dts_out = buf->pts_buffer[dts_idx % gop_buf_size] - 1;
    }
  }
  else {
    return NULL;
  }

  // Index in buf->pic_buffer and buf->pts_buffer.
  int buf_idx = (idx_out + gop_buf_size) % gop_buf_size;

  uvg_picture* next_pic = buf->pic_buffer[buf_idx];
  assert(next_pic != NULL);
  next_pic->dts = dts_out;
  buf->pic_buffer[buf_idx] = NULL;
  state->frame->gop_offset = gop_offset;

  buf->num_out++;
  return next_pic;
}
