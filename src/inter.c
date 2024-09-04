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

#include "inter.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "encoder.h"
#include "imagelist.h"
#include "uvg_math.h"
#include "strategies/generic/picture-generic.h"
#include "strategies/strategies-ipol.h"
#include "videoframe.h"
#include "strategies/strategies-picture.h"


const int8_t uvg_g_imv_to_prec[4] = { 2, 0, -2, 1 }; // Convert AMVR imv to precision

typedef struct {
  const cu_info_t *a[2];
  const cu_info_t *b[3];
  const cu_info_t *c0;
  const cu_info_t *c1;

} merge_candidates_t;


static void inter_recon_frac_luma(const encoder_state_t * const state,
                                  const uvg_picture * const ref,
                                  int32_t xpos,
                                  int32_t ypos,
                                  int32_t block_width,
                                  int32_t block_height,
                                  const mv_t mv_param[2],
                                  yuv_t *out,
                                  unsigned out_stride)
{
  int mv_frac_x = (mv_param[0] & 15);
  int mv_frac_y = (mv_param[1] & 15);

  // Space for extrapolated pixels and the part from the picture.
  // Some extra for AVX2.
  // The extrapolation function will set the pointers and stride.
  uvg_pixel ext_buffer[UVG_IPOL_MAX_INPUT_SIZE_LUMA_SIMD];
  uvg_pixel *ext = NULL;
  uvg_pixel *ext_origin = NULL;
  int ext_s = 0;
  uvg_epol_args epol_args = {
    .src = ref->y,
    .src_w = ref->width,
    .src_h = ref->height,
    .src_s = ref->stride,
    .blk_x = state->tile->offset_x + xpos + (mv_param[0] >> INTERNAL_MV_PREC),
    .blk_y = state->tile->offset_y + ypos + (mv_param[1] >> INTERNAL_MV_PREC),
    .blk_w = block_width,
    .blk_h = block_height,
    .pad_l = UVG_LUMA_FILTER_OFFSET,
    .pad_r = UVG_EXT_PADDING_LUMA - UVG_LUMA_FILTER_OFFSET,
    .pad_t = UVG_LUMA_FILTER_OFFSET,
    .pad_b = UVG_EXT_PADDING_LUMA - UVG_LUMA_FILTER_OFFSET,
    .pad_b_simd = 1 // One row for AVX2
  };

  // Initialize separately. Gets rid of warning
  // about using nonstandard extension.
  epol_args.buf = ext_buffer;
  epol_args.ext = &ext;
  epol_args.ext_origin = &ext_origin;
  epol_args.ext_s = &ext_s;

  uvg_get_extended_block(&epol_args);
  uvg_sample_quarterpel_luma(state->encoder_control,
    ext_origin,
    ext_s,
    block_width,
    block_height,
    out->y,
    out_stride,
    mv_frac_x,
    mv_frac_y,
    mv_param);
}

static void inter_recon_frac_luma_hi(const encoder_state_t *const state,
  const uvg_picture *const ref,
  int32_t xpos,
  int32_t ypos,
  int32_t block_width,
  int32_t block_height,
  const mv_t mv_param[2],
  yuv_im_t *out,
  const unsigned out_stride)
{
  int mv_frac_x = (mv_param[0] & 15);
  int mv_frac_y = (mv_param[1] & 15);

  // Space for extrapolated pixels and the part from the picture.
  // Some extra for AVX2.
  // The extrapolation function will set the pointers and stride.
  uvg_pixel ext_buffer[UVG_IPOL_MAX_INPUT_SIZE_LUMA_SIMD];
  uvg_pixel *ext = NULL;
  uvg_pixel *ext_origin = NULL;
  int ext_s = 0;
  uvg_epol_args epol_args = {
    .src = ref->y,
    .src_w = ref->width,
    .src_h = ref->height,
    .src_s = ref->stride,
    .blk_x = state->tile->offset_x + xpos + (mv_param[0] >> INTERNAL_MV_PREC),
    .blk_y = state->tile->offset_y + ypos + (mv_param[1] >> INTERNAL_MV_PREC),
    .blk_w = block_width,
    .blk_h = block_height,
    .pad_l = UVG_LUMA_FILTER_OFFSET,
    .pad_r = UVG_EXT_PADDING_LUMA - UVG_LUMA_FILTER_OFFSET,
    .pad_t = UVG_LUMA_FILTER_OFFSET,
    .pad_b = UVG_EXT_PADDING_LUMA - UVG_LUMA_FILTER_OFFSET,
    .pad_b_simd = 1 // One row for AVX2
  };

  // Initialize separately. Gets rid of warning
  // about using nonstandard extension.
  epol_args.buf = ext_buffer;
  epol_args.ext = &ext;
  epol_args.ext_origin = &ext_origin;
  epol_args.ext_s = &ext_s;

  uvg_get_extended_block(&epol_args);
  uvg_sample_quarterpel_luma_hi(state->encoder_control,
    ext_origin,
    ext_s,
    block_width,
    block_height,
    out->y,
    out_stride,
    mv_frac_x,
    mv_frac_y,
    mv_param);
}

static void inter_recon_frac_chroma(const encoder_state_t *const state,
  const uvg_picture *const ref,
  int32_t pu_x,
  int32_t pu_y,
  int32_t pu_w,
  int32_t pu_h,
  const mv_t mv_param[2],
  yuv_t *out,
  const unsigned out_stride)
{
  int mv_frac_x = (mv_param[0] & 31);
  int mv_frac_y = (mv_param[1] & 31);

  // Take into account chroma subsampling
  unsigned pb_w = pu_w / 2;
  unsigned pb_h = pu_h / 2;

  // Space for extrapolated pixels and the part from the picture.
  // Some extra for AVX2.
  // The extrapolation function will set the pointers and stride.
  uvg_pixel ext_buffer[UVG_IPOL_MAX_INPUT_SIZE_CHROMA_SIMD];
  uvg_pixel *ext = NULL;
  uvg_pixel *ext_origin = NULL;
  int ext_s = 0;

  // Chroma U
  // Divisions by 2 due to 4:2:0 chroma subsampling
  uvg_epol_args epol_args = {
    .src = ref->u,
    .src_w = ref->width / 2,
    .src_h = ref->height / 2,
    .src_s = ref->stride / 2,
    .blk_x = (state->tile->offset_x + pu_x) / 2 + (mv_param[0] >> (INTERNAL_MV_PREC + 1) ),
    .blk_y = (state->tile->offset_y + pu_y) / 2 + (mv_param[1] >> (INTERNAL_MV_PREC + 1) ),
    .blk_w = pb_w,
    .blk_h = pb_h,
    .pad_l = UVG_CHROMA_FILTER_OFFSET,
    .pad_r = UVG_EXT_PADDING_CHROMA - UVG_CHROMA_FILTER_OFFSET,
    .pad_t = UVG_CHROMA_FILTER_OFFSET,
    .pad_b = UVG_EXT_PADDING_CHROMA - UVG_CHROMA_FILTER_OFFSET,
    .pad_b_simd = 3 // Three rows for AVX2
  };

  // Initialize separately. Gets rid of warning
  // about using nonstandard extension.
  epol_args.buf = ext_buffer;
  epol_args.ext = &ext;
  epol_args.ext_origin = &ext_origin;
  epol_args.ext_s = &ext_s;

  uvg_get_extended_block(&epol_args);
  uvg_sample_octpel_chroma(state->encoder_control,
    ext_origin,
    ext_s,
    pb_w,
    pb_h,
    out->u,
    out_stride,
    mv_frac_x,
    mv_frac_y,
    mv_param);

  // Chroma V
  epol_args.src = ref->v;
  uvg_get_extended_block(&epol_args);
  uvg_sample_octpel_chroma(state->encoder_control,
    ext_origin,
    ext_s,
    pb_w,
    pb_h,
    out->v,
    out_stride,
    mv_frac_x,
    mv_frac_y,
    mv_param);
}

static void inter_recon_frac_chroma_hi(const encoder_state_t *const state,
  const uvg_picture *const ref,
  int32_t pu_x,
  int32_t pu_y,
  int32_t pu_w,
  int32_t pu_h,
  const mv_t mv_param[2],
  yuv_im_t *out,
  const unsigned out_stride)
{
  int mv_frac_x = (mv_param[0] & 31);
  int mv_frac_y = (mv_param[1] & 31);

  // Take into account chroma subsampling
  unsigned pb_w = pu_w / 2;
  unsigned pb_h = pu_h / 2;

  // Space for extrapolated pixels and the part from the picture.
  // Some extra for AVX2.
  // The extrapolation function will set the pointers and stride.
  uvg_pixel ext_buffer[UVG_IPOL_MAX_INPUT_SIZE_CHROMA_SIMD];
  uvg_pixel *ext = NULL;
  uvg_pixel *ext_origin = NULL;
  int ext_s = 0;

  // Chroma U
  // Divisions by 2 due to 4:2:0 chroma subsampling
  uvg_epol_args epol_args = {
    .src = ref->u,
    .src_w = ref->width / 2,
    .src_h = ref->height / 2,
    .src_s = ref->stride / 2,
    .blk_x = (state->tile->offset_x + pu_x) / 2 + (mv_param[0] >> (INTERNAL_MV_PREC + 1) ),
    .blk_y = (state->tile->offset_y + pu_y) / 2 + (mv_param[1] >> (INTERNAL_MV_PREC + 1) ),
    .blk_w = pb_w,
    .blk_h = pb_h,
    .pad_l = UVG_CHROMA_FILTER_OFFSET,
    .pad_r = UVG_EXT_PADDING_CHROMA - UVG_CHROMA_FILTER_OFFSET,
    .pad_t = UVG_CHROMA_FILTER_OFFSET,
    .pad_b = UVG_EXT_PADDING_CHROMA - UVG_CHROMA_FILTER_OFFSET,
    .pad_b_simd = 3 // Three rows for AVX2
  };

  // Initialize separately. Gets rid of warning
  // about using nonstandard extension.
  epol_args.buf = ext_buffer;
  epol_args.ext = &ext;
  epol_args.ext_origin = &ext_origin;
  epol_args.ext_s = &ext_s;

  uvg_get_extended_block(&epol_args);
  uvg_sample_octpel_chroma_hi(state->encoder_control,
    ext_origin,
    ext_s,
    pb_w,
    pb_h,
    out->u,
    out_stride,
    mv_frac_x,
    mv_frac_y,
    mv_param);

  // Chroma V
  epol_args.src = ref->v;
  uvg_get_extended_block(&epol_args);
  uvg_sample_octpel_chroma_hi(state->encoder_control,
    ext_origin,
    ext_s,
    pb_w,
    pb_h,
    out->v,
    out_stride,
    mv_frac_x,
    mv_frac_y,
    mv_param);
}


/**
* \brief Copy from frame with extended border.
*
* \param ref_buf      pointer to the start of ref buffer
* \param ref_stride   stride of ref buffer
* \param ref_width    width of frame
* \param ref_height   height of frame
* \param rec_buf      pointer to the start of pu in rec buffer
* \param rec_stride   stride of rec buffer
* \param width        width of copied block
* \param height       height of copied block
* \param mv_in_frame  coordinates of copied block in frame coordinates
*/
static void inter_cp_with_ext_border(const uvg_pixel *ref_buf, int ref_stride,
                                     int ref_width, int ref_height,
                                     uvg_pixel *rec_buf, int rec_stride,
                                     int width, int height,
                                     const vector2d_t *mv_in_frame)
{
  for (int y = mv_in_frame->y; y < mv_in_frame->y + height; ++y) {
    for (int x = mv_in_frame->x; x < mv_in_frame->x + width; ++x) {
      vector2d_t in_frame = {
        CLIP(0, ref_width - 1, x),
        CLIP(0, ref_height - 1, y),
      };
      vector2d_t in_pu = {
        x - mv_in_frame->x,
        y - mv_in_frame->y,
      };
      int pu_index = in_pu.y * rec_stride + in_pu.x;
      int frame_index = in_frame.y * ref_stride + in_frame.x;
      rec_buf[pu_index] = ref_buf[frame_index];
    }
  }
}


/**
 * \brief Reconstruct an inter PU using uniprediction.
 *
 * \param state          encoder state
 * \param ref            picture to copy the data from
 * \param mv_param       motion vector
 * \param yuv_px         destination buffer for pixel precision
 * \param yuv_im         destination buffer for high-precision, or NULL if not needed
 * \param predict_luma   Enable or disable luma prediction for this call.
 * \param predict_chroma Enable or disable chroma prediction for this call.
 * \param cu_loc         Size and position of current PU/CU
*/
static unsigned inter_recon_unipred(
  const encoder_state_t * const state,
  const uvg_picture * const ref,
  int32_t out_stride_luma,
  const mv_t mv_param[2],
  yuv_t *yuv_px,
  yuv_im_t *yuv_im,
  bool predict_luma,
  bool predict_chroma,
  const cu_loc_t* const cu_loc)
{
  vector2d_t int_mv = { mv_param[0], mv_param[1] };

  uvg_change_precision_vector2d(INTERNAL_MV_PREC, 0, &int_mv);

  const int pu_x = cu_loc->x;
  const int pu_y = cu_loc->y;
  const int pu_w = cu_loc->width;
  const int pu_h = cu_loc->height;

  const vector2d_t int_mv_in_frame = {
    int_mv.x + pu_x + state->tile->offset_x,
    int_mv.y + pu_y + state->tile->offset_y
  };

  const bool int_mv_outside_frame = int_mv_in_frame.x < 0 ||
    int_mv_in_frame.y < 0 ||
    int_mv_in_frame.x + pu_w > ref->width ||
    int_mv_in_frame.y + pu_h > ref->height;

  // With 420, odd coordinates need interpolation.
  const bool fractional_chroma = (int_mv.x & 1) || (int_mv.y & 1);
  const bool fractional_luma = ((mv_param[0] & 15) || (mv_param[1] & 15));

  // Generate prediction for luma.
  if (predict_luma) {
    if (fractional_luma) {
      // With a fractional MV, do interpolation.
      if (state->encoder_control->cfg.bipred && yuv_im) {
        inter_recon_frac_luma_hi(state, ref,
          pu_x, pu_y,
          pu_w, pu_h,
          mv_param, yuv_im, out_stride_luma);
      }
      else {
        inter_recon_frac_luma(state, ref,
          pu_x, pu_y,
          pu_w, pu_h,
          mv_param, yuv_px, out_stride_luma);
      }
    } else {
      // With an integer MV, copy pixels directly from the reference.
      if (int_mv_outside_frame) {
        inter_cp_with_ext_border(ref->y, ref->stride,
          ref->width, ref->height,
          yuv_px->y, out_stride_luma,
          pu_w, pu_h,
          &int_mv_in_frame);
      }
      else {
        const int frame_mv_index = int_mv_in_frame.y * ref->stride + int_mv_in_frame.x;
        uvg_pixels_blit(&ref->y[frame_mv_index],
          yuv_px->y,
          pu_w, pu_h,
          ref->stride, out_stride_luma);
      }
    }
  }

  if (!predict_chroma) {
    return fractional_luma;
  }

  const unsigned out_stride_c = out_stride_luma / 2;

  // Generate prediction for chroma.
  if (fractional_luma || fractional_chroma) {
    // With a fractional MV, do interpolation.
    if (state->encoder_control->cfg.bipred && yuv_im) {
      inter_recon_frac_chroma_hi(state, ref,
                                    pu_x, pu_y,
                                    pu_w, pu_h, 
                                    mv_param, yuv_im, out_stride_c);
    } else {
      inter_recon_frac_chroma(state, ref,
                              pu_x, pu_y,
                              pu_w, pu_h,
                              mv_param, yuv_px, out_stride_c);
    }
  } else {
    // With an integer MV, copy pixels directly from the reference.
    const vector2d_t int_mv_in_frame_c = { int_mv_in_frame.x / 2, int_mv_in_frame.y / 2 };

    if (int_mv_outside_frame) {
      inter_cp_with_ext_border(ref->u, ref->stride / 2,
                               ref->width / 2, ref->height / 2,
                               yuv_px->u, out_stride_c,
                               pu_w / 2, pu_h / 2,
                               &int_mv_in_frame_c);
      inter_cp_with_ext_border(ref->v, ref->stride / 2,
                               ref->width / 2, ref->height / 2,
                               yuv_px->v, out_stride_c,
                               pu_w / 2, pu_h / 2,
                               &int_mv_in_frame_c);
    } else {
      const int frame_mv_index = int_mv_in_frame_c.y * ref->stride / 2 + int_mv_in_frame_c.x;

      uvg_pixels_blit(&ref->u[frame_mv_index],
                      yuv_px->u,
                      pu_w / 2, pu_h / 2,
                      ref->stride / 2, out_stride_c);
      uvg_pixels_blit(&ref->v[frame_mv_index],
                      yuv_px->v,
                      pu_w / 2, pu_h / 2,
                      ref->stride / 2, out_stride_c);
    }
  }

  return fractional_luma | ((fractional_luma || fractional_chroma) << 1);
}
/**
 * \brief Reconstruct bi-pred inter PU
 *
 * \param state          encoder state
 * \param ref1           reference picture to copy the data from
 * \param ref2           other reference picture to copy the data from
 * \param mv_param       motion vectors
 * \param lcu            destination lcu
 * \param predict_luma   Enable or disable luma prediction for this call.
 * \param predict_chroma Enable or disable chroma prediction for this call.
 * \param cu_loc         Size and position of current PU/CU
 */
void uvg_inter_recon_bipred(
  const encoder_state_t *const state,
  const uvg_picture *ref1,
  const uvg_picture *ref2,
  mv_t mv_param[2][2],
  lcu_t *lcu,
  bool predict_luma,
  bool predict_chroma,
  const cu_loc_t* const cu_loc)
{
  // Allocate maximum size arrays for interpolated and copied samples
  ALIGNED(64) uvg_pixel px_buf_L0[LCU_LUMA_SIZE + 2 * LCU_CHROMA_SIZE];
  ALIGNED(64) uvg_pixel px_buf_L1[LCU_LUMA_SIZE + 2 * LCU_CHROMA_SIZE];
  ALIGNED(64) uvg_pixel_im im_buf_L0[LCU_LUMA_SIZE + 2 * LCU_CHROMA_SIZE];
  ALIGNED(64) uvg_pixel_im im_buf_L1[LCU_LUMA_SIZE + 2 * LCU_CHROMA_SIZE];

  const int pu_x = cu_loc->x;
  const int pu_y = cu_loc->y;
  const int pu_w = cu_loc->width;
  const int pu_h = cu_loc->height;

  yuv_t px_L0;
  px_L0.size = pu_w * pu_h;
  px_L0.y = &px_buf_L0[0];
  px_L0.u = &px_buf_L0[LCU_LUMA_SIZE];
  px_L0.v = &px_buf_L0[LCU_LUMA_SIZE + LCU_CHROMA_SIZE];

  yuv_t px_L1;
  px_L1.size = pu_w * pu_h;
  px_L1.y = &px_buf_L1[0];
  px_L1.u = &px_buf_L1[LCU_LUMA_SIZE];
  px_L1.v = &px_buf_L1[LCU_LUMA_SIZE + LCU_CHROMA_SIZE];

  yuv_im_t im_L0;
  im_L0.size = pu_w * pu_h;
  im_L0.y = &im_buf_L0[0];
  im_L0.u = &im_buf_L0[LCU_LUMA_SIZE];
  im_L0.v = &im_buf_L0[LCU_LUMA_SIZE + LCU_CHROMA_SIZE];

  yuv_im_t im_L1;
  im_L1.size = pu_w * pu_h;
  im_L1.y = &im_buf_L1[0];
  im_L1.u = &im_buf_L1[LCU_LUMA_SIZE];
  im_L1.v = &im_buf_L1[LCU_LUMA_SIZE + LCU_CHROMA_SIZE];

  // Sample blocks from both reference picture lists.
  // Flags state if the outputs were written to high-precision / interpolated sample buffers.
  unsigned im_flags_L0 = inter_recon_unipred(state, ref1, pu_w, mv_param[0], &px_L0, &im_L0, predict_luma, predict_chroma,
                                             cu_loc);
  unsigned im_flags_L1 = inter_recon_unipred(state, ref2, pu_w, mv_param[1], &px_L1, &im_L1, predict_luma, predict_chroma,
                                             cu_loc);

  // After reconstruction, merge the predictors by taking an average of each pixel
  uvg_bipred_average(lcu, &px_L0, &px_L1, &im_L0, &im_L1,
                     pu_x, pu_y, pu_w, pu_h,
                     im_flags_L0, im_flags_L1,
                     predict_luma, predict_chroma);
}


/**
 * Reconstruct a single CU.
 *
 * The CU may consist of multiple PUs, each of which can use either
 * uniprediction or biprediction.
 *
 * \param state   encoder state
 * \param lcu     containing LCU
 * \param predict_luma   Enable or disable luma prediction for this call.
 * \param predict_chroma Enable or disable chroma prediction for this call.
 * \param cu_loc         Size and position of current CU
 */
void uvg_inter_recon_cu(
  const encoder_state_t * const state,
  lcu_t *lcu,
  bool predict_luma,
  bool predict_chroma,
  const cu_loc_t* const cu_loc)
{
  uvg_inter_pred_pu(state, lcu, predict_luma, predict_chroma, cu_loc);  
}

static void ibc_recon_cu(const encoder_state_t * const state,
                         lcu_t *lcu,
                         int32_t x,
                         int32_t y,
                         int32_t width,
                         bool predict_luma,
                         bool predict_chroma)
{
  const int x_scu    = SUB_SCU(x);
  const int y_scu    = SUB_SCU(y);
  uint32_t offset = x_scu + y_scu * LCU_WIDTH;
  uint32_t offset_c = x_scu / 2 + y_scu / 2 * LCU_WIDTH_C;
  cu_info_t *cu       = LCU_GET_CU_AT_PX(lcu, x_scu, y_scu);

  int32_t    mv_x     = cu->inter.mv[0][0] >> INTERNAL_MV_PREC;
  int32_t    mv_y     = cu->inter.mv[0][1] >> INTERNAL_MV_PREC;
  uint32_t   ibc_row  = y / LCU_WIDTH;

  int32_t    buffer_x = ((x - x_scu) + LCU_WIDTH <= IBC_BUFFER_WIDTH ?
                          x :
                          x - (((x - x_scu)) - IBC_BUFFER_WIDTH)) + mv_x;
  int32_t buffer_y = y_scu + mv_y;

  // The whole block must be to the left of the current position
  assert((-mv_x >= width || -mv_y >= width) && x >= 0 && y >= 0);

  // Predicted block completely outside of this LCU
  if (mv_x + x_scu + width <= 0) {  
    if(predict_luma) uvg_pixels_blit(&state->tile->frame->ibc_buffer_y[ibc_row][buffer_y * IBC_BUFFER_WIDTH + buffer_x], lcu->rec.y + offset, width, width, IBC_BUFFER_WIDTH, LCU_WIDTH);
    if (predict_chroma) {
      uvg_pixels_blit(&state->tile->frame->ibc_buffer_u[ibc_row][(buffer_y / 2) * IBC_BUFFER_WIDTH_C + (buffer_x / 2)], lcu->rec.u + offset_c, width / 2, width / 2, IBC_BUFFER_WIDTH_C, LCU_WIDTH_C);
      uvg_pixels_blit(&state->tile->frame->ibc_buffer_v[ibc_row][(buffer_y / 2) * IBC_BUFFER_WIDTH_C + (buffer_x / 2)], lcu->rec.v + offset_c, width / 2, width / 2, IBC_BUFFER_WIDTH_C, LCU_WIDTH_C);
    }
  } else if (mv_x + x_scu + width >= width) { // Completely in current LCU
    if(predict_luma) uvg_pixels_blit(&lcu->rec.y[(y_scu + mv_y) * LCU_WIDTH + x_scu + mv_x], lcu->rec.y + offset, width, width, LCU_WIDTH, LCU_WIDTH);
    if (predict_chroma) {
      uvg_pixels_blit(&lcu->rec.u[((y_scu+mv_y) / 2) * LCU_WIDTH_C + (x_scu + mv_x) / 2], lcu->rec.u + offset_c, width / 2, width / 2, LCU_WIDTH_C, LCU_WIDTH_C);
      uvg_pixels_blit(&lcu->rec.v[((y_scu+mv_y) / 2) * LCU_WIDTH_C + (x_scu + mv_x) / 2], lcu->rec.v + offset_c, width / 2, width / 2, LCU_WIDTH_C, LCU_WIDTH_C);
    }
  } else { // Partly on the buffer and party on the current LCU rec

    uint32_t width_buffer = -(mv_x + x_scu);
    uint32_t width_lcu    = width - width_buffer;
    if(predict_luma) uvg_pixels_blit(&state->tile->frame->ibc_buffer_y[ibc_row][buffer_y * IBC_BUFFER_WIDTH + buffer_x], lcu->rec.y + offset, width_buffer, width, IBC_BUFFER_WIDTH, LCU_WIDTH);
    if (predict_chroma) {    
      uvg_pixels_blit(&state->tile->frame->ibc_buffer_u[ibc_row][(buffer_y / 2) * IBC_BUFFER_WIDTH_C + (buffer_x / 2)], lcu->rec.u + offset_c, width_buffer / 2 + (width_buffer&1), width / 2, IBC_BUFFER_WIDTH_C, LCU_WIDTH_C);
      uvg_pixels_blit(&state->tile->frame->ibc_buffer_v[ibc_row][(buffer_y / 2) * IBC_BUFFER_WIDTH_C + (buffer_x / 2)], lcu->rec.v + offset_c, width_buffer / 2 + (width_buffer&1), width / 2, IBC_BUFFER_WIDTH_C, LCU_WIDTH_C);
    }

    offset += width_buffer;
    offset_c += width_buffer/2 + (width_buffer&1);

    if(predict_luma) uvg_pixels_blit(&lcu->rec.y[(y_scu + mv_y) * LCU_WIDTH + x_scu + mv_x + width_buffer], lcu->rec.y + offset, width_lcu, width, LCU_WIDTH, LCU_WIDTH);
    if (predict_chroma && (width_lcu / 2)) {
      uvg_pixels_blit(&lcu->rec.u[((y_scu+mv_y) / 2) * LCU_WIDTH_C + (x_scu + mv_x + width_buffer) / 2], lcu->rec.u + offset_c, width_lcu / 2, width / 2, LCU_WIDTH_C, LCU_WIDTH_C);
      uvg_pixels_blit(&lcu->rec.v[((y_scu+mv_y) / 2) * LCU_WIDTH_C + (x_scu + mv_x + width_buffer) / 2], lcu->rec.v + offset_c, width_lcu / 2, width / 2, LCU_WIDTH_C, LCU_WIDTH_C);
    }
  }
}

/**
 * Predict a single PU.
 *
 * The PU may use either uniprediction or biprediction.
 *
 * \param state          encoder state
 * \param lcu            containing LCU
 * \param predict_luma   Enable or disable luma prediction for this call.
 * \param predict_chroma Enable or disable chroma prediction for this call.
 * \param cu_loc         Size and position of current PU/CU
 */
void uvg_inter_pred_pu(
  const encoder_state_t * const state,
  lcu_t *lcu,
  bool predict_luma,
  bool predict_chroma,
  const cu_loc_t* const cu_loc)

{
  const int x_scu = SUB_SCU(cu_loc->x);
  const int y_scu = SUB_SCU(cu_loc->y);
  cu_info_t *pu = LCU_GET_CU_AT_PX(lcu, x_scu, y_scu);

  if (pu->inter.mv_dir == 3) {
    const uvg_picture *const refs[2] = {
      state->frame->ref->images[
        state->frame->ref_LX[0][
          pu->inter.mv_ref[0]]],
      state->frame->ref->images[
        state->frame->ref_LX[1][
          pu->inter.mv_ref[1]]],
    };
    uvg_inter_recon_bipred(state,
                           refs[0], refs[1],
                           pu->inter.mv, lcu,
                           predict_luma, predict_chroma,
                           cu_loc);
  }
  else if (pu->type == CU_IBC) {
    ibc_recon_cu(state, lcu, cu_loc->x, cu_loc->y, cu_loc->width, predict_luma, predict_chroma);
  } else{
    const int mv_idx = pu->inter.mv_dir - 1;
    const uvg_picture *const ref =
      state->frame->ref->images[
        state->frame->ref_LX[mv_idx][
          pu->inter.mv_ref[mv_idx]]];

    const unsigned offset_luma = SUB_SCU(cu_loc->y) * LCU_WIDTH + SUB_SCU(cu_loc->x);
    const unsigned offset_chroma = SUB_SCU(cu_loc->y) / 2 * LCU_WIDTH_C + SUB_SCU(cu_loc->x) / 2;
    yuv_t lcu_adapter;
    lcu_adapter.size = cu_loc->width * cu_loc->height;
    lcu_adapter.y = lcu->rec.y + offset_luma;
    lcu_adapter.u = lcu->rec.u + offset_chroma;
    lcu_adapter.v = lcu->rec.v + offset_chroma;

    inter_recon_unipred(state,
                        ref,
                        LCU_WIDTH, pu->inter.mv[mv_idx],
                        &lcu_adapter,
                        NULL,
                        predict_luma,
                        predict_chroma,
                        cu_loc);
  }
  if (predict_chroma && state->encoder_control->cfg.jccr) {
    const int offset = x_scu / 2 + y_scu / 2 * LCU_WIDTH_C;
    uvg_pixels_blit(lcu->rec.u + offset, lcu->rec.joint_u + offset, cu_loc->chroma_width, cu_loc->chroma_height, LCU_WIDTH_C, LCU_WIDTH_C);
    uvg_pixels_blit(lcu->rec.v + offset, lcu->rec.joint_v + offset, cu_loc->chroma_width, cu_loc->chroma_height, LCU_WIDTH_C, LCU_WIDTH_C);
  }
}

/**
 * \brief Clear unused L0/L1 motion vectors and reference
 * \param cu coding unit to clear
 */
static void inter_clear_cu_unused(cu_info_t* cu)
{
  for (unsigned i = 0; i < 2; ++i) {
    if (cu->inter.mv_dir & (1 << i)) continue;

    cu->inter.mv[i][0] = 0;
    cu->inter.mv[i][1] = 0;
    cu->inter.mv_ref[i] = 255;
  }
}


static bool is_cand_coded(int cur_x, int cur_y, int cand_x, int cand_y, uint32_t split_tree) 
{
  //Start from the LCU and figure out which sub-blocks the cur and cand cu belong to
  int log2_width = LOG2_LCU_WIDTH; int log2_height = LOG2_LCU_WIDTH;
  int cur_block_x = (cur_x >> log2_width); int cur_block_y = (cur_y >> log2_height);
  int cand_block_x = (cand_x >> log2_width); int cand_block_y = (cand_y >> log2_height);
  
  // If different LCU and cand block is before cur in raster scan order, it is coded.
  if (cur_block_y != cand_block_y) {
    
    return cand_block_y < cur_block_y;
  }
  else if (cur_block_x != cand_block_x) {
    return cand_block_x < cur_block_x;
  }

  struct {
    uint32_t split_tree;
  } split_data = {split_tree};

  int offset_x = 0; int offset_y = 0;
  int cur_cu_index = 0;
  int cand_cu_index = 0;
  int depth = 0;
  do
  {
    uint32_t split = GET_SPLITDATA(&split_data, depth);
    depth++;
    // Figure out which sub-block cand and cur belong to in the current split
    switch (split)
    {
    case QT_SPLIT: // Four way split, index based on both x and y
      log2_width--; log2_height--;

      cur_block_x = ((cur_x + offset_x) >> log2_width) & 1;
      cur_block_y = ((cur_y + offset_y) >> log2_width) & 1;
      cur_cu_index = cur_block_x + 2 * cur_block_y;
      
      cand_block_x = ((cand_x + offset_x) >> log2_height) & 1;
      cand_block_y = ((cand_y + offset_y) >> log2_height) & 1;
      cand_cu_index = cand_block_x + 2 * cand_block_y;
      break;

    case BT_HOR_SPLIT: // Two way split, index only based on y
      log2_height--;

      cur_block_y = ((cur_y + offset_y) >> log2_height) & 1;
      cur_cu_index = cur_block_y;

      cand_block_y = ((cand_y + offset_y) >> log2_height) & 1;
      cand_cu_index = cand_block_y;
      break;

    case BT_VER_SPLIT: // Two way split, index only based on x
      log2_width--;

      cur_block_x = ((cur_x + offset_x) >> log2_width) & 1;
      cur_cu_index = cur_block_x;

      cand_block_x = ((cand_x + offset_x) >> log2_width) & 1;
      cand_cu_index = cand_block_x;
      break;

    case TT_HOR_SPLIT: // Three way split, index only based on y. Need log2_height + 1 and log2_height bit to determine index. block == 0 -> index := 0, block == {1,2} -> index := 1, block == 3 -> index := 2 
      log2_height -= 2; //set to smallest block size

      cur_block_y = ((cur_y + offset_y) >> log2_height) & 3;
      cur_cu_index = (cur_block_y == 0) ? 0 : ((cur_block_y != 3) ? 1 : 2);

      cand_block_y = ((cand_y + offset_y) >> log2_height) & 3;
      cand_cu_index = (cand_block_y == 0) ? 0 : ((cand_block_y != 3) ? 1 : 2);

      if (cur_cu_index == 1) {
        // TT split causes misalignment in the middle block, so need to use offset to get correct alignment (to a (1 << log2_size)-grid) for indexing.
        // Alignment is fixed after a BT split, but no need to reset offset since it does not affect indexing for smaller blocks than current (1 << log2_size). 
        offset_y = (1 << log2_height); 
        log2_height++; //middle block is larger so need to increase size
      }
      break;

    case TT_VER_SPLIT: // Three way split, index only based on x. Need log2_width + 1 and log2_width bit to determine index. block == 0 -> index := 0, block == {1,2} -> index := 1, block == 3 -> index := 2 
      log2_width -= 2; //set to smallest block size

      cur_block_x = ((cur_x + offset_x) >> log2_width) & 3;
      cur_cu_index = (cur_block_x == 0) ? 0 : ((cur_block_x != 3) ? 1 : 2);

      cand_block_x = ((cand_x + offset_x) >> log2_width) & 3;
      cand_cu_index = (cand_block_x == 0) ? 0 : ((cand_block_x != 3) ? 1 : 2);
      
      if (cur_cu_index == 1) {
        // TT split causes misalignment in the middle block, so need to use offset to get correct alignment (to a (1 << log2_size)-grid) for indexing.
        // Alignment is fixed after a BT split, but no need to reset offset since it does not affect indexing for smaller blocks than current (1 << log2_size). 
        offset_x = (1 << log2_width); 
        log2_width++; //middle block is larger so need to increase size
      }
      break;

    default:
      assert(false && "Not a valid split");
    }
    if (cand_cu_index != cur_cu_index) return cand_cu_index < cur_cu_index;

  } while (depth < MAX_DEPTH);

  assert(false && "Either max depth reached or cur and cand are the same block"); // We should not get here 
  return false;
}

/**
 * \brief Check whether a0 mv cand block is coded before the current block.
 * \param x       x-coordinate of the current block (in pixels)
 * \param y       y-coordinate of the current block (in pixels)
 * \param width   width of the current block (in pixels)
 * \param height  height of the current block (in pixels)
 * \return        True, if the a0 mv candidate block is coded before the
 *                current block. Otherwise false.
 */
static bool is_a0_cand_coded(int x, int y, int width, int height)
{
  int size = MIN(width & ~(width - 1), height & ~(height - 1));

  if (height != size) {
    // For SMP and AMP blocks the situation is equivalent to a square block
    // at the lower left corner of the PU.
    y = y + height - size;
  }

  while (size < LCU_WIDTH) {
    const int parent_size = 2 * size;
    const int cu_index    = (x % parent_size != 0) + 2 * (y % parent_size != 0);
    switch (cu_index) {
      case 0:
        // A0 is in the CU directly left of the parent CU so it has been
        // coded already.
        //    +---+---+
        //    | X |   |
        //    |---+---+
        // A0 |   |   |
        //    +---+---+
        return true;

      case 1:
        // A0 is in the CU that will be coded after the current CU.
        //    +---+---+
        //    |   | X |
        //    |---+---+
        //    |A0 |   |
        //    +---+---+
        return false;

      case 2:
        //    +---+---+
        //    |   |   |
        //    |---+---+
        //    | X |   |
        //    +---+---+
        // A0

        // Move to the parent block.
        y -= size;
        size = parent_size;
        break;

      case 3:
        // A0 is in the CU directly down of the parent CU so is has not
        // been coded yet.
        //    +---+---+
        //    |   |   |
        //    |---+---+
        //    |   | X |
        //    +---+---+
        //     A0
        return false;
    }
  }

  // For 64x64 blocks A0 candidate is located outside the LCU.
  return false;
}

/**
 * \brief Check whether b0 mv cand block is coded before the current block.
 * \param x       x-coordinate of the current block (in pixels)
 * \param y       y-coordinate of the current block (in pixels)
 * \param width   width of the current block (in pixels)
 * \param height  height of the current block (in pixels)
 * \return        True, if the b0 mv candidate block is coded before the
 *                current block. Otherwise false.
 */
static bool is_b0_cand_coded(int x, int y, int width, int height)
{
  int size = MIN(width & ~(width - 1), height & ~(height - 1));

  if (width != size) {
    // For SMP and AMP blocks the situation is equivalent to a square block
    // at the upper right corner of the PU.
    x = x + width - size;
  }

  while (size < LCU_WIDTH) {
    const int parent_size = 2 * size;
    const int cu_index    = (x % parent_size != 0) + 2 * (y % parent_size != 0);
    switch (cu_index) {
      case 0:
        // B0 is in the CU directly above the parent CU so it has been
        // coded already.
        //         B0
        //    +---+---+
        //    | X |   |
        //    |---+---+
        //    |   |   |
        //    +---+---+
        return true;

      case 1:
        //             B0
        //    +---+---+
        //    |   | X |
        //    |---+---+
        //    |   |   |
        //    +---+---+

        // Move to the parent block.
        x -= size;
        size = parent_size;
        break;

      case 2:
        //    +---+---+
        //    |   |B0 |
        //    |---+---+
        //    | X |   |
        //    +---+---+
        return true;

      case 3:
        // B0 is in the CU directly right of the parent CU so is has not
        // been coded yet.
        //    +---+---+
        //    |   |   | B0
        //    |---+---+
        //    |   | X |
        //    +---+---+
        return false;
    }
  }

  // The LCU to the right and up of the current LCU has been coded already.
  return true;
}


/**
 * \brief Get merge candidates for current block
 *
 * \param state     encoder control state to use
 * \param cu_loc    Size and position of current CU
 * \param ref_list  which reference list, L0 is 1 and L1 is 2
 * \param ref_idx   index in the reference list
 * \param cand_out  will be filled with C0 and C1 candidates
 */
static void get_temporal_merge_candidates(
  const encoder_state_t * const state,
  const cu_loc_t* const cu_loc,
  uint8_t ref_list,
  uint8_t ref_idx,
  merge_candidates_t *cand_out)
{
  /*
  Predictor block locations
  _________
  |CurrentPU|
  |     __  |
  |    |C1| |
  |_________|_
            |C0|
  */

  cand_out->c0 = cand_out->c1 = NULL;

  // Find temporal reference
  if (state->frame->ref->used_size) {
    uint32_t colocated_ref;

    // Select L0/L1 ref_idx reference
    if (state->frame->ref_LX_size[ref_list-1] > ref_idx) {
      colocated_ref = state->frame->ref_LX[ref_list - 1][ref_idx];
    } else {
      // not found
      return;
    }

    cu_array_t *ref_cu_array = state->frame->ref->cu_arrays[colocated_ref];
    int cu_per_width = ref_cu_array->width / SCU_WIDTH;

    int32_t xColBr = cu_loc->x + cu_loc->width;
    int32_t yColBr = cu_loc->y + cu_loc->height;

    // C0 must be available
    if (xColBr < state->encoder_control->in.width &&
        yColBr < state->encoder_control->in.height) {
      int32_t C0_offset = -1;

      // Y inside the current CTU / LCU
      if (yColBr % LCU_WIDTH != 0) {
        C0_offset = ((xColBr >> 3) << 3) / SCU_WIDTH +
                  (((yColBr >> 3) << 3) / SCU_WIDTH) * cu_per_width;
      }

      if (C0_offset >= 0) {
        // Only use when it's inter block
        if (ref_cu_array->data[C0_offset].type == CU_INTER) {
          cand_out->c0 = &ref_cu_array->data[C0_offset];
        }
      }
    }
    int32_t xColCtr = cu_loc->x + (cu_loc->width / 2);
    int32_t yColCtr = cu_loc->y + (cu_loc->height / 2);

    // C1 must be inside the LCU, in the center position of current CU
    if (xColCtr < state->encoder_control->in.width && yColCtr < state->encoder_control->in.height) {
      uint32_t C1_offset = ((xColCtr >> 3) << 3) / SCU_WIDTH + ((((yColCtr >> 3) << 3) / SCU_WIDTH) * cu_per_width);
      if (ref_cu_array->data[C1_offset].type == CU_INTER) {
        cand_out->c1 = &ref_cu_array->data[C1_offset];
      }
    }
  }
}

static INLINE mv_t get_scaled_mv(mv_t mv, int scale)
{
  int32_t scaled = scale * mv;
  return CLIP(-131072, 131071, (scaled + 127 + (scaled < 0)) >> 8);
}

#define MV_EXPONENT_BITCOUNT 4
#define MV_MANTISSA_BITCOUNT 6
#define MV_MANTISSA_UPPER_LIMIT ((1 << (MV_MANTISSA_BITCOUNT - 1)) - 1)
#define MV_MANTISSA_LIMIT (1 << (MV_MANTISSA_BITCOUNT - 1))
#define MV_EXPONENT_MASK ((1 << MV_EXPONENT_BITCOUNT) - 1)

static int convert_mv_fixed_to_float(int32_t val)
{
  uint32_t sign = val >> 31;
  int scale = uvg_math_floor_log2((val ^ sign) | MV_MANTISSA_UPPER_LIMIT) - (MV_MANTISSA_BITCOUNT - 1);

  int exponent;
  uint32_t mantissa;
  if (scale >= 0)
  {
    int round = (1 << scale) >> 1;
    int n = (val + round) >> scale;
    exponent = scale + ((n ^ sign) >> (MV_MANTISSA_BITCOUNT - 1));
    mantissa = (n & MV_MANTISSA_UPPER_LIMIT) | (sign << (MV_MANTISSA_BITCOUNT - 1));
  }
  else
  {
    exponent = 0;
    mantissa = val;
  }

  return exponent | (mantissa << MV_EXPONENT_BITCOUNT);
}

static int convert_mv_float_to_fixed(int val)
{
  int exponent = val & MV_EXPONENT_MASK;
  uint32_t mantissa = val >> MV_EXPONENT_BITCOUNT;
  return exponent == 0 ? mantissa : (mantissa ^ MV_MANTISSA_LIMIT) << (exponent - 1);
}

static int round_mv_comp(int x)
{
  return convert_mv_float_to_fixed(convert_mv_fixed_to_float(x));
}

static void apply_mv_scaling_pocs(int32_t current_poc,
                                  int32_t current_ref_poc,
                                  int32_t neighbor_poc,
                                  int32_t neighbor_ref_poc,
                                  mv_t mv_cand[2])
{
  int32_t diff_current  = current_poc  - current_ref_poc;
  int32_t diff_neighbor = neighbor_poc - neighbor_ref_poc;

  if (diff_current == diff_neighbor) return;

  diff_current  = CLIP(-128, 127, diff_current);
  diff_neighbor = CLIP(-128, 127, diff_neighbor);

  int scale = CLIP(-4096, 4095,
    (diff_current * ((0x4000 + (abs(diff_neighbor) >> 1)) / diff_neighbor) + 32) >> 6);

  mv_cand[0] = get_scaled_mv(mv_cand[0], scale);
  mv_cand[1] = get_scaled_mv(mv_cand[1], scale);
}

static INLINE void apply_mv_scaling(const encoder_state_t *state,
                                    const cu_info_t *current_cu,
                                    const cu_info_t *neighbor_cu,
                                    int8_t current_reflist,
                                    int8_t neighbor_reflist,
                                    mv_t mv_cand[2])
{
  apply_mv_scaling_pocs(state->frame->poc,
                        state->frame->ref->pocs[
                          state->frame->ref_LX[current_reflist][
                            current_cu->inter.mv_ref[current_reflist]]],
                        state->frame->poc,
                        state->frame->ref->pocs[
                          state->frame->ref_LX[neighbor_reflist][
                            neighbor_cu->inter.mv_ref[neighbor_reflist]]],
                        mv_cand);
}

static INLINE bool add_mvp_candidate(const encoder_state_t *state,
                                     const cu_info_t *cur_cu,
                                     const cu_info_t *cand,
                                     int8_t reflist,
                                     bool scaling,
                                     mv_t mv_cand_out[2])
{
  if (!cand) return false;

  assert(cand->inter.mv_dir != 0);

  for (int i = 0; i < 2; i++) {
    const int cand_list = i == 0 ? reflist : !reflist;

    if ((cand->inter.mv_dir & (1 << cand_list)) == 0) continue;

    if (scaling) {
      mv_cand_out[0] = cand->inter.mv[cand_list][0];
      mv_cand_out[1] = cand->inter.mv[cand_list][1];
      apply_mv_scaling(state, cur_cu, cand, reflist, cand_list, mv_cand_out);
      return true;
    }

    if (state->frame->ref_LX[cand_list][cand->inter.mv_ref[cand_list]] ==
        state->frame->ref_LX[reflist][cur_cu->inter.mv_ref[reflist]])
    {
      mv_cand_out[0] = cand->inter.mv[cand_list][0];
      mv_cand_out[1] = cand->inter.mv[cand_list][1];
      return true;
    }
  }

  return false;
}


static bool is_duplicate_candidate_ibc(const cu_info_t* cu1, const cu_info_t* cu2)
{
  if (!cu2) return false;

  if (cu1->inter.mv[0][0] != cu2->inter.mv[0][0]  ||
      cu1->inter.mv[0][1] != cu2->inter.mv[0][1]) {
    return false;
  }


  return true;
}

/**
 * \brief Get merge candidates for current block.
 *
 * The output parameters b0, b1, b2, a0, a1 are pointed to the
 * corresponding cu_info_t struct in lcu->cu, or set to NULL, if the
 * candidate is not available.
 *
 * \param x               block x position in pixels
 * \param y               block y position in pixels
 * \param width           block width in pixels
 * \param height          block height in pixels
 * \param picture_width   tile width in pixels
 * \param picture_height  tile height in pixels
 * \param lcu             current LCU
 * \param cand_out        will be filled with A and B candidates
 */
static void get_ibc_merge_candidates(const encoder_state_t * const state,
                                     const cu_info_t * const cur_cu,
                                     lcu_t *lcu,
                                     const cu_array_t *cua,
                                     int32_t x,
                                     int32_t y,
                                     int32_t width,
                                     int32_t height,
                                     mv_t mv_cand[IBC_MRG_MAX_NUM_CANDS][2]
  )
{
  /*
  Predictor block locations
  ____      _______
  |B2|______|B1|B0|
     |         |
     |  Cur CU |
   __|         |
  |A1|_________|
  |A0|
  */
  int32_t x_local = SUB_SCU(x); //!< coordinates from top-left of this LCU
  int32_t y_local = SUB_SCU(y);

  cu_info_t *a1      = NULL;
  cu_info_t *b1      = NULL;
    
  uint8_t candidates = 0;

  // A1 availability testing
  if (x != 0) {
    a1 = lcu != NULL?LCU_GET_CU_AT_PX(lcu, x_local - 1, y_local + height - 1): uvg_cu_array_at_const(cua, x - 1, y + height - 1);
    // Do not check a1->coded because the block above is always coded before
    // the current one and the flag is not set when searching an SMP block.
    if (a1->type == CU_IBC) {
      inter_clear_cu_unused(a1);
      mv_cand[candidates][0] = a1->inter.mv[0][0];
      mv_cand[candidates][1] = a1->inter.mv[0][1];
      candidates++;
    } else {
      a1 = NULL;
    }
  }

  // B1 availability testing
  if (y != 0) {
    b1 =  lcu != NULL?LCU_GET_CU_AT_PX(lcu, x_local + width - 1, y_local - 1): uvg_cu_array_at_const(cua, x + width - 1, y - 1);
    // Do not check b1->coded because the block to the left is always coded
    // before the current one and the flag is not set when searching an SMP
    // block.
    if (b1->type == CU_IBC) {
      if(!is_duplicate_candidate_ibc(b1, a1)) {
        inter_clear_cu_unused(b1);
        mv_cand[candidates][0] = b1->inter.mv[0][0];
        mv_cand[candidates][1] = b1->inter.mv[0][1];
        candidates++;
      }
    } else {
      b1 = NULL;
    }
  }  

  if (candidates > 0)
    uvg_round_precision(INTERNAL_MV_PREC, 2, &mv_cand[0][0], &mv_cand[0][1]);
  if (candidates > 1)
    uvg_round_precision(INTERNAL_MV_PREC, 2, &mv_cand[1][0], &mv_cand[1][1]);

  if (candidates < IBC_MRG_MAX_NUM_CANDS)
  {
    const uint32_t ctu_row = (y >> LOG2_LCU_WIDTH);
    const uint32_t ctu_row_mul_five = ctu_row * MAX_NUM_HMVP_CANDS;
    int32_t num_cand = state->tile->frame->hmvp_size_ibc[ctu_row];
    for (int i = 0; i < MIN(MAX_NUM_HMVP_CANDS,num_cand); i++) {
      cu_info_t* cand = &state->tile->frame->hmvp_lut_ibc[ctu_row_mul_five + i];
      bool       duplicate = false;

      // Check that the HMVP candidate is not duplicate
      if (is_duplicate_candidate_ibc(cand, a1)) {
        duplicate = true;
      } else if(is_duplicate_candidate_ibc(cand, b1)) {
        duplicate = true;
      }

      // allow duplicates after the first hmvp lut item
      if (!duplicate || i > 0) {
        mv_cand[candidates][0] = cand->inter.mv[0][0];
        mv_cand[candidates][1] = cand->inter.mv[0][1];
        candidates++;      
        if (candidates == IBC_MRG_MAX_NUM_CANDS) return;
      }
    }
  }

  // Fill with (0,0)
  while (candidates < IBC_MRG_MAX_NUM_CANDS) {
    mv_cand[candidates][0] = 0;
    mv_cand[candidates][1] = 0;
    candidates++;
  }
}


/**
 * \brief Get merge candidates for current block.
 *
 * The output parameters b0, b1, b2, a0, a1 are pointed to the
 * corresponding cu_info_t struct in lcu->cu, or set to NULL, if the
 * candidate is not available.
 *
 * \param x               block x position in pixels
 * \param y               block y position in pixels
 * \param width           block width in pixels
 * \param height          block height in pixels
 * \param picture_width   tile width in pixels
 * \param picture_height  tile height in pixels
 * \param lcu             current LCU
 * \param cand_out        will be filled with A and B candidates
 */
static void get_spatial_merge_candidates(const cu_loc_t* const cu_loc,
                                         int32_t picture_width,
                                         int32_t picture_height,
                                         lcu_t *lcu,
                                         merge_candidates_t *cand_out,
                                         uint8_t parallel_merge_level,
                                         bool wpp
  )
{
  /*
  Predictor block locations
  ____      _______
  |B2|______|B1|B0|
     |         |
     |  Cur CU |
   __|         |
  |A1|_________|
  |A0|
  */
  const int32_t x_local = SUB_SCU(cu_loc->x); //!< coordinates from top-left of this LCU
  const int32_t y_local = SUB_SCU(cu_loc->y);

  const int x = cu_loc->x;
  const int y = cu_loc->y;
  const int width = cu_loc->width;
  const int height = cu_loc->height;
  const cu_info_t* const cur_cu = LCU_GET_CU_AT_PX(lcu, x_local, y_local);

  // A0 and A1 availability testing
  if (x != 0) {
    cu_info_t *a1 = LCU_GET_CU_AT_PX(lcu, x_local - 1, y_local + height - 1);
    // Do not check a1->coded because the block above is always coded before
    // the current one and the flag is not set when searching an SMP block.
    if (a1->type == CU_INTER) {
      inter_clear_cu_unused(a1);
      cand_out->a[1] = a1;
    }

    if (y_local + height < LCU_WIDTH && y + height < picture_height) {
      cu_info_t *a0 = LCU_GET_CU_AT_PX(lcu, x_local - 1, y_local + height);
      if (a0->type == CU_INTER && is_cand_coded(x, y, x - 1, y + height, cur_cu->split_tree)) {
        inter_clear_cu_unused(a0);
        cand_out->a[0] = a0;
      }
    }
  }

  // B0, B1 and B2 availability testing
  if (y != 0) {
    cu_info_t *b0 = NULL;
    if (x + width < picture_width) {
      if (x_local + width < LCU_WIDTH) {
        b0 = LCU_GET_CU_AT_PX(lcu, x_local + width, y_local - 1);
      } else if (!wpp && y_local == 0) {
        // Special case, top-right CU
        b0 = LCU_GET_TOP_RIGHT_CU(lcu);
      }
    }
    if (b0 && b0->type == CU_INTER && is_cand_coded(x, y, x + width, y - 1, cur_cu->split_tree)) {
      inter_clear_cu_unused(b0);
      cand_out->b[0] = b0;
    }

    cu_info_t *b1 = LCU_GET_CU_AT_PX(lcu, x_local + width - 1, y_local - 1);
    // Do not check b1->coded because the block to the left is always coded
    // before the current one and the flag is not set when searching an SMP
    // block.
    if (b1->type == CU_INTER) {
      inter_clear_cu_unused(b1);
      cand_out->b[1] = b1;
    }

    if (x != 0) {
      cu_info_t *b2 = LCU_GET_CU_AT_PX(lcu, x_local - 1, y_local - 1);
      // Do not check b2->coded because the block above and to the left is
      // always coded before the current one.
      if (b2->type == CU_INTER) {
        inter_clear_cu_unused(b2);
        cand_out->b[2] = b2;
      }
    }
  }
}

/**
 * \brief Get merge candidates for current block.
 *
 * The output parameters b0, b1, b2, a0, a1 are pointed to the
 * corresponding cu_info_t struct in lcu->cu, or set to NULL, if the
 * candidate is not available.
 *
 * \param cua             cu information
 * \param x               block x position in pixels
 * \param y               block y position in pixels
 * \param width           block width in pixels
 * \param height          block height in pixels
 * \param picture_width   tile width in pixels
 * \param picture_height  tile height in pixels
 * \param cand_out        will be filled with A and B candidates
 */
static void get_spatial_merge_candidates_cua(
  const cu_array_t *cua,
  int32_t picture_width,
  int32_t picture_height,
  merge_candidates_t *cand_out,
  bool wpp,
  const cu_loc_t* const cu_loc)
{
  /*
  Predictor block locations
  ____      _______
  |B2|______|B1|B0|
     |         |
     |  Cur CU |
   __|         |
  |A1|_________|
  |A0|
  */
  const int x = cu_loc->x;
  const int y = cu_loc->y;
  const int width = cu_loc->width;
  const int height = cu_loc->height;
  const int32_t x_local = SUB_SCU(x); //!< coordinates from top-left of this LCU
  const int32_t y_local = SUB_SCU(y);
  const cu_info_t* const cur_cu = uvg_cu_array_at_const(cua, x, y);

  // A0 and A1 availability testing
  if (x != 0) {
    const cu_info_t *a1 = uvg_cu_array_at_const(cua, x - 1, y + height - 1);
    // The block above is always coded before the current one.
    if (a1->type == CU_INTER) {
      cand_out->a[1] = a1;
    }

    if (y_local + height < LCU_WIDTH && y + height < picture_height) {
      const cu_info_t *a0 = uvg_cu_array_at_const(cua, x - 1, y + height);
      if (a0->type == CU_INTER && is_cand_coded(x, y, x - 1, y + height, cur_cu->split_tree)) {
        cand_out->a[0] = a0;
      }
    }
  }

  // B0, B1 and B2 availability testing
  if (y != 0) {
    if (x + width < picture_width && (x_local + width < LCU_WIDTH || (!wpp && y_local == 0))) {
      const cu_info_t *b0 = uvg_cu_array_at_const(cua, x + width, y - 1);
      if (b0->type == CU_INTER && is_cand_coded(x, y, x + width, y - 1, cur_cu->split_tree)) {
        cand_out->b[0] = b0;
      }
    }

    const cu_info_t* b1 = uvg_cu_array_at_const(cua, x + width - 1, y - 1);
    // The block to the left is always coded before the current one.
    if (b1->type == CU_INTER) {
      cand_out->b[1] = b1;
    }

    if (x != 0) {
      const cu_info_t *b2 = uvg_cu_array_at_const(cua, x - 1, y - 1);
      // The block above and to the left is always coded before the current
      // one.
      if (b2->type == CU_INTER) {
        cand_out->b[2] = b2;
      }
    }
  }
}

/**
 * \brief Try to add a temporal MVP or merge candidate.
 *
 * \param state         encoder state
 * \param current_ref   index of the picture referenced by the current CU
 * \param colocated     colocated CU
 * \param reflist       either 0 (for L0) or 1 (for L1)
 * \param[out] mv_out   Returns the motion vector
 *
 * \return Whether a temporal candidate was added or not.
 */
static bool add_temporal_candidate(const encoder_state_t *state,
                                   uint8_t current_ref,
                                   const cu_info_t *colocated,
                                   int32_t reflist,
                                   mv_t mv_out[2])
{
  if (!colocated) return false;

  int colocated_ref;
  if (state->frame->ref_LX_size[0] > 0) {
    // get the first reference from L0 if it exists
    colocated_ref = state->frame->ref_LX[0][0];
  } else {
    // otherwise no candidate added
    return false;
  }

  // When there are reference pictures from the future (POC > current POC)
  // in L0 or L1, the primary list for the colocated PU is the inverse of
  // collocated_from_l0_flag. Otherwise it is equal to reflist.
  //
  // uvg266 always sets collocated_from_l0_flag so the list is L1 when
  // there are future references.
  int col_list = reflist;
  for (uint32_t i = 0; i < state->frame->ref->used_size; i++) {
    if (state->frame->ref->pocs[i] > state->frame->poc) {
      col_list = 1;
      break;
    }
  }
  
  if ((colocated->inter.mv_dir & (col_list + 1)) == 0) {
    // Use the other list if the colocated PU does not have a MV for the
    // primary list.
    col_list = 1 - col_list;
  }

  mv_out[0] = colocated->inter.mv[col_list][0];
  mv_out[1] = colocated->inter.mv[col_list][1];

  mv_out[0] = round_mv_comp(mv_out[0]);
  mv_out[1] = round_mv_comp(mv_out[1]);

  apply_mv_scaling_pocs(
    state->frame->poc,
    state->frame->ref->pocs[current_ref],
    state->frame->ref->pocs[colocated_ref],
    state->frame->ref->images[colocated_ref]->ref_pocs[
      state->frame->ref->ref_LXs[colocated_ref]
        [col_list][colocated->inter.mv_ref[col_list]]],
    mv_out
  );

  return true;
}

/**
 * \brief Pick two mv candidates from the spatial and temporal candidates.
 */
static void get_mv_cand_from_candidates(
  const encoder_state_t * const state,
  const merge_candidates_t *merge_cand,
  const cu_info_t * const cur_cu,
  int8_t reflist,
  mv_t mv_cand[2][2],
  int ctu_row)
{
  const cu_info_t *const *a = merge_cand->a;
  const cu_info_t *const *b = merge_cand->b;
  const cu_info_t *c0 = merge_cand->c0;
  const cu_info_t *c1  = merge_cand->c1;

  uint8_t candidates = 0;
  uint8_t b_candidates = 0;

  // Left predictors without scaling
  if (add_mvp_candidate(state, cur_cu, a[0], reflist, false, mv_cand[candidates])) {
    candidates++;
  } else if (add_mvp_candidate(state, cur_cu, a[1], reflist, false, mv_cand[candidates])) {
    candidates++;
  }


  // Top predictors without scaling  
  if (add_mvp_candidate(state, cur_cu, b[0], reflist, false, mv_cand[candidates])) {
    b_candidates++;
  } else if (add_mvp_candidate(state, cur_cu, b[1], reflist, false, mv_cand[candidates])) {
    b_candidates++;
  }
  else if (add_mvp_candidate(state, cur_cu, b[2], reflist, false, mv_cand[candidates])) {
    b_candidates++;
  }

  candidates += b_candidates;

  if (candidates > 0)
    uvg_round_precision(INTERNAL_MV_PREC, 2, &mv_cand[0][0], &mv_cand[0][1]);
  if (candidates > 1)
    uvg_round_precision(INTERNAL_MV_PREC, 2, &mv_cand[1][0], &mv_cand[1][1]);

  // Remove identical candidate
  if (candidates == 2 && mv_cand[0][0] == mv_cand[1][0] && mv_cand[0][1] == mv_cand[1][1]) {
    candidates = 1;
  }

  // Use Temporal Motion Vector Prediction when enabled.
  // TMVP required at least two sequential P/B-frames.
  bool can_use_tmvp =
    state->encoder_control->cfg.tmvp_enable &&
    state->frame->poc > 1 &&
    state->frame->ref->used_size &&
    candidates < AMVP_MAX_NUM_CANDS &&
    (c0 != NULL || c1 != NULL);

  if (can_use_tmvp && add_temporal_candidate(state,
                                             state->frame->ref_LX[reflist][cur_cu->inter.mv_ref[reflist]],
                                             (c0 != NULL) ? c0 : c1,
                                             reflist,
                                             mv_cand[candidates]))
  {
    candidates++;
  }

  if (candidates < AMVP_MAX_NUM_CANDS)
  {
    const uint32_t ctu_row_mul_five = ctu_row * MAX_NUM_HMVP_CANDS;
    int32_t num_cand = state->tile->frame->hmvp_size[ctu_row];
    for (int i = 0; i < MIN(/*MAX_NUM_HMVP_AVMPCANDS*/4,num_cand); i++) {
      for (int pred_source = 0; pred_source < 2; pred_source++) {
        const int cand_list = pred_source == 0 ? reflist : !reflist;
        cu_info_t* cand = &state->tile->frame->hmvp_lut[ctu_row_mul_five + num_cand - 1 - i];
        if ((cand->inter.mv_dir & (1 << cand_list)) == 0) continue;
        // Make sure the candidate points to the same reference frame
        if (state->frame->ref_LX[cand_list][cand->inter.mv_ref[cand_list]] ==
            state->frame->ref_LX[reflist][cur_cu->inter.mv_ref[reflist]])
        {
          mv_cand[candidates][0] = cand->inter.mv[cand_list][0];
          mv_cand[candidates][1] = cand->inter.mv[cand_list][1];
          candidates++;
          if (candidates == AMVP_MAX_NUM_CANDS) return;
        }
      }
    }
  }


  // Fill with (0,0)
  while (candidates < AMVP_MAX_NUM_CANDS) {
    mv_cand[candidates][0] = 0;
    mv_cand[candidates][1] = 0;
    candidates++;
  }
}

/**
 * \brief Get MV prediction for current block.
 *
 * \param state     encoder state
 * \param mv_cand   Return the motion vector candidates.
 * \param cur_cu    current CU
 * \param lcu       current LCU
 * \param reflist   reflist index (either 0 or 1)
 * \param cu_loc    Size and position of current CU
 */
void uvg_inter_get_mv_cand(
  const encoder_state_t * const state,
  mv_t mv_cand[2][2],
  const cu_info_t  * const cur_cu,
  lcu_t *lcu,
  int8_t reflist,
  const cu_loc_t* const cu_loc)
{
  merge_candidates_t merge_cand = { 0 };
  const uint8_t parallel_merge_level = state->encoder_control->cfg.log2_parallel_merge_level;
  if (cur_cu->type == CU_IBC) {
    mv_t ibc_mv_cand[IBC_MRG_MAX_NUM_CANDS][2];
    get_ibc_merge_candidates(state, cur_cu,lcu,NULL, cu_loc->x, cu_loc->y, cu_loc->width, cu_loc->height,ibc_mv_cand);
    memcpy(mv_cand[0], ibc_mv_cand[0], sizeof(mv_t) * 2);
    memcpy(mv_cand[1], ibc_mv_cand[1], sizeof(mv_t) * 2);
  } else { 
    get_spatial_merge_candidates(cu_loc, state->tile->frame->width, state->tile->frame->height, lcu,
                                 &merge_cand,
                                 parallel_merge_level,
                                 state->encoder_control->cfg.wpp);
    get_temporal_merge_candidates(state, cu_loc, 1, 0, &merge_cand);
    get_mv_cand_from_candidates(state, &merge_cand, cur_cu, reflist, mv_cand, cu_loc->y >> LOG2_LCU_WIDTH);
  }
    
  uvg_round_precision(INTERNAL_MV_PREC, 2, &mv_cand[0][0], &mv_cand[0][1]);
  uvg_round_precision(INTERNAL_MV_PREC, 2, &mv_cand[1][0], &mv_cand[1][1]);
}

/**
 * \brief Get MV prediction for current block using state->tile->frame->cu_array.
 *
 * \param state     encoder state
 * \param mv_cand   Return the motion vector candidates.
 * \param cur_cu    current CU
 * \param reflist   reflist index (either 0 or 1)
 * \param cu_loc    Size and position of current PU/CU
 */
void uvg_inter_get_mv_cand_cua(
  const encoder_state_t * const state,
  mv_t mv_cand[2][2],
  const cu_info_t* cur_cu,
  int8_t reflist,
  const cu_loc_t* const cu_loc)
{
  merge_candidates_t merge_cand = { 0 };

  const cu_array_t *cua = state->tile->frame->cu_array;
  if (cur_cu->type == CU_IBC) {
    mv_t ibc_mv_cand[IBC_MRG_MAX_NUM_CANDS][2];
    get_ibc_merge_candidates(state, cur_cu, NULL,cua,cu_loc->x, cu_loc->y, cu_loc->width, cu_loc->height,ibc_mv_cand);
    memcpy(mv_cand[0], ibc_mv_cand[0], sizeof(mv_t) * 2);
    memcpy(mv_cand[1], ibc_mv_cand[1], sizeof(mv_t) * 2);    
  } else {
    get_spatial_merge_candidates_cua(cua,
                                     state->tile->frame->width, state->tile->frame->height, &merge_cand, state->encoder_control->cfg.wpp,
                                     cu_loc);
    get_temporal_merge_candidates(state, cu_loc, 1, 0, &merge_cand);
    get_mv_cand_from_candidates(state, &merge_cand, cur_cu, reflist, mv_cand, cu_loc->y >> LOG2_LCU_WIDTH);
  }

  uvg_round_precision(INTERNAL_MV_PREC, 2, &mv_cand[0][0], &mv_cand[0][1]);
  uvg_round_precision(INTERNAL_MV_PREC, 2, &mv_cand[1][0], &mv_cand[1][1]);
}

/**
•	\brief Checks if two CUs have similar motion vectors. The function takes two CUs and compares their motion vectors.
•	\param cu1   first CU
•	\param cu2   second CU
•	\return returns 0 if the two CUs have dissimilar motion vectors, and 1 if the motions are similar. 
*/

static bool is_duplicate_candidate(const cu_info_t* cu1, const cu_info_t* cu2)
{
  if (!cu2) return false;
  if (cu1->inter.mv_dir != cu2->inter.mv_dir) return false;

  for (int reflist = 0; reflist < 2; reflist++) {
    if (cu1->inter.mv_dir & (1 << reflist)) {
      if (cu1->inter.mv[reflist][0]  != cu2->inter.mv[reflist][0]  ||
          cu1->inter.mv[reflist][1]  != cu2->inter.mv[reflist][1]  ||
          cu1->inter.mv_ref[reflist] != cu2->inter.mv_ref[reflist]) {
        return false;
      }
    }
  }

  return true;
}

/**
* Adds a merge candidate to the list of possible candidates, if it is not a duplicate.
*
* \param cand The candidate to be added.
* \param possible_duplicate1 The first possible duplicate candidate to check for duplication.
* \param possible_duplicate2 The second possible duplicate candidate to check for duplication.
* \param merge_cand_out The output parameter to store the merge candidate information.
*
* @return Returns true if the merge candidate was added successfully, false otherwise.
*/
static bool add_merge_candidate(const cu_info_t *cand,
                                const cu_info_t *possible_duplicate1,
                                const cu_info_t *possible_duplicate2,
                                inter_merge_cand_t *merge_cand_out)
{
  if (!cand ||
      is_duplicate_candidate(cand, possible_duplicate1) ||
      is_duplicate_candidate(cand, possible_duplicate2)) {
    return false;
  }

  merge_cand_out->mv[0][0] = cand->inter.mv[0][0];
  merge_cand_out->mv[0][1] = cand->inter.mv[0][1];
  merge_cand_out->mv[1][0] = cand->inter.mv[1][0];
  merge_cand_out->mv[1][1] = cand->inter.mv[1][1];
  merge_cand_out->ref[0]   = cand->inter.mv_ref[0]; // L0/L1 references
  merge_cand_out->ref[1]   = cand->inter.mv_ref[1];
  merge_cand_out->dir      = cand->inter.mv_dir;
  return true;
}


static void hmvp_shift_lut(cu_info_t* lut, int32_t size, int32_t start, int32_t end) {

  if (end > MAX_NUM_HMVP_CANDS-1) end = MAX_NUM_HMVP_CANDS-1;
  if (end == 0 && size == 1) end = 1;
  for (int i = end-1; i >= start; i--) {
    memcpy(&lut[i + 1], &lut[i], sizeof(cu_info_t));
  }
}

static bool hmvp_push_lut_item(cu_info_t* lut, int32_t size, const cu_info_t* cu, bool ibc) {

  int8_t duplicate = -1;

  if (ibc) {
    for (int i = 0; i < size; i++) {
      if (is_duplicate_candidate_ibc(cu, (const cu_info_t *)&lut[i])) {
        duplicate = i;
        break;
      }
    }
  } else {
    for (int i = 0; i < size; i++) {
      if (is_duplicate_candidate(cu, (const cu_info_t *)&lut[i])) {
        duplicate = i;
        break;
      }
    }
  }
  // If duplicate found, shift the whole lut up to the duplicate, otherwise to the end
  if(duplicate != 0) hmvp_shift_lut(lut, size, 0, duplicate == -1 ? MAX_NUM_HMVP_CANDS : duplicate);

  memcpy(lut, cu, sizeof(cu_info_t));

  return duplicate == -1;
}

/**
 * \brief Add processed CU motion vector to the HMVP look-up table
 *
 * \param state         encoder state
 * \param pic_x         block x position in pixels
 * \param pic_y         block y position in pixels
 * \param block_width   block width in pixels
 * \param block_height  block height in pixels
 * \param cu            current CU
 */
void uvg_hmvp_add_mv(const encoder_state_t* const state, uint32_t pic_x, uint32_t pic_y, uint32_t block_width, uint32_t block_height, const cu_info_t* cu)
{
  //if (!cu.geoFlag && !cu.affine)
  if(cu->type != CU_INTRA)
  {    
    assert((cu->type != CU_IBC || block_width * block_height > 16) && "Do not add IBC hmvp for small blocks");
    const uint8_t parallel_merge_level = state->encoder_control->cfg.log2_parallel_merge_level;
    const uint32_t xBr = block_width + pic_x;
    const uint32_t yBr = block_height + pic_y;
    bool hmvp_possible = ((xBr >> parallel_merge_level) > (pic_x >> parallel_merge_level)) && ((yBr >> parallel_merge_level) > (pic_y >> parallel_merge_level));
    if (hmvp_possible || cu->type == CU_IBC) {
      const uint32_t ctu_row = (pic_y >> LOG2_LCU_WIDTH);
      const uint32_t ctu_row_mul_five = ctu_row * MAX_NUM_HMVP_CANDS;

      
      if (cu->type == CU_IBC) {
        bool add_row = hmvp_push_lut_item(&state->tile->frame->hmvp_lut_ibc[ctu_row_mul_five], state->tile->frame->hmvp_size_ibc[ctu_row], cu, true);
        if(add_row && state->tile->frame->hmvp_size_ibc[ctu_row] < MAX_NUM_HMVP_CANDS) {
          state->tile->frame->hmvp_size_ibc[ctu_row]++;
        }
      } else {
        bool add_row = hmvp_push_lut_item(&state->tile->frame->hmvp_lut[ctu_row_mul_five], state->tile->frame->hmvp_size[ctu_row], cu, false);
        if(add_row && state->tile->frame->hmvp_size[ctu_row] < MAX_NUM_HMVP_CANDS) {
          state->tile->frame->hmvp_size[ctu_row]++;
        }
      }
    }
  }
}

static void round_avg_mv(mv_t* mvx, mv_t* mvy, int nShift)
{
  const int nOffset = 1 << (nShift - 1);
  *mvx = (*mvx + nOffset - (*mvx >= 0)) >> nShift;
  *mvy = (*mvy + nOffset - (*mvy >= 0)) >> nShift;
}

static bool different_mer(int32_t x, int32_t y, int32_t x2, int32_t y2, uint8_t parallel_merge_level) {

  if ((x >> parallel_merge_level) != (x2 >> parallel_merge_level) || (y >> parallel_merge_level) != (y2 >> parallel_merge_level))
  {
    return true;
  }

  return false;
}



void uvg_change_precision(int src, int dst, mv_t* hor, mv_t* ver) {

  const int shift = (int)dst - (int)src;
  if (shift >= 0)
  {
    uint32_t* hor_unsigned = (uint32_t *)hor;
    uint32_t* ver_unsigned = (uint32_t *)ver;

    *hor_unsigned <<= shift;
    *ver_unsigned <<= shift;
  }
  else
  {
    const int right_shift = -shift;
    const int offset = 1 << (right_shift - 1);
    *hor = *hor >= 0 ? (*hor + offset - 1) >> right_shift : (*hor + offset) >> right_shift;
    *ver = *ver >= 0 ? (*ver + offset - 1) >> right_shift : (*ver + offset) >> right_shift;
  }
}

void uvg_change_precision_vector2d(int src, int dst, vector2d_t *mv) {

  const int shift = (int)dst - (int)src;
  if (shift >= 0)
  {
    int* hor_unsigned = &mv->x;
    int* ver_unsigned = &mv->y;

    *hor_unsigned <<= shift;
    *ver_unsigned <<= shift;
  }
  else
  {
    const int right_shift = -shift;
    const int offset = 1 << (right_shift - 1);
    mv->x = mv->x >= 0 ? (mv->x + offset - 1) >> right_shift : (mv->x + offset) >> right_shift;
    mv->y = mv->y >= 0 ? (mv->y + offset - 1) >> right_shift : (mv->y + offset) >> right_shift;
  }
}

void uvg_round_precision(int src, int dst, mv_t* hor, mv_t* ver) {
  uvg_change_precision(src, dst, hor, ver);
  uvg_change_precision(dst, src, hor, ver);
}

void uvg_round_precision_vector2d(int src, int dst, vector2d_t* mv) {
  mv_t hor = mv->x;
  mv_t ver = mv->y;
  uvg_change_precision(src, dst, &hor, &ver);
  uvg_change_precision(dst, src, &hor, &ver);
  mv->x = hor;
  mv->y = ver;
}

/**
 * \brief Get merge predictions for current block
 * \param state     the encoder state
 * \param cu_loc    Size and position of current PU/CU
 * \param mv_cand   Returns the merge candidates.
 * \param lcu       lcu containing the block
 * \return          number of merge candidates
 */
uint8_t uvg_inter_get_merge_cand(
  const encoder_state_t * const state,
  const cu_loc_t* const cu_loc,
  inter_merge_cand_t mv_cand[MRG_MAX_NUM_CANDS],
  lcu_t *lcu)
{
  uint8_t candidates = 0;
  int8_t zero_idx = 0;
  const uint8_t parallel_merge_level = state->encoder_control->cfg.log2_parallel_merge_level;
  merge_candidates_t merge_cand = { 0 };
  const uint8_t max_num_cands = state->encoder_control->cfg.max_merge;
  // Current CU
  cu_info_t         *cur_cu = LCU_GET_CU_AT_PX(lcu, SUB_SCU(cu_loc->x), SUB_SCU(cu_loc->y));

  if(cur_cu->type == CU_IBC) {
    mv_t ibc_mv_cand[IBC_MRG_MAX_NUM_CANDS][2];
    get_ibc_merge_candidates(state, cur_cu,lcu,NULL, cu_loc->x, cu_loc->y, cu_loc->width, cu_loc->height,ibc_mv_cand);
    for (int i = 0; i < IBC_MRG_MAX_NUM_CANDS; i++) {
      mv_cand[i].dir = 1;
      mv_cand[i].mv[0][0] = ibc_mv_cand[i][0];
      mv_cand[i].mv[0][1] = ibc_mv_cand[i][1];
    }
    return IBC_MRG_MAX_NUM_CANDS;
  }
  get_spatial_merge_candidates(cu_loc, state->tile->frame->width, state->tile->frame->height, lcu,
                               &merge_cand,
                               parallel_merge_level,
                               state->encoder_control->cfg.wpp);

  const cu_info_t **a = merge_cand.a;
  const cu_info_t **b = merge_cand.b;

  const int x = cu_loc->x;
  const int y = cu_loc->y;

  if (different_mer(x, y, x, y - 1, parallel_merge_level) && add_merge_candidate(b[1], NULL, NULL, &mv_cand[candidates])) candidates++;
  if (different_mer(x, y, x - 1, y, parallel_merge_level) && add_merge_candidate(a[1], b[1], NULL, &mv_cand[candidates])) candidates++;
  if (different_mer(x, y, x + 1, y - 1, parallel_merge_level) && add_merge_candidate(b[0], b[1], NULL, &mv_cand[candidates])) candidates++;
  if (different_mer(x, y, x - 1, y + 1, parallel_merge_level) && add_merge_candidate(a[0], a[1], NULL, &mv_cand[candidates])) candidates++;
  if (candidates < 4 &&
     different_mer(x, y, x - 1, y - 1, parallel_merge_level) && add_merge_candidate(b[2], a[1], b[1], &mv_cand[candidates])) candidates++;

  bool can_use_tmvp =
    state->encoder_control->cfg.tmvp_enable &&
    candidates < max_num_cands &&
    state->frame->ref->used_size;

  if (can_use_tmvp) {
    mv_cand[candidates].dir = 0;

    const int max_reflist = (state->frame->slicetype == UVG_SLICE_B ? 1 : 0);
    for (int reflist = 0; reflist <= max_reflist; reflist++) {
      // Fetch temporal candidates for the current CU
      // ToDo: change collocated_from_l0_flag to allow L1 ref
      get_temporal_merge_candidates(state, cu_loc, 1, 0, &merge_cand);
      // TODO: enable L1 TMVP candidate
      // get_temporal_merge_candidates(state, x, y, width, height, 2, 0, &merge_cand);

      const cu_info_t *temporal_cand =
        (merge_cand.c0 != NULL) ? merge_cand.c0 : merge_cand.c1;

      if (add_temporal_candidate(state,
                                 // Reference index 0 is always used for
                                 // the temporal merge candidate.
                                 state->frame->ref_LX[0][0],
                                 temporal_cand,
                                 reflist,
                                 mv_cand[candidates].mv[reflist])) {
        mv_cand[candidates].ref[reflist] = 0;
        mv_cand[candidates].dir |= (1 << reflist);

        if (state->frame->ref->pocs[state->frame->ref_LX[reflist][0]] > state->frame->poc) {
          mv_cand[candidates].mv[reflist][0] *= -1;
          mv_cand[candidates].mv[reflist][1] *= -1;
        }
      }
    }
    
    if (mv_cand[candidates].dir != 0) {
      candidates++;
    }
  }

  if (candidates == max_num_cands) return candidates;

  if (candidates != max_num_cands - 1) {
    const uint32_t ctu_row = (cu_loc->y >> LOG2_LCU_WIDTH);
    const uint32_t ctu_row_mul_five = ctu_row * MAX_NUM_HMVP_CANDS;
    int32_t num_cand = state->tile->frame->hmvp_size[ctu_row];

    for (int i = 0; i < num_cand; i++) {
      const cu_info_t* hmvp_cand = &state->tile->frame->hmvp_lut[ctu_row_mul_five + i];
      if (i > 1 || ((!is_duplicate_candidate(hmvp_cand, a[1]))
                 && (!is_duplicate_candidate(hmvp_cand, b[1]))) ) {
        mv_cand[candidates].mv[0][0] = state->tile->frame->hmvp_lut[ctu_row_mul_five + i].inter.mv[0][0];
        mv_cand[candidates].mv[0][1] = state->tile->frame->hmvp_lut[ctu_row_mul_five + i].inter.mv[0][1];
        mv_cand[candidates].dir = state->tile->frame->hmvp_lut[ctu_row_mul_five + i].inter.mv_dir;
        mv_cand[candidates].ref[0] = state->tile->frame->hmvp_lut[ctu_row_mul_five + i].inter.mv_ref[0];
        if (state->frame->slicetype == UVG_SLICE_B) {
          mv_cand[candidates].mv[1][0] = state->tile->frame->hmvp_lut[ctu_row_mul_five + i].inter.mv[1][0];
          mv_cand[candidates].mv[1][1] = state->tile->frame->hmvp_lut[ctu_row_mul_five + i].inter.mv[1][1];
          mv_cand[candidates].ref[1] = state->tile->frame->hmvp_lut[ctu_row_mul_five + i].inter.mv_ref[1];
        }
        candidates++;
        if (candidates == max_num_cands - 1) break;
      }
    }
  }

  // pairwise-average candidates
  if (candidates > 1 && candidates < max_num_cands)
  {
    // calculate average MV for L0 and L1 seperately
    uint8_t inter_dir = 0;

    for (int reflist = 0; reflist < (state->frame->slicetype == UVG_SLICE_B ? 2 : 1); reflist++)
    {
      const int16_t ref_i = mv_cand[0].dir & (reflist + 1) ? mv_cand[0].ref[reflist] : -1;
      const int16_t ref_j = mv_cand[1].dir & (reflist + 1) ? mv_cand[1].ref[reflist] : -1;

      // both MVs are invalid, skip
      if ((ref_i == -1) && (ref_j == -1))
      {
        continue;
      }

      inter_dir += 1 << reflist;
      // both MVs are valid, average these two MVs
      if ((ref_i != -1) && (ref_j != -1))
      {
        mv_t mv_i[2] = { mv_cand[0].mv[reflist][0], mv_cand[0].mv[reflist][1] };
        mv_t mv_j[2] = { mv_cand[1].mv[reflist][0], mv_cand[1].mv[reflist][1] };

        // average two MVs
        mv_t avg_mv[2] = { mv_i[0] + mv_j[0], mv_i[1] + mv_j[1] };
        
        round_avg_mv(&avg_mv[0], &avg_mv[1], 1);

        mv_cand[candidates].mv[reflist][0] = avg_mv[0];
        mv_cand[candidates].mv[reflist][1] = avg_mv[1];
        mv_cand[candidates].ref[reflist] = (uint8_t)ref_i;
      }
      // only one MV is valid, take the only one MV
      else if (ref_i != -1)
      {
        mv_t mv_i[2] = { mv_cand[0].mv[reflist][0], mv_cand[0].mv[reflist][1] };

        mv_cand[candidates].mv[reflist][0] = mv_i[0];
        mv_cand[candidates].mv[reflist][1] = mv_i[1];
        mv_cand[candidates].ref[reflist] = (uint8_t)ref_i;
      }
      else if (ref_j != -1)
      {
        mv_t mv_j[2] = { mv_cand[1].mv[reflist][0], mv_cand[1].mv[reflist][1] };

        mv_cand[candidates].mv[reflist][0] = mv_j[0];
        mv_cand[candidates].mv[reflist][1] = mv_j[1];
        mv_cand[candidates].ref[reflist] = (uint8_t)ref_j;
      }
    }

    mv_cand[candidates].dir = inter_dir;
    if (inter_dir > 0)
    {
      candidates++;
    }
  }

  if (candidates == max_num_cands) return candidates;


  int num_ref = state->frame->ref->used_size;

  if (candidates < max_num_cands && state->frame->slicetype == UVG_SLICE_B) {
    int ref_negative = 0;
    int ref_positive = 0;
    for (uint32_t j = 0; j < state->frame->ref->used_size; j++) {
      if (state->frame->ref->pocs[j] < state->frame->poc) {
        ref_negative++;
      } else {
        ref_positive++;
      }
    }
    num_ref = MIN(ref_negative, ref_positive);
  }

  // Add (0,0) prediction
  while (candidates != max_num_cands) {
    mv_cand[candidates].mv[0][0] = 0;
    mv_cand[candidates].mv[0][1] = 0;
    mv_cand[candidates].ref[0] = (zero_idx >= num_ref - 1) ? 0 : zero_idx;
    mv_cand[candidates].dir = 1;
    if (state->frame->slicetype == UVG_SLICE_B) {
      mv_cand[candidates].ref[1] = mv_cand[candidates].ref[0];
      mv_cand[candidates].mv[1][0] = 0;
      mv_cand[candidates].mv[1][1] = 0;
      mv_cand[candidates].dir = 3;
    }
    zero_idx++;
    candidates++;
  }

  return candidates;
}
