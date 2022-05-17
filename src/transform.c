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

#include "transform.h"

#include "image.h"
#include "uvg266.h"
#include "lfnst_tables.h"
#include "rdo.h"
#include "strategies/strategies-dct.h"
#include "strategies/strategies-quant.h"
#include "strategies/strategies-picture.h"
#include "tables.h"
#include "reshape.h"

/**
 * \brief RDPCM direction.
 */
typedef enum rdpcm_dir {
  RDPCM_VER = 0, // vertical
  RDPCM_HOR = 1, // horizontal
} rdpcm_dir;

//////////////////////////////////////////////////////////////////////////
// INITIALIZATIONS
//


const uint8_t uvg_g_chroma_scale[58]=
{
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,
  17,18,19,20,21,22,23,24,25,26,27,28,29,29,30,31,32,
  33,33,34,34,35,35,36,36,37,37,38,39,40,41,42,43,44,
  45,46,47,48,49,50,51
};

//////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//

/**
 * \brief Bypass transform and quantization.
 *
 * Copies the reference pixels directly to reconstruction and the residual
 * directly to coefficients. Used when cu_transquant_bypass_flag is set.
 * Parameters pred_in and rec_out may be aliased.
 *
 * \param width       Transform width.
 * \param in_stride   Stride for ref_in and pred_in
 * \param out_stride  Stride for rec_out.
 * \param ref_in      Reference pixels.
 * \param pred_in     Predicted pixels.
 * \param rec_out     Returns the reconstructed pixels.
 * \param coeff_out   Returns the coefficients used for reconstruction of rec_out.
 *
 * \returns  Whether coeff_out contains any non-zero coefficients.
 */
static bool bypass_transquant(const int width,
                              const int in_stride,
                              const int out_stride,
                              const uvg_pixel *const ref_in,
                              const uvg_pixel *const pred_in,
                              uvg_pixel *rec_out,
                              coeff_t *coeff_out)
{
  bool nonzero_coeffs = false;

  for (int y = 0; y < width; ++y) {
    for (int x = 0; x < width; ++x) {
      int32_t in_idx    = x + y * in_stride;
      int32_t out_idx   = x + y * out_stride;
      int32_t coeff_idx = x + y * width;

      // The residual must be computed before writing to rec_out because
      // pred_in and rec_out may point to the same array.
      coeff_t coeff        = (coeff_t)(ref_in[in_idx] - pred_in[in_idx]);
      coeff_out[coeff_idx] = coeff;
      rec_out[out_idx]     = ref_in[in_idx];

      nonzero_coeffs |= (coeff != 0);
    }
  }

  return nonzero_coeffs;
}

/**
 * Apply DPCM to residual.
 *
 * \param width   width of the block
 * \param dir     RDPCM direction
 * \param coeff   coefficients (residual) to filter
 */
static void rdpcm(const int width,
                  const rdpcm_dir dir,
                  coeff_t *coeff)
{
  const int offset = (dir == RDPCM_HOR) ? 1 : width;
  const int min_x  = (dir == RDPCM_HOR) ? 1 : 0;
  const int min_y  = (dir == RDPCM_HOR) ? 0 : 1;

  for (int y = width - 1; y >= min_y; y--) {
    for (int x = width - 1; x >= min_x; x--) {
      const int index = x + y * width;
      coeff[index] -= coeff[index - offset];
    }
  }
}

/**
 * \brief Get scaled QP used in quantization
 *
 */
int32_t uvg_get_scaled_qp(color_t color, int8_t qp, int8_t qp_offset, int8_t const * const chroma_scale)
{
  int32_t qp_scaled = 0;
  if(color == 0) {
    qp_scaled = qp + qp_offset;
  } else {
    qp_scaled = CLIP(-qp_offset, 57, qp);
    if (chroma_scale) {
      qp_scaled = chroma_scale[qp] + qp_offset;
    }
    else {
      qp_scaled = qp_scaled + qp_offset;
    } 
  }
  return qp_scaled;
}

/**
 * \brief NxN inverse transform (2D)
 * \param coeff input data (transform coefficients)
 * \param block output data (residual)
 * \param block_size input data (width of transform)
 */
void uvg_transformskip(const encoder_control_t * const encoder, int16_t *block,int16_t *coeff, int8_t block_size)
{
  int32_t  j,k;
  for (j = 0; j < block_size; j++) {
    for(k = 0; k < block_size; k ++) {
      // Casting back and forth to make UBSan not trigger due to left-shifting negatives
      coeff[j * block_size + k] = (int16_t)((uint16_t)(block[j * block_size + k]));
    }
  }
}

/**
 * \brief inverse transform skip
 * \param coeff input data (transform coefficients)
 * \param block output data (residual)
 * \param block_size width of transform
 */
void uvg_itransformskip(const encoder_control_t * const encoder, int16_t *block,int16_t *coeff, int8_t block_size)
{
  int32_t  j,k;
  for ( j = 0; j < block_size; j++ ) {
    for(k = 0; k < block_size; k ++) {
      block[j * block_size + k] =  coeff[j * block_size + k];
    }
  }
}

/**
 * \brief forward transform (2D)
 * \param block input residual
 * \param coeff transform coefficients
 * \param block_size width of transform
 */
void uvg_transform2d(const encoder_control_t * const encoder,
                     int16_t *block,
                     int16_t *coeff,
                     int8_t block_size,
                     color_t color,
                     const cu_info_t *tu)
{
  if (encoder->cfg.mts)
  {
    uvg_mts_dct(encoder->bitdepth, color, tu, block_size, block, coeff, encoder->cfg.mts);
  }
  else
  {
    dct_func *dct_func = uvg_get_dct_func(block_size, color, tu->type);
    dct_func(encoder->bitdepth, block, coeff);
  }
}

void uvg_itransform2d(const encoder_control_t * const encoder,
                      int16_t *block,
                      int16_t *coeff,
                      int8_t block_size,
                      color_t color,
                      const cu_info_t *tu)
{
  if (encoder->cfg.mts)
  {
    uvg_mts_idct(encoder->bitdepth, color, tu, block_size, coeff, block, encoder->cfg.mts);
  }
  else
  {
    dct_func *idct_func = uvg_get_idct_func(block_size, color, tu->type);
    idct_func(encoder->bitdepth, coeff, block);
  }
}

void uvg_fwd_lfnst_NxN(coeff_t *src, coeff_t *dst, const int8_t mode, const int8_t index, const int8_t size, int zero_out_size)
{
  const int8_t *tr_mat = (size > 4) ? uvg_lfnst_8x8[mode][index][0] : uvg_lfnst_4x4[mode][index][0];
  const int     tr_size = (size > 4) ? 48 : 16;
  int coef;
  coeff_t *out = dst;
  assert(index < 3 && "LFNST index must be in [0, 2]");

  for (int j = 0; j < zero_out_size; j++)
  {
    coeff_t *src_ptr = src;
    const int8_t* tr_mat_tmp = tr_mat;
    coef = 0;
    for (int i = 0; i < tr_size; i++)
    {
      coef += *src_ptr++ * *tr_mat_tmp++;
    }
    *out++ = (coeff_t)((coef + 64) >> 7);
    tr_mat += tr_size;
  }

  // Possible tr_size values 16, 48. Possible zero_out_size values 8, 16
  switch (tr_size - zero_out_size) {
    case 0:
      break;
    case 8:
      FILL_ARRAY(out, 0, 8);
      break;
    case 32:
      FILL_ARRAY(out, 0, 32);
      break;
    case 40:
      FILL_ARRAY(out, 0, 40);
      break;
    default:
      assert(false && "LFNST: This should never trip.");
  }
}

static inline bool get_transpose_flag(const int8_t intra_mode)
{
  return ((intra_mode >= NUM_LUMA_MODE) && (intra_mode >= (NUM_LUMA_MODE + (NUM_EXT_LUMA_MODE >> 1)))) ||
         ((intra_mode < NUM_LUMA_MODE) && (intra_mode > DIA_IDX));
}

void uvg_fwd_lfnst(const cu_info_t* const cur_cu,
  const int width, const int height,
  const uint8_t color,
  const uint16_t lfnst_idx,
  coeff_t *coeffs)
{
  const uint16_t lfnst_index = lfnst_idx;
  int8_t intra_mode = (color == COLOR_Y) ? cur_cu->intra.mode : cur_cu->intra.mode_chroma;
  bool mts_skip = cur_cu->tr_idx == MTS_SKIP;
  const int depth = cur_cu->depth;
  bool is_separate_tree = depth == 4; // TODO: proper dual tree check when that structure is implemented
  bool is_cclm_mode = (intra_mode >= 81 && intra_mode <= 83); // CCLM modes are in [81, 83]

  bool is_mip = cur_cu->type == CU_INTRA ? cur_cu->intra.mip_flag : false;
  bool is_wide_angle = false; // TODO: get wide angle mode when implemented

  const int cu_type = cur_cu->type;

  const int scan_order = uvg_get_scan_order(cu_type, intra_mode, depth);

  if (lfnst_index && !mts_skip && (is_separate_tree || color == COLOR_Y))
  {
    const uint32_t log2_block_size = uvg_g_convert_to_bit[width] + 2;
    assert(log2_block_size != -1 && "LFNST: invalid block width.");
    const bool whge3 = width >= 8 && height >= 8;
    const uint32_t* scan = whge3 ? uvg_coef_top_left_diag_scan_8x8[log2_block_size] : uvg_g_sig_last_scan[scan_order][log2_block_size - 1];

    if (is_cclm_mode) {
      intra_mode = cur_cu->intra.mode;
    }
    if (is_mip) {
      intra_mode = 0; // Set to planar mode
    }
    assert(intra_mode < NUM_INTRA_MODE && "LFNST: Invalid intra mode.");
    assert(lfnst_index < 3 && "LFNST: Invalid LFNST index. Must be in [0, 2]");

    if (is_wide_angle) {
      // Transform wide angle mode to intra mode
      intra_mode = intra_mode; // TODO: wide angle modes not implemented yet. Do nothing.
    }

    bool transpose = get_transpose_flag(intra_mode);
    const int sb_size = whge3 ? 8 : 4;
    bool tu_4x4 = (width == 4 && height == 4);
    bool tu_8x8 = (width == 8 && height == 8);
 
    coeff_t tmp_in_matrix[48];
    coeff_t tmp_out_matrix[48];
    coeff_t *lfnst_tmp = tmp_in_matrix; // forward low frequency non-separable transform
      
    coeff_t *coeff_tmp = coeffs;

    int y;
    if (transpose) {
      if (sb_size == 4) {
        for (y = 0; y < 4; y++) {
          lfnst_tmp[0] = coeff_tmp[0];
          lfnst_tmp[4] = coeff_tmp[1];
          lfnst_tmp[8] = coeff_tmp[2];
          lfnst_tmp[12] = coeff_tmp[3];
          lfnst_tmp++;
          coeff_tmp += width;
        }
      }
      else { // ( sb_size == 8 )
        for (y = 0; y < 8; y++) {
          lfnst_tmp[0] = coeff_tmp[0];
          lfnst_tmp[8] = coeff_tmp[1];
          lfnst_tmp[16] = coeff_tmp[2];
          lfnst_tmp[24] = coeff_tmp[3];
          if (y < 4) {
            lfnst_tmp[32] = coeff_tmp[4];
            lfnst_tmp[36] = coeff_tmp[5];
            lfnst_tmp[40] = coeff_tmp[6];
            lfnst_tmp[44] = coeff_tmp[7];
          }
          lfnst_tmp++;
          coeff_tmp += width;
        }
      }
    }
    else {
      for (y = 0; y < sb_size; y++) {
        uint32_t stride = (y < 4) ? sb_size : 4;
        memcpy(lfnst_tmp, coeff_tmp, stride * sizeof(coeff_t));
        lfnst_tmp += stride;
        coeff_tmp += width;
      }
    }

    uvg_fwd_lfnst_NxN(tmp_in_matrix, tmp_out_matrix, uvg_lfnst_lut[intra_mode], lfnst_index - 1, sb_size,
      (tu_4x4 || tu_8x8) ? 8 : 16);

    lfnst_tmp = tmp_out_matrix;   // forward spectral rearrangement
    coeff_tmp = coeffs;
    int lfnst_coeff_num = (sb_size == 4) ? sb_size * sb_size : 48;

    const uint32_t *scan_ptr = scan;

    for (y = 0; y < lfnst_coeff_num; y++) {
      coeff_tmp[*scan_ptr] = *lfnst_tmp++;
      scan_ptr++;
    }
  }
}

void uvg_inv_lfnst_NxN(coeff_t *src, coeff_t *dst, const uint32_t mode, const uint32_t index, const uint32_t size, int zero_out_size, const int max_log2_tr_dyn_range)
{
  const coeff_t output_min = -(1 << max_log2_tr_dyn_range);
  const coeff_t output_max = (1 << max_log2_tr_dyn_range) - 1;
  const int8_t *tr_mat = (size > 4) ? uvg_lfnst_8x8[mode][index][0] : uvg_lfnst_4x4[mode][index][0];
  const int tr_size = (size > 4) ? 48 : 16;
  int resi;
  coeff_t *out = dst;
  assert(index < 3);

  for (int j = 0; j < tr_size; j++)
  {
    resi = 0;
    const int8_t* tr_mat_tmp = tr_mat;
    coeff_t *src_ptr = src;
    for (int i = 0; i < zero_out_size; i++)
    {
      resi += *src_ptr++ * *tr_mat_tmp;
      tr_mat_tmp += tr_size;
    }
    *out++ = CLIP(output_min, output_max, (coeff_t)((resi + 64) >> 7));
    tr_mat++;
  }
}

void uvg_inv_lfnst(const cu_info_t *cur_cu, 
                   const int width, const int height,
                   const uint8_t color,
                   const uint16_t lfnst_idx,
                   coeff_t *coeffs)
{
  // In VTM, max log2 dynamic range is something in range [15, 20] depending on whether extended precision processing is enabled
  // Such is not yet present in uvg266 so use 15 for now
  const int max_log2_dyn_range = 15;
  const uint32_t  lfnst_index = lfnst_idx;
  int8_t intra_mode = (color == COLOR_Y) ? cur_cu->intra.mode : cur_cu->intra.mode_chroma;
  bool mts_skip = cur_cu->tr_idx == MTS_SKIP;
  const int depth = cur_cu->depth;
  bool is_separate_tree = depth == 4; // TODO: proper dual tree check when that structure is implemented
  bool is_cclm_mode = (intra_mode >= 81 && intra_mode <= 83); // CCLM modes are in [81, 83]

  bool is_mip = cur_cu->type == CU_INTRA ? cur_cu->intra.mip_flag : false;
  bool is_wide_angle = false; // TODO: get wide angle mode when implemented

  const int cu_type = cur_cu->type;

  const int scan_order = uvg_get_scan_order(cu_type, intra_mode, depth);
  
  if (lfnst_index && !mts_skip && (is_separate_tree || color == COLOR_Y)) {
    const uint32_t log2_block_size = uvg_g_convert_to_bit[width] + 2;
    const bool whge3 = width >= 8 && height >= 8;
    const uint32_t* scan = whge3 ? uvg_coef_top_left_diag_scan_8x8[log2_block_size] : uvg_g_sig_last_scan[scan_order][log2_block_size - 1];
    
    if (is_cclm_mode) {
      intra_mode = cur_cu->intra.mode;
    }
    if (is_mip) {
      intra_mode = 0; // Set to planar mode
    }
    assert(intra_mode < NUM_INTRA_MODE && "LFNST: Invalid intra mode.");
    assert(lfnst_index < 3 && "LFNST: Invalid LFNST index. Must be in [0, 2]");

    if (is_wide_angle) {
      // Transform wide angle mode to intra mode
      intra_mode = intra_mode; // TODO: wide angle modes not implemented yet. Do nothing.
    }

    bool          transpose_flag = get_transpose_flag(intra_mode);
    const int     sb_size = whge3 ? 8 : 4;
    bool          tu_4x4_flag = (width == 4 && height == 4);
    bool          tu_8x8_flag = (width == 8 && height == 8);
    coeff_t tmp_in_matrix[48];
    coeff_t tmp_out_matrix[48];
    coeff_t *lfnst_tmp;
    coeff_t *coeff_tmp;
    int           y;
    lfnst_tmp = tmp_in_matrix;   // inverse spectral rearrangement
    coeff_tmp = coeffs;
    coeff_t *dst = lfnst_tmp;

    const uint32_t *scan_ptr = scan;
    for (y = 0; y < 16; y++) {
      *dst++ = coeff_tmp[*scan_ptr];
      scan_ptr++;
    }

    uvg_inv_lfnst_NxN(tmp_in_matrix, tmp_out_matrix, uvg_lfnst_lut[intra_mode], lfnst_index - 1, sb_size,
      (tu_4x4_flag || tu_8x8_flag) ? 8 : 16, max_log2_dyn_range);
    lfnst_tmp = tmp_out_matrix;   // inverse low frequency non-separale transform

    if (transpose_flag) {
      if (sb_size == 4) {
        for (y = 0; y < 4; y++) {
          coeff_tmp[0] = lfnst_tmp[0];
          coeff_tmp[1] = lfnst_tmp[4];
          coeff_tmp[2] = lfnst_tmp[8];
          coeff_tmp[3] = lfnst_tmp[12];
          lfnst_tmp++;
          coeff_tmp += width;
        }
      }
      else { // ( sb_size == 8 )
        for (y = 0; y < 8; y++) {
          coeff_tmp[0] = lfnst_tmp[0];
          coeff_tmp[1] = lfnst_tmp[8];
          coeff_tmp[2] = lfnst_tmp[16];
          coeff_tmp[3] = lfnst_tmp[24];
          if (y < 4) {
            coeff_tmp[4] = lfnst_tmp[32];
            coeff_tmp[5] = lfnst_tmp[36];
            coeff_tmp[6] = lfnst_tmp[40];
            coeff_tmp[7] = lfnst_tmp[44];
          }
          lfnst_tmp++;
          coeff_tmp += width;
        }
      }
    }
    else {
      for (y = 0; y < sb_size; y++) {
        uint32_t uiStride = (y < 4) ? sb_size : 4;
        memcpy(coeff_tmp, lfnst_tmp, uiStride * sizeof(coeff_t));
        lfnst_tmp += uiStride;
        coeff_tmp += width;
      }
    }
  }
}

/**
 * \brief Like uvg_quantize_residual except that this uses trskip if that is better.
 *
 * Using this function saves one step of quantization and inverse quantization
 * compared to doing the decision separately from the actual operation.
 *
 * \param width  Transform width.
 * \param color  Color.
 * \param scan_order  Coefficient scan order.
 * \param trskip_out  Whether transform skip is used.
 * \param stride  Stride for ref_in, pred_in and rec_out.
 * \param ref_in  Reference pixels.
 * \param pred_in  Predicted pixels.
 * \param rec_out  Reconstructed pixels.
 * \param coeff_out  Coefficients used for reconstruction of rec_out.
 *
 * \returns  Whether coeff_out contains any non-zero coefficients.
 */
int uvg_quantize_residual_trskip(
    encoder_state_t *const state,
    const cu_info_t *const cur_cu, const int width, const color_t color,
    const coeff_scan_order_t scan_order, int8_t *trskip_out, 
    const int in_stride, const int out_stride,
    const uvg_pixel *const ref_in, const uvg_pixel *const pred_in, 
    uvg_pixel *rec_out, coeff_t *coeff_out, int lmcs_chroma_adj)
{
  struct {
    uvg_pixel rec[LCU_WIDTH * LCU_WIDTH];
    coeff_t coeff[LCU_WIDTH * LCU_WIDTH];
    double cost;
    int has_coeffs;
  } skip, *best;
  
  //noskip.has_coeffs = uvg_quantize_residual(
  //    state, cur_cu, width, color, scan_order,
  //    0, in_stride, 4,
  //    ref_in, pred_in, noskip.rec, noskip.coeff, false);
  //noskip.cost = uvg_pixels_calc_ssd(ref_in, noskip.rec, in_stride, 4, 4);
  //noskip.cost += uvg_get_coeff_cost(state, noskip.coeff, 4, 0, scan_order) * bit_cost;

  skip.has_coeffs = uvg_quantize_residual(
    state, cur_cu, width, color, scan_order,
    1, in_stride, width,
    ref_in, pred_in, skip.rec, skip.coeff, false, lmcs_chroma_adj);
  skip.cost = uvg_pixels_calc_ssd(ref_in, skip.rec, in_stride, width, width);
  skip.cost += uvg_get_coeff_cost(state, skip.coeff, NULL, width, 0, scan_order, 1) * state->frame->lambda;

/*  if (noskip.cost <= skip.cost) {
    *trskip_out = 0;
    best = &noskip;
  } else */{
    *trskip_out = 1;
    best = &skip;
  }

  if (best->has_coeffs || rec_out != pred_in) {
    // If there is no residual and reconstruction is already in rec_out, 
    // we can skip this.
    uvg_pixels_blit(best->rec, rec_out, width, width, width, out_stride);
  }
  copy_coeffs(best->coeff, coeff_out, width);

  return best->has_coeffs;
}

/**
 * Calculate the residual coefficients for a single TU.
 *
 * \param early_skip if this is used for early skip, bypass IT and IQ
 */
static void quantize_tr_residual(encoder_state_t * const state,
                                 const color_t color,
                                 const int32_t x,
                                 const int32_t y,
                                 const uint8_t depth,
                                 cu_info_t *cur_pu,
                                 lcu_t* lcu,
                                 bool early_skip)
{
  const uvg_config *cfg    = &state->encoder_control->cfg;
  const int32_t shift      = color == COLOR_Y ? 0 : 1;
  const vector2d_t lcu_px  = { SUB_SCU(x) >> shift, SUB_SCU(y) >> shift};

  // If luma is 4x4, do chroma for the 8x8 luma area when handling the top
  // left PU because the coordinates are correct.
  bool handled_elsewhere = color != COLOR_Y &&
                           depth == MAX_DEPTH &&
                           (x % 4 != 0 || y % 4 != 0);
  if (handled_elsewhere) {
    return;
  }

  // Clear coded block flag structures for depths lower than current depth.
  // This should ensure that the CBF data doesn't get corrupted if this function
  // is called more than once.

  int32_t tr_width;
  if (color == COLOR_Y) {
    tr_width = LCU_WIDTH >> depth;
  } else {
    const int chroma_depth = (depth == MAX_PU_DEPTH ? depth - 1 : depth);
    tr_width = LCU_WIDTH_C >> chroma_depth;
  }
  const int32_t lcu_width = LCU_WIDTH >> shift;
  const int8_t mode =
    (color == COLOR_Y) ? cur_pu->intra.mode : cur_pu->intra.mode_chroma;
  const coeff_scan_order_t scan_idx =
    uvg_get_scan_order(cur_pu->type, mode, depth);
  const int offset = lcu_px.x + lcu_px.y * lcu_width;
  const int z_index = xy_to_zorder(lcu_width, lcu_px.x, lcu_px.y);

  // Pointers to current location in arrays with prediction. The
  // reconstruction will be written to this array.
  uvg_pixel *pred = NULL;
  // Pointers to current location in arrays with reference.
  const uvg_pixel *ref = NULL;
  // Pointers to current location in arrays with quantized coefficients.
  coeff_t *coeff = NULL;

  switch (color) {
    case COLOR_Y:
      pred  = &lcu->rec.y[offset];
      ref   = &lcu->ref.y[offset];
      coeff = &lcu->coeff.y[z_index];
      break;
    case COLOR_U:
      pred = &lcu->rec.u[offset];
      ref  = &lcu->ref.u[offset];
      coeff = &lcu->coeff.u[z_index];
      break;
    case COLOR_V:
      pred = &lcu->rec.v[offset];
      ref  = &lcu->ref.v[offset];
      coeff = &lcu->coeff.v[z_index];
      break;
    default:
      break;
  }

  const bool can_use_trskip = tr_width <= (1 << state->encoder_control->cfg.trskip_max_size) &&
                              color == COLOR_Y &&
                              cfg->trskip_enable && 
                              cur_pu->tr_idx == 1;

  uint8_t has_coeffs;


  int lmcs_chroma_adj = 0;
  if (state->tile->frame->lmcs_aps->m_sliceReshapeInfo.enableChromaAdj && color != COLOR_Y) {
    lmcs_chroma_adj = uvg_calculate_lmcs_chroma_adj_vpdu_nei(state, state->tile->frame->lmcs_aps, x, y);
  }

  if (cfg->lossless) {
    has_coeffs = bypass_transquant(tr_width,
                                   lcu_width, // in stride
                                   lcu_width, // out stride
                                   ref,
                                   pred,
                                   pred,
                                   coeff);
    if (cfg->implicit_rdpcm && cur_pu->type == CU_INTRA) {
      // implicit rdpcm for horizontal and vertical intra modes
      if (mode == 18) {
        rdpcm(tr_width, RDPCM_HOR, coeff);
      } else if (mode == 50) {
        rdpcm(tr_width, RDPCM_VER, coeff);
      }
    }

  } else if (can_use_trskip) {
    int8_t tr_skip = 0;

    // Try quantization with trskip and use it if it's better.
    has_coeffs = uvg_quantize_residual_trskip(state,
                                              cur_pu,
                                              tr_width,
                                              color,
                                              scan_idx,
                                              &tr_skip,
                                              lcu_width,
                                              lcu_width,
                                              ref,
                                              pred,
                                              pred,
                                              coeff,
                                              lmcs_chroma_adj);
    cur_pu->tr_skip = tr_skip;
  } else {
    if(color == COLOR_UV) {
      has_coeffs = uvg_quant_cbcr_residual(
        state,
        cur_pu,
        tr_width,
        scan_idx,
        lcu_width,
        lcu_width,
        &lcu->ref.u[offset], &lcu->ref.v[offset],
        &lcu->rec.joint_u[offset], &lcu->rec.joint_v[offset],
        &lcu->rec.joint_u[offset], &lcu->rec.joint_v[offset],
        &lcu->coeff.joint_uv[z_index],
        early_skip,
        lmcs_chroma_adj
      );
      cur_pu->joint_cb_cr = has_coeffs;
      return;
    }

    has_coeffs = uvg_quantize_residual(state,
                                       cur_pu,
                                       tr_width,
                                       color,
                                       scan_idx,
                                       false, // tr skip
                                       lcu_width,
                                       lcu_width,
                                       ref,
                                       pred,
                                       pred,
                                       coeff,
                                       early_skip,
                                       lmcs_chroma_adj);
    
  }

  cbf_clear(&cur_pu->cbf, depth, color);
  if (has_coeffs) {
    cbf_set(&cur_pu->cbf, depth, color);
  }

}

/**
 * This function calculates the residual coefficients for a region of the LCU
 * (defined by x, y and depth) and updates the reconstruction with the
 * kvantized residual. Processes the TU tree recursively.
 *
 * Inputs are:
 * - lcu->rec   pixels after prediction for the area
 * - lcu->ref   reference pixels for the area
 * - lcu->cu    for the area
 * - early_skip if this is used for early skip, bypass IT and IQ
 *
 * Outputs are:
 * - lcu->rec               reconstruction after quantized residual
 * - lcu->coeff             quantized coefficients for the area
 * - lcu->cbf               coded block flags for the area
 * - lcu->cu.intra.tr_skip  tr skip flags for the area (in case of luma)
 */
void uvg_quantize_lcu_residual(
  encoder_state_t * const state,
  const bool luma,
  const bool chroma,
  const bool jccr,
  const int32_t x,
  const int32_t y,
  const uint8_t depth,
  cu_info_t *cur_pu,
  lcu_t* lcu,
  bool early_skip)
{
  const int32_t width = LCU_WIDTH >> depth;
  const vector2d_t lcu_px  = { SUB_SCU(x), SUB_SCU(y) };

  if (cur_pu == NULL) {
    cur_pu = LCU_GET_CU_AT_PX(lcu, lcu_px.x, lcu_px.y);
  }

  // Tell clang-analyzer what is up. For some reason it can't figure out from
  // asserting just depth.
  assert(width ==  4 ||
         width ==  8 ||
         width == 16 ||
         width == 32 ||
         width == 64);

  // Reset CBFs because CBFs might have been set
  // for depth earlier
  if (luma) {
    cbf_clear(&cur_pu->cbf, depth, COLOR_Y);
  }
  if (chroma || jccr) {
    cbf_clear(&cur_pu->cbf, depth, COLOR_U);
    cbf_clear(&cur_pu->cbf, depth, COLOR_V);
  }

  if (depth == 0 || cur_pu->tr_depth > depth) {

    // Split transform and increase depth
    const int offset = width / 2;
    const int32_t x2 = x + offset;
    const int32_t y2 = y + offset;

    // jccr is currently not supported if transform is split
    uvg_quantize_lcu_residual(state, luma, chroma, 0,  x,  y, depth + 1, NULL, lcu, early_skip);
    uvg_quantize_lcu_residual(state, luma, chroma, 0, x2,  y, depth + 1, NULL, lcu, early_skip);
    uvg_quantize_lcu_residual(state, luma, chroma, 0,  x, y2, depth + 1, NULL, lcu, early_skip);
    uvg_quantize_lcu_residual(state, luma, chroma, 0, x2, y2, depth + 1, NULL, lcu, early_skip);

    // Propagate coded block flags from child CUs to parent CU.
    uint16_t child_cbfs[3] = {
      LCU_GET_CU_AT_PX(lcu, lcu_px.x + offset, lcu_px.y         )->cbf,
      LCU_GET_CU_AT_PX(lcu, lcu_px.x,          lcu_px.y + offset)->cbf,
      LCU_GET_CU_AT_PX(lcu, lcu_px.x + offset, lcu_px.y + offset)->cbf,
    };

    if (depth <= MAX_DEPTH) {
      cbf_set_conditionally(&cur_pu->cbf, child_cbfs, depth, COLOR_Y);
      cbf_set_conditionally(&cur_pu->cbf, child_cbfs, depth, COLOR_U);
      cbf_set_conditionally(&cur_pu->cbf, child_cbfs, depth, COLOR_V);
    }

  } else {
    // Process a leaf TU.
    if (luma) {
      quantize_tr_residual(state, COLOR_Y, x, y, depth, cur_pu, lcu, early_skip);
    }
    if (chroma) {
      quantize_tr_residual(state, COLOR_U, x, y, depth, cur_pu, lcu, early_skip);
      quantize_tr_residual(state, COLOR_V, x, y, depth, cur_pu, lcu, early_skip);   
    }
    if (jccr && cur_pu->tr_depth == cur_pu->depth) {
      quantize_tr_residual(state, COLOR_UV, x, y, depth, cur_pu, lcu, early_skip);
    }
  }
}
