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

#include "strategies/generic/quant-generic.h"

#include <stdlib.h>

#include "encoder.h"
#include "uvg_math.h"
#include "rdo.h"
#include "scalinglist.h"
#include "strategies/strategies-quant.h"
#include "strategyselector.h"
#include "transform.h"
#include "fast_coeff_cost.h"
#include "reshape.h"

/**
* \brief quantize transformed coefficents
*
*/
void uvg_quant_generic(
  const encoder_state_t * const state,
  coeff_t *coef,
  coeff_t *q_coef,
  int32_t width,
  int32_t height,
  color_t color,
  int8_t scan_idx,
  int8_t block_type,
  int8_t transform_skip,
  uint8_t lfnst_idx)
{
  const encoder_control_t * const encoder = state->encoder_control;
  const uint32_t log2_tr_width  = uvg_g_convert_to_log2[width];
  const uint32_t log2_tr_height = uvg_g_convert_to_log2[height];
  const uint32_t * const scan = uvg_get_scan_order_table(SCAN_GROUP_4X4, scan_idx, log2_tr_width, log2_tr_height, 0);

  int32_t qp_scaled = uvg_get_scaled_qp(color, state->qp, (encoder->bitdepth - 8) * 6, encoder->qp_map[0]);
  qp_scaled = transform_skip ? MAX(qp_scaled, 4 + 6 * MIN_QP_PRIME_TS) : qp_scaled;
  bool needs_block_size_trafo_scale = !transform_skip && ((log2_tr_height + log2_tr_width) % 2 == 1);
  needs_block_size_trafo_scale |= 0; // Non log2 block size
    
  const int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)color;
  const int32_t *quant_coeff = encoder->scaling_list.quant_coeff[log2_tr_width][log2_tr_height][scalinglist_type][qp_scaled % 6];
  const int32_t transform_shift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - ((log2_tr_height + log2_tr_width) >> 1) - needs_block_size_trafo_scale; //!< Represents scaling through forward transform
  const int64_t q_bits = QUANT_SHIFT + qp_scaled / 6 + (transform_skip ? 0 : transform_shift );
  const int32_t add = ((state->frame->slicetype == UVG_SLICE_I) ? 171 : 85) << (q_bits - 9);
  const int32_t q_bits8 = q_bits - 8;

  const int32_t default_quant_coeff = uvg_g_quant_scales[needs_block_size_trafo_scale][qp_scaled % 6];

  uint32_t ac_sum = 0;

  const bool use_scaling_list = state->encoder_control->cfg.scaling_list != UVG_SCALING_LIST_OFF;

  if(lfnst_idx == 0){
    for (int32_t n = 0; n < width * height; n++) {
      int32_t level = coef[n];
      int64_t abs_level = (int64_t)abs(level);
      int32_t sign;

      sign = (level < 0 ? -1 : 1);

      int32_t curr_quant_coeff = use_scaling_list ? quant_coeff[n] : default_quant_coeff;
      level = (int32_t)((abs_level * curr_quant_coeff + add) >> q_bits);
      ac_sum += level;

      level *= sign;
      q_coef[n] = (coeff_t)(CLIP(-32768, 32767, level));

    }
  }
  else {
    const int max_number_of_coeffs = ((width == 4 && height == 4) || (width == 8 && height == 8)) ? 8 : 16;
    memset(q_coef, 0, width * height * sizeof(coeff_t));
    for (int32_t n = 0; n < max_number_of_coeffs; n++) {
      const uint32_t idx = scan[n];
      int32_t level = coef[idx];
      int64_t abs_level = (int64_t)abs(level);
      int32_t sign;

      sign = (level < 0 ? -1 : 1);

      int32_t curr_quant_coeff = quant_coeff[n];
      level = (abs_level * curr_quant_coeff + add) >> q_bits;
      ac_sum += level;

      level *= sign;
      q_coef[idx] = (coeff_t)(CLIP(-32768, 32767, level));
    }
  }

  // Signhiding
  if (!encoder->cfg.signhide_enable || ac_sum < 2) return;

  int32_t delta_u[LCU_WIDTH*LCU_WIDTH >> 2];

  if(lfnst_idx == 0) {
    for (int32_t n = 0; n < width * height; n++) {
      int32_t level = coef[n];
      int64_t abs_level = (int64_t)abs(level);
      int32_t curr_quant_coeff = quant_coeff[n];

      level = (int32_t)((abs_level * curr_quant_coeff + add) >> q_bits);
      delta_u[n] = (int32_t)((abs_level * curr_quant_coeff - (level << q_bits)) >> q_bits8);
    }
  }
  else {
    const int max_number_of_coeffs = ((width == 4 && height == 4) || (width == 8 && height == 8)) ? 8 : 16;
    for (int32_t n = 0; n < max_number_of_coeffs; n++) {
      const uint32_t idx = scan[n];
      int32_t level = coef[idx];
      int64_t abs_level = (int64_t)abs(level);
      int32_t curr_quant_coeff = quant_coeff[idx];

      level = (abs_level * curr_quant_coeff + add) >> q_bits;
      delta_u[idx] = (int32_t)((abs_level * curr_quant_coeff - (level << q_bits)) >> q_bits8);
    }
  }

  if (ac_sum >= 2) {
#define SCAN_SET_SIZE 16
#define LOG2_SCAN_SET_SIZE 4
    int32_t n, last_cg = -1, abssum = 0, subset, subpos;
    for (subset = (width*height - 1) >> LOG2_SCAN_SET_SIZE; subset >= 0; subset--) {
      int32_t first_nz_pos_in_cg = SCAN_SET_SIZE, last_nz_pos_in_cg = -1;
      subpos = subset << LOG2_SCAN_SET_SIZE;
      abssum = 0;

      // Find last coeff pos
      for (n = SCAN_SET_SIZE - 1; n >= 0; n--)  {
        if (q_coef[scan[n + subpos]])  {
          last_nz_pos_in_cg = n;
          break;
        }
      }

      // First coeff pos
      for (n = 0; n <SCAN_SET_SIZE; n++) {
        if (q_coef[scan[n + subpos]]) {
          first_nz_pos_in_cg = n;
          break;
        }
      }

      // Sum all uvg_quant coeffs between first and last
      for (n = first_nz_pos_in_cg; n <= last_nz_pos_in_cg; n++) {
        abssum += q_coef[scan[n + subpos]];
      }

      if (last_nz_pos_in_cg >= 0 && last_cg == -1) {
        last_cg = 1;
      }

      if (last_nz_pos_in_cg - first_nz_pos_in_cg >= 4) {
        int32_t signbit = (q_coef[scan[subpos + first_nz_pos_in_cg]] > 0 ? 0 : 1);
        if (signbit != (abssum & 0x1)) { // compare signbit with sum_parity
          int32_t min_cost_inc = 0x7fffffff, min_pos = -1, cur_cost = 0x7fffffff;
          int16_t final_change = 0, cur_change = 0;
          for (n = (last_cg == 1 ? last_nz_pos_in_cg : SCAN_SET_SIZE - 1); n >= 0; n--) {
            uint32_t blkPos = scan[n + subpos];
            if (q_coef[blkPos] != 0) {
              if (delta_u[blkPos] > 0) {
                cur_cost = -delta_u[blkPos];
                cur_change = 1;
              }
              else if (n == first_nz_pos_in_cg && abs(q_coef[blkPos]) == 1) {
                cur_cost = 0x7fffffff;
              }
              else {
                cur_cost = delta_u[blkPos];
                cur_change = -1;
              }
            }
            else if (n < first_nz_pos_in_cg && ((coef[blkPos] >= 0) ? 0 : 1) != signbit) {
              cur_cost = 0x7fffffff;
            }
            else {
              cur_cost = -delta_u[blkPos];
              cur_change = 1;
            }

            if (cur_cost < min_cost_inc) {
              min_cost_inc = cur_cost;
              final_change = cur_change;
              min_pos = blkPos;
            }
          } // CG loop

          if (q_coef[min_pos] == 32767 || q_coef[min_pos] == -32768) {
            final_change = -1;
          }

          if (coef[min_pos] >= 0) q_coef[min_pos] += final_change;
          else q_coef[min_pos] -= final_change;
        } // Hide
      }
      if (last_cg == 1) last_cg = 0;
    }

#undef SCAN_SET_SIZE
#undef LOG2_SCAN_SET_SIZE
  }
}

static INLINE int64_t square(int x) {
  return x * (int64_t)x;
}


int uvg_quant_cbcr_residual_generic(
  encoder_state_t* const state, 
  const cu_info_t* const cur_cu,
  const int width,
  const int height,
  const coeff_scan_order_t scan_order,
  const int in_stride, const int out_stride,
  const uvg_pixel* const u_ref_in, 
  const uvg_pixel* const v_ref_in, 
  const uvg_pixel* const u_pred_in,
  const uvg_pixel* const v_pred_in,
  uvg_pixel* u_rec_out,
  uvg_pixel* v_rec_out,
  coeff_t* coeff_out,
  bool early_skip, 
  int lmcs_chroma_adj, enum uvg_tree_type tree_type) 
{
  ALIGNED(64) int16_t u_residual[TR_MAX_WIDTH * TR_MAX_WIDTH];
  ALIGNED(64) int16_t v_residual[TR_MAX_WIDTH * TR_MAX_WIDTH];
  ALIGNED(64) int16_t combined_residual[TR_MAX_WIDTH * TR_MAX_WIDTH];
  ALIGNED(64) coeff_t coeff[TR_MAX_WIDTH * TR_MAX_WIDTH];
  // TODO: this function is not fully converted to handle non-square blocks
  {
    int y, x;
    for (y = 0; y < height; ++y) {
      for (x = 0; x < width; ++x) {
        u_residual[x + y * width] = (int16_t)(u_ref_in[x + y * in_stride] - u_pred_in[x + y * in_stride]);
        v_residual[x + y * width] = (int16_t)(v_ref_in[x + y * in_stride] - v_pred_in[x + y * in_stride]);
      }
    }
  }
  uvg_generate_residual(u_ref_in, u_pred_in, u_residual, width, height, in_stride, in_stride);
  uvg_generate_residual(v_ref_in, v_pred_in, v_residual, width, height, in_stride, in_stride);
  
  
  const int cbf_mask = cur_cu->joint_cb_cr * (state->frame->jccr_sign ? -1 : 1);
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      const int16_t cbx = u_residual[x + y * width], crx = v_residual[x + y * width];
      if (cbf_mask == 2)
      {
        combined_residual[x + y * width] = (4 * cbx + 2 * crx) / 5;
      }
      else if (cbf_mask == -2)
      {
        combined_residual[x + y * width] = (4 * cbx - 2 * crx) / 5;
      }
      else if (cbf_mask == 3)
      {
        combined_residual[x + y * width] = (cbx + crx) / 2;
      }
      else if (cbf_mask == -3)
      {
        combined_residual[x + y * width] = (cbx - crx) / 2;
      }
      else if (cbf_mask == 1)
      {
        combined_residual[x + y * width] = (4 * crx + 2 * cbx) / 5;
      }
      else if (cbf_mask == -1)
      {
        combined_residual[x + y * width] = (4 * crx - 2 * cbx) / 5;
      }
      else
      {
        assert(0);
      }
    }
  }


  uvg_transform2d(state->encoder_control, combined_residual, coeff, width, height, cur_cu->joint_cb_cr == 1 ? COLOR_V : COLOR_U, cur_cu);
  uint8_t lfnst_idx = tree_type == UVG_CHROMA_T ? cur_cu->cr_lfnst_idx : cur_cu->lfnst_idx;
  if(lfnst_idx) {
    uvg_fwd_lfnst(cur_cu, width, height, COLOR_UV, lfnst_idx, coeff, tree_type, state->collocated_luma_mode);
  }
  int abs_sum = 0;
  if (!false && state->encoder_control->cfg.dep_quant) {
    uvg_dep_quant(
      state,
      cur_cu,
      width,
      height,
      coeff,
      coeff_out,
      COLOR_U,
      tree_type,
      &abs_sum,
      state->encoder_control->cfg.scaling_list);
  }
  else if (state->encoder_control->cfg.rdoq_enable &&
    (width > 4 || !state->encoder_control->cfg.rdoq_skip))
  {
    uvg_rdoq(state, coeff, coeff_out, width, height, cur_cu->joint_cb_cr == 1 ? COLOR_V : COLOR_U,
             scan_order, cur_cu->type, cur_cu->cbf, lfnst_idx, 0);
  }
  else if (state->encoder_control->cfg.rdoq_enable && false) {
    uvg_ts_rdoq(state, coeff, coeff_out, width, height, cur_cu->joint_cb_cr == 2 ? COLOR_V : COLOR_U,
      scan_order);
  }
  else {
    uvg_quant(state, coeff, coeff_out, width, height, cur_cu->joint_cb_cr == 1 ? COLOR_V : COLOR_U,
      scan_order, cur_cu->type, cur_cu->tr_idx == MTS_SKIP && false, lfnst_idx);
  }

  int8_t has_coeffs = 0;
  {
    int i;
    for (i = 0; i < width * height; ++i) {
      if (coeff_out[i] != 0) {
        has_coeffs = 1;
        break;
      }
    }
  }

  if (has_coeffs && !early_skip) {

    // Get quantized residual. (coeff_out -> coeff -> residual)
    uvg_dequant(state, coeff_out, coeff, width, height, cur_cu->joint_cb_cr == 1 ? COLOR_V : COLOR_U,
      cur_cu->type, cur_cu->tr_idx == MTS_SKIP && false);
    if (lfnst_idx) {
      uvg_inv_lfnst(cur_cu, width, height, COLOR_UV, lfnst_idx, coeff, tree_type, state->collocated_luma_mode);
    }
    
    uvg_itransform2d(state->encoder_control, combined_residual, coeff, width, height, cur_cu->joint_cb_cr == 1 ? COLOR_V : COLOR_U, cur_cu);
    

    //if (state->tile->frame->lmcs_aps->m_sliceReshapeInfo.enableChromaAdj && color != COLOR_Y) {
    //  int y, x;
    //  int sign, absval;
    //  int maxAbsclipBD = (1 << UVG_BIT_DEPTH) - 1;
    //  for (y = 0; y < width; ++y) {
    //    for (x = 0; x < width; ++x) {
    //      residual[x + y * width] = (int16_t)CLIP((int16_t)(-maxAbsclipBD - 1), (int16_t)maxAbsclipBD, residual[x + y * width]);
    //      sign = residual[x + y * width] >= 0 ? 1 : -1;
    //      absval = sign * residual[x + y * width];
    //      int val = sign * ((absval * lmcs_chroma_adj + (1 << (CSCALE_FP_PREC - 1))) >> CSCALE_FP_PREC);
    //      if (sizeof(uvg_pixel) == 2) // avoid overflow when storing data
    //      {
    //        val = CLIP(-32768, 32767, val);
    //      }
    //      residual[x + y * width] = (int16_t)val;
    //    }
    //  }
    //}
    const int temp = cur_cu->joint_cb_cr * (state->frame->jccr_sign ? -1 : 1);
    // Get quantized reconstruction. (residual + pred_in -> rec_out)
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        if (temp == 2) {
          u_residual[x + y * width] = combined_residual[x + y * width];
          v_residual[x + y * width] = combined_residual[x + y * width] >> 1;
        }
        else if (temp == -2) {
          u_residual[x + y * width] = combined_residual[x + y * width];
          v_residual[x + y * width] = -combined_residual[x + y * width] >> 1;
        }
        else if (temp == 3) {
          u_residual[x + y * width] = combined_residual[x + y * width];
          v_residual[x + y * width] = combined_residual[x + y * width];
        }
        else if (temp == -3) {
          // non-normative clipping to prevent 16-bit overflow
          u_residual[x + y * width] = combined_residual[x + y * width]; // == -32768 && sizeof(Pel) == 2) ? 32767 : -v1_residual[best_cbf_mask][x];
          v_residual[x + y * width] = -combined_residual[x + y * width];
        }
        else if (temp == 1) {
          u_residual[x + y * width] = combined_residual[x + y * width] >> 1;
          v_residual[x + y * width] = combined_residual[x + y * width];
        }
        else if (temp == -1) {
          u_residual[x + y * width] = -combined_residual[x + y * width] >> 1;
          v_residual[x + y * width] = combined_residual[x + y * width];
        }
      }
    }
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        int16_t u_val = u_residual[x + y * width] + u_pred_in[x + y * in_stride];
        u_rec_out[x + y * out_stride] = (uvg_pixel)CLIP(0, PIXEL_MAX, u_val);
        int16_t v_val = v_residual[x + y * width] + v_pred_in[x + y * in_stride];
        v_rec_out[x + y * out_stride] = (uvg_pixel)CLIP(0, PIXEL_MAX, v_val);
      }
    }
  }
  else/* if (rec_out != pred_in)*/ {
    // With no coeffs and rec_out == pred_int we skip copying the coefficients
    // because the reconstruction is just the prediction.

    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        u_rec_out[x + y * out_stride] = u_pred_in[x + y * in_stride];
        v_rec_out[x + y * out_stride] = v_pred_in[x + y * in_stride];
      }
    }
  }
  
  return has_coeffs ? cur_cu->joint_cb_cr : 0;
}

/**
* \brief Quantize residual and get both the reconstruction and coeffs.
*
* \param width  Transform width.
* \param color  Color.
* \param scan_order  Coefficient scan order.
* \param use_trskip  Whether transform skip is used.
* \param stride  Stride for ref_in, pred_in and rec_out.
* \param ref_in  Reference pixels.
* \param pred_in  Predicted pixels.
* \param rec_out  Reconstructed pixels.
* \param coeff_out  Coefficients used for reconstruction of rec_out.
* \param early_skip if this is used for early skip, bypass IT and IQ
*
* \returns  Whether coeff_out contains any non-zero coefficients.
*/
int uvg_quantize_residual_generic(encoder_state_t *const state,
  const cu_info_t *const cur_cu, const int width, const int height, const color_t color,
  const coeff_scan_order_t scan_order, const int use_trskip,
  const int in_stride, const int out_stride,
  const uvg_pixel *const ref_in, const uvg_pixel *const pred_in,
  uvg_pixel *rec_out, coeff_t *coeff_out,
  bool early_skip, int lmcs_chroma_adj, enum uvg_tree_type tree_type)
{
  // Temporary arrays to pass data to and from uvg_quant and transform functions.
  ALIGNED(64) int16_t residual[TR_MAX_WIDTH * TR_MAX_WIDTH];
  ALIGNED(64) coeff_t coeff[TR_MAX_WIDTH * TR_MAX_WIDTH];

  int has_coeffs = 0;

  // With ISP these checks no longer apply, since width and height 2 is now possible
  // With MTT even 1x16 and 16x1 ISP splits are possible
  //assert(width <= TR_MAX_WIDTH && height <= TR_MAX_WIDTH);
  //assert(width >= TR_MIN_WIDTH && height >= TR_MIN_WIDTH);

  // Get residual. (ref_in - pred_in -> residual)
  uvg_generate_residual(ref_in, pred_in, residual, width, height, in_stride, in_stride);

  if (state->tile->frame->lmcs_aps->m_sliceReshapeInfo.enableChromaAdj && color != COLOR_Y) {
    int y, x;
    int sign, absval;
    int maxAbsclipBD = (1 << UVG_BIT_DEPTH) - 1;
    for (y = 0; y < height; ++y) {
      for (x = 0; x < width; ++x) {
        sign = residual[x + y * width] >= 0 ? 1 : -1;
        absval = sign * residual[x + y * width];
        residual[x + y * width] = (int16_t)CLIP(-maxAbsclipBD, maxAbsclipBD, sign * (((absval << CSCALE_FP_PREC) + (lmcs_chroma_adj >> 1)) / lmcs_chroma_adj));
      }
    }
  }

  // Transform residual. (residual -> coeff)
  if (use_trskip) {
    uvg_transformskip(state->encoder_control, residual, coeff, width, height);
  }
  else {
    uvg_transform2d(state->encoder_control, residual, coeff, width, height, color, cur_cu);
  }

  const uint8_t lfnst_index = tree_type != UVG_CHROMA_T || color == COLOR_Y ? cur_cu->lfnst_idx : cur_cu->cr_lfnst_idx;

  if (state->encoder_control->cfg.lfnst && cur_cu->type == CU_INTRA) {
    // Forward low frequency non-separable transform
    uvg_fwd_lfnst(cur_cu, width, height, color, lfnst_index, coeff, tree_type, state->collocated_luma_mode);
  }
  

  // Quantize coeffs. (coeff -> coeff_out)
  
  int abs_sum = 0;
  if (!use_trskip && state->encoder_control->cfg.dep_quant) {
    uvg_dep_quant(
      state,
      cur_cu,
      width,
      height,
      coeff,
      coeff_out,
      color,
      tree_type,
      &abs_sum,
      state->encoder_control->cfg.scaling_list);
  }
  else if (state->encoder_control->cfg.rdoq_enable &&
      (width > 4 || !state->encoder_control->cfg.rdoq_skip) && !use_trskip)
  {
    uvg_rdoq(state, coeff, coeff_out, width, height, color,
             scan_order, cur_cu->type, cur_cu->cbf, lfnst_index, color == 0 ? cur_cu->tr_idx : 0);
  } else if(state->encoder_control->cfg.rdoq_enable && use_trskip) {
    uvg_ts_rdoq(state, coeff, coeff_out, width, height, color,
      scan_order);
  } else {
  
    uvg_quant(state, coeff, coeff_out, width, height, color,
      scan_order, cur_cu->type, cur_cu->tr_idx == MTS_SKIP && color == COLOR_Y, lfnst_index);
  }

  // Check if there are any non-zero coefficients.
  {
    int i;
    for (i = 0; i < width * height; ++i) {
      if (coeff_out[i] != 0) {
        has_coeffs = 1;
        break;
      }
    }
  }

  // Do the inverse quantization and transformation and the reconstruction to
  // rec_out.
  if (has_coeffs && !early_skip) {
    int y, x;

    // Get quantized residual. (coeff_out -> coeff -> residual)
    uvg_dequant(state, coeff_out, coeff, width, height, color,
      cur_cu->type, cur_cu->tr_idx == MTS_SKIP && color == COLOR_Y);
    
    if (state->encoder_control->cfg.lfnst && cur_cu->type == CU_INTRA) {
      // Inverse low frequency non-separable transform
      uvg_inv_lfnst(cur_cu, width, height, color, lfnst_index, coeff, tree_type, state->collocated_luma_mode);
    }
    if (use_trskip) {
      uvg_itransformskip(state->encoder_control, residual, coeff, width, height);
    }
    else {
      uvg_itransform2d(state->encoder_control, residual, coeff, width, height, color, cur_cu);
    }
    
    if (state->tile->frame->lmcs_aps->m_sliceReshapeInfo.enableChromaAdj && color != COLOR_Y) {
      int y, x;
      int sign, absval;
      int maxAbsclipBD = (1 << UVG_BIT_DEPTH) - 1;
      for (y = 0; y < height; ++y) {
        for (x = 0; x < width; ++x) {
          residual[x + y * width] = (int16_t)CLIP((int16_t)(-maxAbsclipBD - 1), (int16_t)maxAbsclipBD, residual[x + y * width]);
          sign = residual[x + y * width] >= 0 ? 1 : -1;
          absval = sign * residual[x + y * width];
          int val = sign * ((absval * lmcs_chroma_adj + (1 << (CSCALE_FP_PREC - 1))) >> CSCALE_FP_PREC);
          if (sizeof(uvg_pixel) == 2) // avoid overflow when storing data
          {
            val = CLIP(-32768, 32767, val);
          }
          residual[x + y * width] = (int16_t)val;
        }
      }
    }

    // Get quantized reconstruction. (residual + pred_in -> rec_out)
    for (y = 0; y < height; ++y) {
      for (x = 0; x < width; ++x) {
        int16_t val = residual[x + y * width] + pred_in[x + y * in_stride];
        rec_out[x + y * out_stride] = (uvg_pixel)CLIP(0, PIXEL_MAX, val);
      }
    }
  }
  else if (rec_out != pred_in) {
    // With no coeffs and rec_out == pred_int we skip copying the coefficients
    // because the reconstruction is just the prediction.
    int y, x;

    for (y = 0; y < height; ++y) {
      for (x = 0; x < width; ++x) {
        rec_out[x + y * out_stride] = pred_in[x + y * in_stride];
      }
    }
  }

  return has_coeffs;
}

/**
 * \brief inverse quantize transformed and quantized coefficents
 *
 */
void uvg_dequant_generic(const encoder_state_t * const state, coeff_t *q_coef, coeff_t *coef, int32_t width, int32_t height,color_t color, int8_t block_type, int8_t transform_skip)
{
  const encoder_control_t * const encoder = state->encoder_control;
  if(encoder->cfg.dep_quant && !transform_skip) {
    uvg_dep_quant_dequant(state, block_type, width, height, color, q_coef, coef, encoder->cfg.scaling_list);
    return;
  }
  int32_t shift,add,coeff_q;
  int32_t n;
  const uint32_t log2_tr_width  = uvg_g_convert_to_log2[width];
  const uint32_t log2_tr_height = uvg_g_convert_to_log2[height];
  int32_t transform_shift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - ((log2_tr_width + log2_tr_height) >> 1); // Represents scaling through forward transform

  bool needs_block_size_trafo_scale = !transform_skip && ((log2_tr_height + log2_tr_width) % 2 == 1);
  needs_block_size_trafo_scale |= 0; // Non log2 block size

  int32_t qp_scaled = uvg_get_scaled_qp(color, state->qp, (encoder->bitdepth-8)*6, encoder->qp_map[0]);
  qp_scaled = transform_skip ? MAX(qp_scaled, 4 + 6 * MIN_QP_PRIME_TS) : qp_scaled;

  shift = 20 - QUANT_SHIFT - (transform_skip ? 0 : transform_shift - needs_block_size_trafo_scale);

  if (encoder->scaling_list.enable)
  {
    int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)(color);

    const int32_t *dequant_coef = encoder->scaling_list.de_quant_coeff[log2_tr_width][log2_tr_height][scalinglist_type][qp_scaled%6];
    shift += 4;

    if (shift >qp_scaled / 6) {
      add = 1 << (shift - qp_scaled/6 - 1);

      for (n = 0; n < width * height; n++) {
        coeff_q = ((q_coef[n] * dequant_coef[n]) + add ) >> (shift -  qp_scaled/6);
        coef[n] = (coeff_t)CLIP(-32768,32767,coeff_q);
      }
    } else {
      for (n = 0; n < width * height; n++) {
        // Clip to avoid possible overflow in following shift left operation
        coeff_q   = CLIP(-32768, 32767, q_coef[n] * dequant_coef[n]);
        coef[n] = (coeff_t)CLIP(-32768, 32767, coeff_q << (qp_scaled/6 - shift));
      }
    }
  } else {
    int32_t scale = uvg_g_inv_quant_scales[needs_block_size_trafo_scale][qp_scaled%6] << (qp_scaled/6);
    add = 1 << (shift-1);

    for (n = 0; n < width * height; n++) {
      coeff_q   = (q_coef[n] * scale + add) >> shift;
      coef[n] = (coeff_t)CLIP(-32768, 32767, coeff_q);
    }
  }
}

static uint32_t coeff_abs_sum_generic(const coeff_t *coeffs, size_t length)
{
  uint32_t sum = 0;
  for (int i = 0; i < length; i++) {
    sum += abs(coeffs[i]);
  }
  return sum;
}

static INLINE void get_coeff_weights(uint64_t wts_packed, uint16_t *weights)
{
  weights[0] = (wts_packed >>  0) & 0xffff;
  weights[1] = (wts_packed >> 16) & 0xffff;
  weights[2] = (wts_packed >> 32) & 0xffff;
  weights[3] = (wts_packed >> 48) & 0xffff;
}

static uint32_t fast_coeff_cost_generic(const coeff_t *coeff, int32_t width, int32_t height, uint64_t weights)
{
  assert((width == height) && "Non-square block handling not implemented for this function.");
  uint32_t sum = 0;
  uint16_t weights_unpacked[4];

  get_coeff_weights(weights, weights_unpacked);

  for (int32_t i = 0; i < width * height; i++) {
     int16_t curr = coeff[i];
    uint32_t curr_abs = abs(curr);
    if (curr_abs > 3) {
      curr_abs = 3;
    }
    sum += weights_unpacked[curr_abs];
  }
  return (sum + (1 << 7)) >> 8;
}

int uvg_strategy_register_quant_generic(void* opaque, uint8_t bitdepth)
{
  bool success = true;

  success &= uvg_strategyselector_register(opaque, "quant", "generic", 0, &uvg_quant_generic);
  success &= uvg_strategyselector_register(opaque, "quant_cbcr_residual", "generic", 0, &uvg_quant_cbcr_residual_generic);
  success &= uvg_strategyselector_register(opaque, "quantize_residual", "generic", 0, &uvg_quantize_residual_generic);
  success &= uvg_strategyselector_register(opaque, "dequant", "generic", 0, &uvg_dequant_generic);
  success &= uvg_strategyselector_register(opaque, "coeff_abs_sum", "generic", 0, &coeff_abs_sum_generic);
  success &= uvg_strategyselector_register(opaque, "fast_coeff_cost", "generic", 0, &fast_coeff_cost_generic);

  return success;
}
