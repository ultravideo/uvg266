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

#include "strategyselector.h"

#include "cabac.h"
#include "context.h"
#include "encode_coding_tree-generic.h"
#include "encode_coding_tree.h"


 /**
 * \brief Encode block coefficients
 *
 * \param state     current encoder state
 * \param cabac     current cabac state
 * \param coeff     Input coefficients
 * \param width     Block width
 * \param color     plane type / luminance or chrominance
 * \param scan_mode    scan type (diag, hor, ver) DEPRECATED?
 *
 * This method encodes coefficients of a block
 *
 */
void uvg_encode_coeff_nxn_generic(encoder_state_t * const state,
  cabac_data_t * const cabac,
  const coeff_t *coeff,
  const cu_loc_t * const cu_loc,
  uint8_t color,
  int8_t scan_mode,
  cu_info_t* cur_cu,
  double* bits_out) 
{
  const int x = cu_loc->x;
  const int y = cu_loc->y;
  const int width  = color == COLOR_Y ? cu_loc->width  : cu_loc->chroma_width;
  const int height = color == COLOR_Y ? cu_loc->height : cu_loc->chroma_height;

  //const encoder_control_t * const encoder = state->encoder_control;
  //int c1 = 1;
  uint8_t last_coeff_x = 0;
  uint8_t last_coeff_y = 0;
  int32_t i;
  // ToDo: large block support in VVC?
  uint32_t sig_coeffgroup_flag[32 * 32] = { 0 };

  int32_t scan_pos;
  //int32_t next_sig_pos;
  uint32_t blk_pos, pos_y, pos_x, sig, ctx_sig;
  double bits = 0;

  // CONSTANTS

  const uint8_t log2_block_width =  uvg_g_convert_to_log2[width];
  const uint8_t log2_block_height = uvg_g_convert_to_log2[height];
  
  const uint32_t log2_cg_size = uvg_g_log2_sbb_size[log2_block_width][log2_block_height][0] + uvg_g_log2_sbb_size[log2_block_width][log2_block_height][1];
  const uint32_t* const scan = uvg_get_scan_order_table(SCAN_GROUP_4X4, scan_mode, log2_block_width, log2_block_height, 0);
  const uint32_t* const scan_cg = uvg_get_scan_order_table(SCAN_GROUP_UNGROUPED, scan_mode, log2_block_width, log2_block_height, 0);


  // Init base contexts according to block type
  cabac_ctx_t *base_coeff_group_ctx = &(cabac->ctx.sig_coeff_group_model[(color == 0 ? 0 : 1) * 2]);
  

  unsigned scan_cg_last = (unsigned)-1;
  unsigned scan_pos_last = (unsigned)-1;

  for (int i = 0; i < (width * height); ++i) {
    if (coeff[scan[i]]) {
      scan_pos_last = i;
      sig_coeffgroup_flag[scan_cg[i >> log2_cg_size]] = 1;
    }
  }

  scan_cg_last = scan_pos_last >> log2_cg_size;

  int pos_last = scan[scan_pos_last];

  last_coeff_y = (uint8_t)(pos_last / width);
  last_coeff_x = (uint8_t)(pos_last - (last_coeff_y * width));
  bool is_chroma = color != COLOR_Y;

  if (cur_cu != NULL && /*cur_cu->tr_idx != MTS_SKIP &&*/ height >= 4 && width >= 4) {
    const unsigned max_lfnst_pos = ((height == 4 && width == 4) || (height == 8 && width == 8)) ? 7 : 15;
    if(!is_chroma) {
      cur_cu->violates_lfnst_constrained_luma |= scan_pos_last > max_lfnst_pos;
    }
    else {
      cur_cu->violates_lfnst_constrained_chroma |= scan_pos_last > max_lfnst_pos;
    }
    cur_cu->lfnst_last_scan_pos |= scan_pos_last >= 1;
  }

  // Code last_coeff_x and last_coeff_y
  uvg_encode_last_significant_xy(cabac,
    last_coeff_x,
    last_coeff_y,
    width,
    height,
    color,
    scan_mode,
    bits_out);



  uint32_t quant_state_transition_table = state->encoder_control->cfg.dep_quant ? 32040 : 0; 
  int32_t quant_state = 0;
  uint8_t  ctx_offset[16];
  int32_t temp_diag = -1;
  int32_t temp_sum = -1;

  int32_t reg_bins = (width * height * 28) >> 4; //8 for 2x2

  // significant_coeff_flag
  for (i = scan_cg_last; i >= 0; i--) {

    //int32_t abs_coeff[64*64];
    const uint32_t log2_cg_width = uvg_g_log2_sbb_size[log2_block_width][log2_block_height][0];
    const uint32_t log2_cg_height = uvg_g_log2_sbb_size[log2_block_width][log2_block_height][1];
    const uint32_t cg_width = (MIN((uint8_t)TR_MAX_WIDTH, width) >> log2_cg_width);
    const uint32_t cg_height = (MIN((uint8_t)TR_MAX_WIDTH, height) >> log2_cg_height);
    int32_t cg_blk_pos = scan_cg[i];
    int32_t cg_pos_y = cg_blk_pos / (MIN((uint8_t)32, width) >> log2_cg_width);
    int32_t cg_pos_x = cg_blk_pos - (cg_pos_y * (MIN((uint8_t)32, width) >> log2_cg_width));


    // !!! residual_coding_subblock() !!!

    // Encode significant coeff group flag when not the last or the first
    if (i == scan_cg_last || i == 0) {
      sig_coeffgroup_flag[cg_blk_pos] = 1;
    } else {
      uint32_t sig_coeff_group = (sig_coeffgroup_flag[cg_blk_pos] != 0);
      uint32_t ctx_sig = uvg_context_get_sig_coeff_group(sig_coeffgroup_flag, cg_pos_x,
        cg_pos_y, cg_width, cg_height);
      CABAC_FBITS_UPDATE(cabac, &base_coeff_group_ctx[ctx_sig], sig_coeff_group, bits, "significant_coeffgroup_flag");
    }


    if (sig_coeffgroup_flag[cg_blk_pos]) {

      int32_t min_sub_pos = i << log2_cg_size; // LOG2_SCAN_SET_SIZE;
      int32_t first_sig_pos = (i == scan_cg_last) ? scan_pos_last : (min_sub_pos + (1 << log2_cg_size) - 1);
      int32_t next_sig_pos = first_sig_pos;

      int32_t infer_sig_pos = (next_sig_pos != scan_pos_last) ? ((i != 0) ? min_sub_pos : -1) : next_sig_pos;
      int32_t num_non_zero = 0;
      int32_t last_nz_pos_in_cg = -1;
      int32_t first_nz_pos_in_cg = next_sig_pos;
      int32_t remainder_abs_coeff = -1;
      uint32_t coeff_signs = 0;


      /*
         ****  FIRST PASS ****
      */
      for (next_sig_pos = first_sig_pos; next_sig_pos >= min_sub_pos && reg_bins >= 4; next_sig_pos--) {


        blk_pos = scan[next_sig_pos];
        pos_y = blk_pos / width;
        pos_x = blk_pos - (pos_y * width);

        sig = (coeff[blk_pos] != 0) ? 1 : 0;
        if (num_non_zero || next_sig_pos != infer_sig_pos) {
          ctx_sig = uvg_context_get_sig_ctx_idx_abs(coeff, pos_x, pos_y, width, height, color, &temp_diag, &temp_sum);
          cabac_ctx_t* sig_ctx_luma = &(cabac->ctx.cu_sig_model_luma[MAX(0, (quant_state - 1))][ctx_sig]);
          cabac_ctx_t* sig_ctx_chroma = &(cabac->ctx.cu_sig_model_chroma[MAX(0, (quant_state - 1))][MIN(ctx_sig,7)]);

          CABAC_FBITS_UPDATE(cabac, (color == 0 ? sig_ctx_luma : sig_ctx_chroma), sig, bits, "sig_coeff_flag");
          reg_bins--;

        } else if (next_sig_pos != scan_pos_last) {
          ctx_sig = uvg_context_get_sig_ctx_idx_abs(coeff, pos_x, pos_y, width, height, color, &temp_diag, &temp_sum);
        }


        if (sig) {
          assert(next_sig_pos - min_sub_pos >= 0 && next_sig_pos - min_sub_pos < 16);
          uint8_t* offset = &ctx_offset[next_sig_pos - min_sub_pos];
          num_non_zero++;
          // ctxOffsetAbs()
          {
            *offset = 0;
            if (temp_diag != -1) {
              *offset = MIN(temp_sum, 4) + 1;
              *offset += (!temp_diag ? (color == COLOR_Y ? 15 : 5) : color == COLOR_Y ? temp_diag < 3 ? 10 : (temp_diag < 10 ? 5 : 0) : 0);
            }
          }


          last_nz_pos_in_cg = MAX(last_nz_pos_in_cg, next_sig_pos);
          first_nz_pos_in_cg = next_sig_pos;

          remainder_abs_coeff = abs(coeff[blk_pos]) - 1;

          // If shift sign pattern and add current sign
          coeff_signs = (next_sig_pos != scan_pos_last ? 2 * coeff_signs : coeff_signs) + (coeff[blk_pos] < 0);



          // Code "greater than 1" flag
          uint8_t gt1 = remainder_abs_coeff ? 1 : 0;
          CABAC_FBITS_UPDATE(cabac, (color == 0) ? &(cabac->ctx.cu_gtx_flag_model_luma[1][*offset]) :
            &(cabac->ctx.cu_gtx_flag_model_chroma[1][*offset]),
            gt1, bits, "abs_level_gtx_flag");
          reg_bins--;

          if (gt1) {
            remainder_abs_coeff -= 1;

            // Code coeff parity
            CABAC_FBITS_UPDATE(cabac, (color == 0) ? &(cabac->ctx.cu_parity_flag_model_luma[*offset]) :
              &(cabac->ctx.cu_parity_flag_model_chroma[*offset]),
              remainder_abs_coeff & 1, bits, "par_flag");
            remainder_abs_coeff >>= 1;

            reg_bins--;
            uint8_t gt2 = remainder_abs_coeff ? 1 : 0;
            CABAC_FBITS_UPDATE(cabac, (color == 0) ? &(cabac->ctx.cu_gtx_flag_model_luma[0][*offset]) :
              &(cabac->ctx.cu_gtx_flag_model_chroma[0][*offset]),
              gt2, bits, "gt2_flag");
            reg_bins--;
          }
        }

        quant_state = (quant_state_transition_table >> ((quant_state << 2) + ((coeff[blk_pos] & 1) << 1))) & 3;
      }


      /*
      ****  SECOND PASS: Go-rice  ****
      */
      uint32_t rice_param = 0;
      uint32_t pos0 = 0;
      for (scan_pos = first_sig_pos; scan_pos > next_sig_pos; scan_pos--) {
        blk_pos = scan[scan_pos];
        pos_y = blk_pos / width;
        pos_x = blk_pos - (pos_y * width);
        int32_t abs_sum = uvg_abs_sum(coeff, pos_x, pos_y, width, height, 4);

        rice_param = g_go_rice_pars[abs_sum];
        uint32_t second_pass_abs_coeff = abs(coeff[blk_pos]);
        if (second_pass_abs_coeff >= 4) {
          uint32_t remainder = (second_pass_abs_coeff - 4) >> 1;
          bits += uvg_cabac_write_coeff_remain(cabac, remainder, rice_param, 5);
        }
      }

      /*
      ****  coeff bypass  ****
      */
      for (scan_pos = next_sig_pos; scan_pos >= min_sub_pos; scan_pos--) {
        blk_pos = scan[scan_pos];
        pos_y = blk_pos / width;
        pos_x = blk_pos - (pos_y * width);
        uint32_t coeff_abs = abs(coeff[blk_pos]);
        int32_t abs_sum = uvg_abs_sum(coeff, pos_x, pos_y, width, height, 0);
        rice_param = g_go_rice_pars[abs_sum];        
        pos0 = ((quant_state<2)?1:2) << rice_param;
        uint32_t remainder = (coeff_abs == 0 ? pos0 : coeff_abs <= pos0 ? coeff_abs - 1 : coeff_abs);
        bits += uvg_cabac_write_coeff_remain(cabac, remainder, rice_param, 5);
        quant_state = (quant_state_transition_table >> ((quant_state << 2) + ((coeff_abs & 1) << 1))) & 3;
        if (coeff_abs) {
          num_non_zero++;
          first_nz_pos_in_cg = scan_pos;
          last_nz_pos_in_cg = MAX(last_nz_pos_in_cg, scan_pos);
          coeff_signs <<= 1;
          if (coeff[blk_pos] < 0) coeff_signs++;
        }
      }

      uint32_t num_signs = num_non_zero;

      if (state->encoder_control->cfg.signhide_enable && !state->encoder_control->cfg.dep_quant && (last_nz_pos_in_cg - first_nz_pos_in_cg >= 4)) {
        num_signs--;
        coeff_signs >>= 1;
      }

      if (color == COLOR_Y && cur_cu != NULL && cur_cu->tr_idx != MTS_SKIP)
      {
        cur_cu->mts_last_scan_pos |= first_sig_pos > 0;
      }

      CABAC_BINS_EP(cabac, coeff_signs, num_signs, "coeff_signs");
      if (cabac->only_count) bits += num_signs;
    }

    if (color == COLOR_Y && cur_cu != NULL && (cg_pos_y > 3 || cg_pos_x > 3) && sig_coeffgroup_flag[cg_blk_pos] != 0)
    {
      cur_cu->violates_mts_coeff_constraint = true;
    }
  }
  if (cabac->only_count && bits_out) *bits_out += bits;
}


int uvg_strategy_register_encode_generic(void* opaque, uint8_t bitdepth)
{
  bool success = true;

  success &= uvg_strategyselector_register(opaque, "encode_coeff_nxn", "generic", 0, &uvg_encode_coeff_nxn_generic);

  return success;
}
