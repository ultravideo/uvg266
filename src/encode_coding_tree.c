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

#include "encode_coding_tree.h"

#include "cabac.h"
#include "context.h"
#include "cu.h"
#include "debug.h"
#include "encoder.h"
#include "global.h"
#include "imagelist.h"
#include "inter.h"
#include "intra.h"
#include "uvg266.h"
#include "uvg_math.h"
#include "strategyselector.h"
#include "tables.h"
#include "videoframe.h"

bool uvg_is_mts_allowed(const encoder_state_t * const state, cu_info_t *const pred_cu, const cu_loc_t*
                        const cu_loc)
{
  uint32_t ts_max_size = 1 << state->encoder_control->cfg.trskip_max_size; 
  const uint32_t max_size = 32; // CU::isIntra(cu) ? MTS_INTRA_MAX_CU_SIZE : MTS_INTER_MAX_CU_SIZE;
  const uint32_t cu_width    = 1 << pred_cu->log2_width;
  const uint32_t cu_height   = 1 << pred_cu->log2_height;
  //bool mts_allowed = cu.chType == CHANNEL_TYPE_LUMA && compID == COMPONENT_Y;

  uint8_t mts_type = state->encoder_control->cfg.mts;
  bool mts_allowed = mts_type == UVG_MTS_BOTH || (pred_cu->type == CU_INTRA ? mts_type == UVG_MTS_INTRA : pred_cu->type == CU_INTER && mts_type == UVG_MTS_INTER);
  mts_allowed &= cu_width <= max_size && cu_height <= max_size;
  mts_allowed &= pred_cu->type == CU_INTRA ? !pred_cu->intra.isp_mode : true;
  //mts_allowed &= !cu.sbtInfo;
  mts_allowed &= !(pred_cu->bdpcmMode && cu_width <= ts_max_size && cu_height <= ts_max_size);
  mts_allowed &= pred_cu->tr_idx != MTS_SKIP && !pred_cu->violates_mts_coeff_constraint && pred_cu->mts_last_scan_pos ;
  mts_allowed &= pred_cu->lfnst_idx == 0;
  return mts_allowed;
}

static void encode_mts_idx(
  encoder_state_t * const state,
  cabac_data_t * const cabac,
  const cu_info_t *const pred_cu,
  const cu_loc_t* const cu_loc)
{
  //TransformUnit &tu = *cu.firstTU;
  int mts_idx = pred_cu->tr_idx;

  if (uvg_is_mts_allowed(state, (cu_info_t* const )pred_cu, cu_loc) && mts_idx != MTS_SKIP
       && !pred_cu->violates_mts_coeff_constraint
       && pred_cu->mts_last_scan_pos       
    )
  {
    int symbol = mts_idx != MTS_DCT2_DCT2 ? 1 : 0;
    int ctx_idx = 0;

    cabac->cur_ctx = &(cabac->ctx.mts_idx_model[ctx_idx]);
    CABAC_BIN(cabac, symbol, "mts_idx");

    if (symbol)
    {
      ctx_idx = 1;
      for (int i = 0; i < 3; i++, ctx_idx++)
      {
        symbol = mts_idx > i + MTS_DST7_DST7 ? 1 : 0;
        cabac->cur_ctx = &(cabac->ctx.mts_idx_model[ctx_idx]);
        CABAC_BIN(cabac, symbol, "mts_idx");

        if (!symbol)
        {
          break;
        }
      }
    }
  }
}


bool uvg_is_lfnst_allowed(
  const encoder_state_t* const state,
  const cu_info_t* const pred_cu,
  enum uvg_tree_type tree_type,
  const color_t color,
  const cu_loc_t* const cu_loc, const lcu_t* const lcu) 
{
  if (state->encoder_control->cfg.lfnst && pred_cu->type == CU_INTRA && PU_IS_TU(pred_cu)) {
    const int isp_mode = pred_cu->intra.isp_mode;
    const int cu_width  = tree_type != UVG_CHROMA_T ? 1 << pred_cu->log2_width : 1 << pred_cu->log2_chroma_width;
    const int cu_height = tree_type != UVG_CHROMA_T ? 1 << pred_cu->log2_height : 1 << pred_cu->log2_chroma_height;
    bool can_use_lfnst_with_mip = (cu_width >= 16 && cu_height >= 16);
    bool is_sep_tree = tree_type != UVG_BOTH_T;
    bool mip_flag = pred_cu->type == CU_INTRA && color == COLOR_Y ? pred_cu->intra.mip_flag : false;

    if ((isp_mode && !uvg_can_use_isp_with_lfnst(cu_width, cu_height, isp_mode, tree_type) && color == COLOR_Y) ||
      (pred_cu->type == CU_INTRA && mip_flag && !can_use_lfnst_with_mip && color == COLOR_Y) ||
      (is_sep_tree && MIN(cu_width, cu_height) < 4) || 
      (cu_width > (TR_MAX_WIDTH >> (tree_type == UVG_CHROMA_T)) || cu_height > (TR_MAX_WIDTH >> (tree_type == UVG_CHROMA_T)))) {
      return false;
    }
    bool luma_flag = tree_type != UVG_CHROMA_T;
    bool chroma_flag = tree_type != UVG_LUMA_T;
    bool non_zero_coeff_non_ts_corner_8x8 = false;
    bool last_scan_pos = false;
    bool is_tr_skip = false;
    
    int split_num = color == COLOR_Y && isp_mode ? uvg_get_isp_split_num(cu_width, cu_height, isp_mode, false) : 0;
    const videoframe_t* const frame = state->tile->frame;

    if (split_num) {
      // Constraints for ISP split blocks
      for (int i = 0; i < split_num; ++i) {
        cu_loc_t split_loc;
        uvg_get_isp_split_loc(&split_loc, cu_loc->x, cu_loc->y, cu_width, cu_height, i, isp_mode, false);
        int local_split_x = lcu ? split_loc.local_x : split_loc.x;
        int local_split_y = lcu ? split_loc.local_y : split_loc.y;
        uvg_get_isp_cu_arr_coords(&local_split_x, &local_split_y, MAX(cu_width, cu_height));
        const cu_info_t* split_cu = lcu ? LCU_GET_CU_AT_PX(lcu, local_split_x, local_split_y) :
          uvg_cu_array_at_const(frame->cu_array, local_split_x, local_split_y);

        //if (cbf_is_set(split_cu->cbf, depth, COLOR_Y)) {
        // ISP_TODO: remove this if clause altogether if it seems it is not needed
        if (true) {
          non_zero_coeff_non_ts_corner_8x8 |= (luma_flag && split_cu->violates_lfnst_constrained_luma) || (chroma_flag && split_cu->violates_lfnst_constrained_chroma);
          //last_scan_pos |= split_cu->lfnst_last_scan_pos;
          last_scan_pos |= true;
        }
      }
    }
    else {
      non_zero_coeff_non_ts_corner_8x8 |= (luma_flag && pred_cu->violates_lfnst_constrained_luma) || (chroma_flag && pred_cu->violates_lfnst_constrained_chroma);
      last_scan_pos |= pred_cu->lfnst_last_scan_pos;
    }

    if (color == COLOR_Y && pred_cu->tr_idx == MTS_SKIP) {
      is_tr_skip = true;
    }

    if ((!last_scan_pos) || non_zero_coeff_non_ts_corner_8x8 || is_tr_skip) {
      return false;
    }
    return true;
  }
  else {
    return false;
  }
}

static bool encode_lfnst_idx(
  const encoder_state_t* const state,
  cabac_data_t * const cabac,
  const cu_info_t * const pred_cu,
  enum uvg_tree_type tree_type,
  const color_t color,
  const cu_loc_t* const cu_loc)
{
  
  if (uvg_is_lfnst_allowed(state, pred_cu, tree_type, color, cu_loc, NULL)) {
    // Getting separate tree bool from block size is a temporary fix until a proper dual tree check is possible (there is no dual tree structure at time of writing this).
    // VTM seems to force explicit dual tree structure for small 4x4 blocks
    bool is_separate_tree = tree_type != UVG_BOTH_T;

    const int lfnst_index = !is_separate_tree || color == COLOR_Y ? pred_cu->lfnst_idx : pred_cu->cr_lfnst_idx;
    assert((lfnst_index >= 0 && lfnst_index < 3) && "Invalid LFNST index.");

    uint16_t ctx_idx = 0;
    if (is_separate_tree) ctx_idx++;

    cabac->cur_ctx = &(cabac->ctx.lfnst_idx_model[ctx_idx]);
    CABAC_BIN(cabac, lfnst_index ? 1 : 0, "lfnst_idx");

    if (lfnst_index) {
      cabac->cur_ctx = &(cabac->ctx.lfnst_idx_model[2]);
      CABAC_BIN(cabac, (lfnst_index - 1) ? 1 : 0, "lfnst_idx");
    }
    return true;
  }
  else {
    if(color == COLOR_Y) {
      assert(pred_cu->lfnst_idx == 0);
    }
    if(tree_type == UVG_CHROMA_T && color != COLOR_Y) {
      assert(pred_cu->cr_lfnst_idx == 0);
    }
    return false;
  }
}

void uvg_encode_ts_residual(encoder_state_t* const state,
  cabac_data_t* const cabac,
  const coeff_t* coeff,
  uint32_t width,
  uint32_t height,
  uint8_t type,
  int8_t scan_mode,
  double* bits_out) 
{
  //const encoder_control_t * const encoder = state->encoder_control;
  //int c1 = 1;
  uint32_t i;
  int32_t blk_pos;
  // ToDo: large block support in VVC?
  uint32_t sig_coeffgroup_flag[32 * 32] = { 0 };
  

  // CONSTANTS

  const uint32_t log2_block_width  = uvg_g_convert_to_log2[width];
  const uint32_t log2_block_height = uvg_g_convert_to_log2[height];
  // TODO: log2_cg_size is wrong if width != height
  const uint32_t log2_cg_size = uvg_g_log2_sbb_size[log2_block_width][log2_block_width][0] + uvg_g_log2_sbb_size[log2_block_width][log2_block_height][1];
  
  const uint32_t* const scan = uvg_get_scan_order_table(SCAN_GROUP_4X4, scan_mode, log2_block_width, log2_block_height);
  const uint32_t* const scan_cg = uvg_get_scan_order_table(SCAN_GROUP_UNGROUPED, scan_mode, log2_block_width, log2_block_height);

  double bits = 0;

  // Init base contexts according to block type
  cabac_ctx_t* base_coeff_group_ctx = &(cabac->ctx.transform_skip_sig_coeff_group[0]);

  cabac->cur_ctx = base_coeff_group_ctx;
  
  int maxCtxBins = (width * height * 7) >> 2;
  unsigned scan_cg_last = (unsigned )-1;
  //unsigned scan_pos_last = (unsigned )-1;

  for (i = 0; i < width * height; i++) {
    if (coeff[scan[i]]) {
      sig_coeffgroup_flag[scan_cg[i >> log2_cg_size]] = 1;
    }
  }
  // TODO: this won't work with non-square blocks
  scan_cg_last = (width * height - 1) >> log2_cg_size;
  const uint32_t cg_width = (MIN((uint8_t)32, width) >> (log2_cg_size / 2));

  bool no_sig_group_before_last = true;

  for (i = 0; i <= scan_cg_last; i++) {
    if (!((width == 4 && height == 4) || (i ==scan_cg_last && no_sig_group_before_last))) {
      uint32_t cg_blkpos = scan_cg[i];
      uint32_t cg_pos_y = cg_blkpos / cg_width;
      uint32_t cg_pos_x = cg_blkpos - (cg_pos_y * cg_width);

      uint32_t ctx_sig = uvg_context_get_sig_coeff_group_ts(sig_coeffgroup_flag, cg_pos_x, cg_pos_y, cg_width);

      cabac->cur_ctx = &base_coeff_group_ctx[ctx_sig];

      if(!sig_coeffgroup_flag[scan_cg[i]]) {
        CABAC_FBITS_UPDATE(cabac, &base_coeff_group_ctx[ctx_sig], 0, bits, "ts_sigGroup");
        continue;
      }
      CABAC_FBITS_UPDATE(cabac, &base_coeff_group_ctx[ctx_sig], 1, bits, "ts_sigGroup");
      no_sig_group_before_last = false;
    }

    int firstSigPos = i << log2_cg_size;
    int min_sub_pos = firstSigPos + (1 << log2_cg_size) - 1;
    int nextSigPos = firstSigPos;

    //===== encode absolute values =====
    const int inferSigPos = min_sub_pos;
    int       remAbsLevel = -1;
    int       numNonZero = 0;

    int rightPixel, belowPixel, modAbsCoeff;

    int lastScanPosPass1 = -1;
    int lastScanPosPass2 = -1;
    for (; nextSigPos <= min_sub_pos && maxCtxBins >= 4; nextSigPos++)
    {
      blk_pos = scan[nextSigPos];
      int pos_x = blk_pos % width;
      int pos_y = blk_pos / width;
      coeff_t curr_coeff = coeff[blk_pos];
      unsigned  sigFlag = (curr_coeff != 0);
      if (numNonZero || nextSigPos != inferSigPos)
      {
        CABAC_FBITS_UPDATE(cabac, &cabac->ctx.transform_skip_sig[
          uvg_context_get_sig_ctx_idx_abs_ts(coeff, pos_x, pos_y, width)
        ], sigFlag, bits, "sig_coeff_flag");
        maxCtxBins--;
      }

      if (sigFlag)
      {
        //===== encode sign's =====
        int sign = curr_coeff < 0;
        CABAC_FBITS_UPDATE(cabac, &cabac->ctx.transform_skip_res_sign[
          uvg_sign_ctx_id_abs_ts(coeff, pos_x, pos_y, width, 0)
        ], sign, bits, "coeff_sign_flag");
        maxCtxBins--;
        numNonZero++;

        rightPixel = pos_x > 0 ? coeff[pos_x + pos_y * width - 1] : 0;
        belowPixel = pos_y > 0 ? coeff[pos_x + (pos_y - 1) * width] : 0;

        modAbsCoeff = uvg_derive_mod_coeff(rightPixel, belowPixel, abs(curr_coeff), 0);

        remAbsLevel = modAbsCoeff - 1;

        unsigned gt1 = !!remAbsLevel;
        CABAC_FBITS_UPDATE(cabac, &cabac->ctx.transform_skip_gt1[
          uvg_lrg1_ctx_id_abs_ts(coeff, pos_x, pos_y, width, 0)
        ], gt1, bits, "abs_level_gtx_flag");
        maxCtxBins--;

        if (gt1)
        {
          remAbsLevel -= 1;
          CABAC_FBITS_UPDATE(cabac, &cabac->ctx.transform_skip_par, remAbsLevel & 1, bits, "par_level_flag");
          maxCtxBins--;
        }
      }
      lastScanPosPass1 = nextSigPos;
    }


    uint32_t cutoffVal = 2;
    uint32_t numGtBins = 4;
    for (int scanPos = firstSigPos; scanPos <= min_sub_pos && maxCtxBins >= 4; scanPos++)
    {
      blk_pos = scan[scanPos];
      int pos_x = blk_pos % width;
      int pos_y = blk_pos / width;
      unsigned absLevel;
      rightPixel = pos_x > 0 ? coeff[pos_x + pos_y * width - 1] : 0;
      belowPixel = pos_y > 0 ? coeff[pos_x + (pos_y - 1) * width] : 0;
      absLevel = uvg_derive_mod_coeff(rightPixel, belowPixel, abs(coeff[blk_pos]), 0);
      cutoffVal = 2;
      for (int j = 0; j < numGtBins; j++)
      {
        if (absLevel >= cutoffVal)
        {
          unsigned gt2 = (absLevel >= (cutoffVal + 2));
          CABAC_FBITS_UPDATE(cabac, &cabac->ctx.transform_skip_gt2[cutoffVal >> 1], gt2, bits, "abs_level_gtx_flag");
          maxCtxBins--;
        }
        cutoffVal += 2;
      }
      lastScanPosPass2 = scanPos;
    }

    //===== coeff bypass ====
    for (int scanPos = firstSigPos; scanPos <= min_sub_pos; scanPos++)
    {
      blk_pos = scan[scanPos];
      int pos_x = blk_pos % width;
      int pos_y = blk_pos / width;
      unsigned absLevel;
      rightPixel = pos_x > 0 ? coeff[pos_x + pos_y * width - 1] : 0;
      belowPixel = pos_y > 0 ? coeff[pos_x + (pos_y - 1) * width] : 0;
      cutoffVal = (scanPos <= lastScanPosPass2 ? 10 : (scanPos <= lastScanPosPass1 ? 2 : 0));
      absLevel = uvg_derive_mod_coeff(rightPixel, belowPixel, abs(coeff[blk_pos]), 0 || !cutoffVal);
      
      if (absLevel >= cutoffVal)
      {
        int       rice = 1;
        unsigned  rem = scanPos <= lastScanPosPass1 ? (absLevel - cutoffVal) >> 1 : absLevel;
        bits += uvg_cabac_write_coeff_remain(cabac, rem, rice, 5);

        if (absLevel && scanPos > lastScanPosPass1)
        {
          int sign = coeff[blk_pos] < 0;
          CABAC_BIN_EP(cabac, sign, "coeff_sign_flag");
          bits += 1;
        }
      }
    }
  }
  if (bits_out && cabac->only_count) *bits_out += bits;
}

/**
 * \brief Encode (X,Y) position of the last significant coefficient
 *
 * \param lastpos_x   X component of last coefficient
 * \param lastpos_y   Y component of last coefficient
 * \param width       Block width
 * \param height      Block height
 * \param type        plane type / luminance or chrominance
 * \param scan        scan type (diag, hor, ver) DEPRECATED?
 *
 * This method encodes the X and Y component within a block of the last
 * significant coefficient.
 */
void uvg_encode_last_significant_xy(cabac_data_t * const cabac,
                                    uint8_t lastpos_x, uint8_t lastpos_y,
                                    uint8_t width, uint8_t height,
                                    uint8_t type, uint8_t scan, double* bits_out)
{
  const int index_x = uvg_math_floor_log2(width);
  const int index_y = uvg_math_floor_log2(height);
  const int prefix_ctx[8] = { 0, 0, 0, 3, 6, 10, 15, 21 };
  //ToDo: own ctx_offset and shift for X and Y 
  uint8_t ctx_offset_x = type ? 0 : prefix_ctx[index_x];
  uint8_t ctx_offset_y = type ? 0 : prefix_ctx[index_y];
  uint8_t shift_x = type ? CLIP(0, 2, width >> 3) : (index_x + 1) >> 2;
  uint8_t shift_y = type ? CLIP(0, 2, height >> 3) : (index_y + 1) >> 2;
  double bits = 0;

  cabac_ctx_t *base_ctx_x = (type ? cabac->ctx.cu_ctx_last_x_chroma : cabac->ctx.cu_ctx_last_x_luma);
  cabac_ctx_t *base_ctx_y = (type ? cabac->ctx.cu_ctx_last_y_chroma : cabac->ctx.cu_ctx_last_y_luma);

  const int group_idx_x = g_group_idx[lastpos_x];
  const int group_idx_y = g_group_idx[lastpos_y];

  // x prefix
  int last_x = 0;
  for (; last_x < group_idx_x; last_x++) {
    CABAC_FBITS_UPDATE(cabac, &base_ctx_x[ctx_offset_x + (last_x >> shift_x)], 1, bits, "last_sig_coeff_x_prefix");
  }
  if (group_idx_x < ( /*width == 32 ? g_group_idx[15] : */g_group_idx[MIN(32, (int32_t)width) - 1])) {
    CABAC_FBITS_UPDATE(cabac, &base_ctx_x[ctx_offset_x + (last_x >> shift_x)], 0, bits, "last_sig_coeff_x_prefix");
  }

  // y prefix
  int last_y = 0;
  for (; last_y < group_idx_y; last_y++) {
    CABAC_FBITS_UPDATE(cabac, &base_ctx_y[ctx_offset_y + (last_y >> shift_y)], 1, bits, "last_sig_coeff_y_prefix");
  }
  if (group_idx_y < (/* height == 32 ? g_group_idx[15] : */g_group_idx[MIN(32, (int32_t)height) - 1])) {
    CABAC_FBITS_UPDATE(cabac, &base_ctx_y[ctx_offset_y + (last_y >> shift_y)], 0, bits, "last_sig_coeff_y_prefix");
  }

  // last_sig_coeff_x_suffix
  if (group_idx_x > 3) {
    const int suffix = lastpos_x - g_min_in_group[group_idx_x];
    const int write_bits = (group_idx_x - 2) / 2;
    CABAC_BINS_EP(cabac, suffix, write_bits, "last_sig_coeff_x_suffix");
    if (cabac->only_count) bits += write_bits;
  }

  // last_sig_coeff_y_suffix
  if (group_idx_y > 3) {
    const int suffix = lastpos_y - g_min_in_group[group_idx_y];
    const int write_bits = (group_idx_y - 2) / 2;
    CABAC_BINS_EP(cabac, suffix, write_bits, "last_sig_coeff_y_suffix");
    if (cabac->only_count) bits += write_bits;
  }
  if (cabac->only_count && bits_out) *bits_out += bits;
}

static void encode_chroma_tu(
  encoder_state_t* const state,
  const cu_loc_t * const cu_loc,
  cu_info_t* cur_pu,
  int8_t* scan_idx,
  lcu_coeff_t* coeff,
  uint8_t joint_chroma)
{
  int width_c = cu_loc->chroma_width;
  int height_c = cu_loc->chroma_height;
  int x_local = (cu_loc->x >> 1) % LCU_WIDTH_C;
  int y_local = (cu_loc->y >> 1) % LCU_WIDTH_C;
  cabac_data_t* const cabac = &state->cabac;
  *scan_idx = SCAN_DIAG;
  if(!joint_chroma){
    // const coeff_t *coeff_u = &coeff->u[xy_to_zorder(LCU_WIDTH_C, x_local, y_local)];
    // const coeff_t *coeff_v = &coeff->v[xy_to_zorder(LCU_WIDTH_C, x_local, y_local)];
    coeff_t coeff_u[TR_MAX_WIDTH * TR_MAX_WIDTH];
    coeff_t coeff_v[TR_MAX_WIDTH * TR_MAX_WIDTH];
    uvg_get_sub_coeff(coeff_u, coeff->u, x_local, y_local, cu_loc->chroma_width, cu_loc->chroma_height, LCU_WIDTH_C);
    uvg_get_sub_coeff(coeff_v, coeff->v, x_local, y_local, cu_loc->chroma_width, cu_loc->chroma_height, LCU_WIDTH_C);

    if (cbf_is_set(cur_pu->cbf, COLOR_U)) {
      if(state->encoder_control->cfg.trskip_enable 
        && width_c <= (1 << state->encoder_control->cfg.trskip_max_size)
        && height_c <= (1 << state->encoder_control->cfg.trskip_max_size)){
        cabac->cur_ctx = &cabac->ctx.transform_skip_model_chroma;
        // HEVC only supports transform_skip for Luma
        // TODO: transform skip for chroma blocks
        CABAC_BIN(cabac, (cur_pu->tr_skip >> COLOR_U) & 1, "transform_skip_flag");
      }
      uvg_encode_coeff_nxn(state, &state->cabac, coeff_u, cu_loc, COLOR_U, *scan_idx, cur_pu, NULL);
    }

    if (cbf_is_set(cur_pu->cbf, COLOR_V)) {
      if (state->encoder_control->cfg.trskip_enable 
        && width_c <= (1 << state->encoder_control->cfg.trskip_max_size)
        && height_c <= (1 << state->encoder_control->cfg.trskip_max_size)) {
        cabac->cur_ctx = &cabac->ctx.transform_skip_model_chroma;
        CABAC_BIN(cabac, (cur_pu->tr_skip >> COLOR_V) & 1, "transform_skip_flag");
      }
      uvg_encode_coeff_nxn(state, &state->cabac, coeff_v, cu_loc, COLOR_V, *scan_idx, cur_pu, NULL);
    }
  }
  else {
    coeff_t coeff_uv[TR_MAX_WIDTH * TR_MAX_WIDTH];
    uvg_get_sub_coeff(coeff_uv, coeff->joint_uv, x_local, y_local, cu_loc->chroma_width, cu_loc->chroma_height, LCU_WIDTH_C);
    if (state->encoder_control->cfg.trskip_enable 
      && width_c <= (1 << state->encoder_control->cfg.trskip_max_size)
      && height_c <= (1 << state->encoder_control->cfg.trskip_max_size)) {
      cabac->cur_ctx = &cabac->ctx.transform_skip_model_chroma;
      CABAC_BIN(cabac, 0, "transform_skip_flag");
    }
    uvg_encode_coeff_nxn(state, &state->cabac, coeff_uv, cu_loc, COLOR_V, *scan_idx, cur_pu, NULL);
    
  }
}

static void encode_transform_unit(
  encoder_state_t * const state,
  const cu_loc_t *cu_loc,
  const cu_info_t* cur_pu,
  lcu_coeff_t* coeff,
  bool only_chroma,
  enum uvg_tree_type tree_type,
  bool last_split,
  const cu_loc_t *original_loc,
  const cu_loc_t* const chroma_loc)               // Original cu dimensions, before CU split
{
  const videoframe_t * const frame = state->tile->frame;
  cabac_data_t* const cabac = &state->cabac;
  const int x = cu_loc->x;
  const int y = cu_loc->y;
  const uint8_t width = cu_loc->width;
  const uint8_t height = cu_loc->height;
  const uint8_t width_c  = cu_loc->chroma_width;
  const uint8_t height_c = cu_loc->chroma_height;

  cu_array_t* used_cu_array = tree_type != UVG_CHROMA_T ? frame->cu_array : frame->chroma_cu_array;
  int isp_x = x;
  int isp_y = y;
  uvg_get_isp_cu_arr_coords(&isp_x, &isp_y, MAX(width, height));
  if(cur_pu == NULL) {
    cur_pu = uvg_cu_array_at_const(used_cu_array, isp_x, isp_y);
  }

  int8_t scan_idx = SCAN_DIAG;

  int cbf_y = cbf_is_set(cur_pu->cbf, COLOR_Y);

  if (cbf_y && !only_chroma) {
    int x_local = x % LCU_WIDTH;
    int y_local = y % LCU_WIDTH;
    // const coeff_t *coeff_y = &coeff->y[xy_to_zorder(LCU_WIDTH, x_local, y_local)];
    coeff_t coeff_y[TR_MAX_WIDTH * TR_MAX_WIDTH];
    uvg_get_sub_coeff(coeff_y, coeff->y, x_local, y_local, width, height, LCU_WIDTH);

    // CoeffNxN
    // Residual Coding

    if(state->encoder_control->cfg.trskip_enable 
      && width <= (1 << state->encoder_control->cfg.trskip_max_size) 
      && height <= (1 << state->encoder_control->cfg.trskip_max_size)
      && !(cur_pu->type == CU_INTRA && cur_pu->intra.isp_mode != ISP_MODE_NO_ISP)) {
      cabac->cur_ctx = &cabac->ctx.transform_skip_model_luma;
      CABAC_BIN(cabac, cur_pu->tr_idx == MTS_SKIP, "transform_skip_flag");
      DBG_YUVIEW_VALUE(state->frame->poc, DBG_YUVIEW_TR_SKIP, x, y, width, height, (cur_pu->tr_idx == MTS_SKIP) ? 1 : 0);
    }
    if(cur_pu->tr_idx == MTS_SKIP) {
      uvg_encode_ts_residual(state, cabac, coeff_y, width, height, 0, scan_idx, NULL);      
    }
    else {
      uvg_encode_coeff_nxn(state,
                           cabac,
                           coeff_y,
                           cu_loc,
                           0,
                           scan_idx,
                           (cu_info_t * )cur_pu,
                           NULL);
    }
    if (tree_type == UVG_LUMA_T) return;
  }

  bool joint_chroma = cur_pu->joint_cb_cr != 0;
  if (cur_pu->log2_height + cur_pu->log2_width < 6 && tree_type != UVG_CHROMA_T && !only_chroma) {
    // For size 4x4 luma transform the corresponding chroma transforms are
    // also of size 4x4 covering 8x8 luma pixels. The residual is coded in
    // the last transform unit.
    if ((x % 8 == 0 || y % 8 == 0) || !only_chroma) {
      // Not the last luma transform block so there is nothing more to do.
      return;
    } else {
      // Time to to code the chroma transform blocks. Move to the top-left
      // corner of the block.
      cur_pu = uvg_cu_array_at_const((const cu_array_t *)used_cu_array, x, y);
    }
  }

  bool chroma_cbf_set = cbf_is_set(cur_pu->cbf, COLOR_U) ||
                        cbf_is_set(cur_pu->cbf, COLOR_V);
  if ((chroma_cbf_set || joint_chroma) && last_split && chroma_loc) {
    //Need to drop const to get lfnst constraints
    // Use original dimensions instead of ISP split dimensions
    encode_chroma_tu(state, chroma_loc, (cu_info_t*)cur_pu, &scan_idx, coeff, joint_chroma);
  }
}

/**
 * \param encoder
 * \param x_pu            Prediction units' x coordinate.
 * \param y_pu            Prediction units' y coordinate.
 * \param depth           Depth from LCU.
 * \param parent_coeff_u  What was signaled at previous level for cbf_cb.
 * \param parent_coeff_v  What was signlaed at previous level for cbf_cr.
 */
static void encode_transform_coeff(
  encoder_state_t * const state,
  const cu_loc_t * cu_loc,
  bool only_chroma,
  lcu_coeff_t* coeff,
  const cu_info_t* cur_tu,
  enum uvg_tree_type tree_type,
  bool last_split,
  bool can_skip_last_cbf,
  int *luma_cbf_ctx,
  // Always true except when writing sub partition coeffs (ISP)
  const cu_loc_t * const original_loc,
  const cu_loc_t* const chroma_loc)       // Original dimensions before ISP split
{
  cabac_data_t * const cabac = &state->cabac;

  bool isp_split = cu_loc->x != original_loc->x || cu_loc->y != original_loc->y;
  int x = cu_loc->x;
  int y = cu_loc->y;
  if (isp_split) {
    uvg_get_isp_cu_arr_coords(&x, &y, MAX(cu_loc->width, cu_loc->height));
  }

  //const encoder_control_t *const ctrl = state->encoder_control;
  const videoframe_t * const frame = state->tile->frame;
  const cu_array_t* used_array = tree_type != UVG_CHROMA_T ? frame->cu_array : frame->chroma_cu_array;
  if(cur_tu == NULL) {
    cur_tu = uvg_cu_array_at_const(used_array, x, y);
  }

  const int tr_limit = TR_MAX_WIDTH;
  const bool ver_split = cu_loc->height > tr_limit;
  const bool hor_split = cu_loc->width > tr_limit;

  const int cb_flag_y = tree_type != UVG_CHROMA_T ? cbf_is_set(cur_tu->cbf, COLOR_Y) : 0;
  const int cb_flag_u = tree_type != UVG_LUMA_T ?(cur_tu->joint_cb_cr ? (cur_tu->joint_cb_cr >> 1) & 1 : cbf_is_set(cur_tu->cbf, COLOR_U)) : 0;
  const int cb_flag_v = tree_type != UVG_LUMA_T ? (cur_tu->joint_cb_cr ? cur_tu->joint_cb_cr & 1 : cbf_is_set(cur_tu->cbf, COLOR_V)) : 0;


  if (hor_split || ver_split) {
    enum split_type split;
    if (cu_loc->width > tr_limit && cu_loc->height > tr_limit) {
      split = QT_SPLIT;
    }
    else if (cu_loc->width > tr_limit) {
      split = BT_VER_SPLIT;
    }
    else {
      split = BT_HOR_SPLIT;
    }

    cu_loc_t split_cu_loc[4];
    const int split_count = uvg_get_split_locs(cu_loc, split, split_cu_loc,NULL);
    for (int i = 0; i < split_count; ++i) {
      encode_transform_coeff(state, &split_cu_loc[i], only_chroma,
        coeff, NULL, tree_type, true, false, luma_cbf_ctx, &split_cu_loc[i], chroma_loc ? &split_cu_loc[i] : NULL);
    }
    return;
  }


  // Chroma cb flags are not signaled when one of the following:
  // No chroma.
  // Not the last CU for area of 64 pixels cowered by more than one luma CU.
  // Not the last ISP Split
  if (state->encoder_control->chroma_format != UVG_CSP_400
    && (chroma_loc || only_chroma)
    && tree_type != UVG_LUMA_T
    && last_split) {
    cabac->cur_ctx = &(cabac->ctx.qt_cbf_model_cb[0]);
    CABAC_BIN(cabac, cb_flag_u, "cbf_cb");
    cabac->cur_ctx = &(cabac->ctx.qt_cbf_model_cr[cb_flag_u ? 1 : 0]);
    CABAC_BIN(cabac, cb_flag_v, "cbf_cr");
  }
  // Luma coded block flag is signaled when one of the following:
  // - prediction mode is intra
  // - transform depth > 0
  // - we have chroma coefficients at this level
  // When it is not present, it is inferred to be 1.
  if ((cur_tu->type == CU_INTRA || !PU_IS_TU(cur_tu) || cb_flag_u || cb_flag_v) && !only_chroma && tree_type != UVG_CHROMA_T) {
    if (can_skip_last_cbf && isp_split && last_split) {
      // Do not write luma cbf if first three isp splits have luma cbf 0
    } else {
      cabac->cur_ctx = &(cabac->ctx.qt_cbf_model_luma[*luma_cbf_ctx]);
      CABAC_BIN(cabac, cb_flag_y, "cbf_luma");
      if (PU_IS_TU(cur_tu)) {
        *luma_cbf_ctx = 2 + cb_flag_y;
      }
    }
  }

  if (cb_flag_y | cb_flag_u | cb_flag_v) {
    if (state->must_code_qp_delta && (only_chroma || cb_flag_y || chroma_loc) ) {
      const int qp_pred      = uvg_get_cu_ref_qp(state, cu_loc->x, cu_loc->y, state->last_qp);
      const int qp_delta     = cur_tu->qp - qp_pred;
      // Possible deltaQP range depends on bit depth as stated in HEVC specification.
      assert(qp_delta >= UVG_QP_DELTA_MIN && qp_delta <= UVG_QP_DELTA_MAX && "QP delta not in valid range.");

      const int qp_delta_abs = ABS(qp_delta);
      cabac_data_t* cabac    = &state->cabac;

      // cu_qp_delta_abs prefix
      uvg_cabac_write_unary_max_symbol(cabac, cabac->ctx.cu_qp_delta_abs, MIN(qp_delta_abs, 5), 1, 5, NULL);

      if (qp_delta_abs >= 5) {
        // cu_qp_delta_abs suffix
        uvg_cabac_write_ep_ex_golomb(state, cabac, qp_delta_abs - 5, 0);
      }

      if (qp_delta != 0) {
        CABAC_BIN_EP(cabac, (qp_delta >= 0 ? 0 : 1), "qp_delta_sign_flag");
      }

      state->must_code_qp_delta = false;
    }
    if((
        ((cb_flag_u || cb_flag_v ) 
          && cur_tu->type == CU_INTRA)
        || (cb_flag_u && cb_flag_v)) 
      && (chroma_loc || only_chroma || tree_type == UVG_CHROMA_T)
      && state->encoder_control->cfg.jccr
      && last_split
      ) {
      assert(cur_tu->joint_cb_cr < 4 && "JointCbCr is in search state.");
      cabac->cur_ctx = &cabac->ctx.joint_cb_cr[cb_flag_u * 2 + cb_flag_v - 1];
      CABAC_BIN(cabac, cur_tu->joint_cb_cr != 0, "tu_joint_cbcr_residual_flag");
    }

    encode_transform_unit(state, cu_loc, only_chroma ? cur_tu : NULL, coeff, only_chroma, tree_type, last_split, original_loc, chroma_loc);
  }
}

/**
 * \brief Writes inter block parameters to the bitstream
 * \param state           Encoder state in use
 * \param x               Slice x coordinate.
 * \param y               Slice y coordinate.
 * \param depth           Depth from LCU.
 * \return if non-zero mvd is coded
 */
int uvg_encode_inter_prediction_unit(
  encoder_state_t * const state,
  cabac_data_t * const cabac,
  const cu_info_t * const cur_cu,
  lcu_t* lcu,
  double* bits_out,
  const cu_loc_t* const cu_loc)
{
  // Mergeflag
  int16_t num_cand = 0;
  bool non_zero_mvd = false;
  double bits = 0;

  CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_merge_flag_ext_model), cur_cu->merged, bits, "MergeFlag");

  num_cand = state->encoder_control->cfg.max_merge;
  if (cur_cu->merged) { //merge
    if (num_cand > 1) {
      int32_t ui;
      for (ui = 0; ui < num_cand - 1; ui++) {
        int32_t symbol = (ui != cur_cu->merge_idx);
        if (ui == 0) {
          CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_merge_idx_ext_model), symbol, bits, "MergeIndex");
        } else {
          CABAC_BIN_EP(cabac,symbol,"MergeIndex");
          if(cabac->only_count) bits += 1;
        }
        if (symbol == 0) break;
      }
    }
#ifdef UVG_DEBUG_PRINT_YUVIEW_CSV
    int abs_x = x + state->tile->offset_x;
    int abs_y = y + state->tile->offset_y;
    if (cur_cu->inter.mv_dir & 1) DBG_YUVIEW_MV(state->frame->poc, DBG_YUVIEW_MVMERGE_L0, abs_x, abs_y, width, height, cur_cu->inter.mv[0][0], cur_cu->inter.mv[0][1]);
    if (cur_cu->inter.mv_dir & 2) DBG_YUVIEW_MV(state->frame->poc, DBG_YUVIEW_MVMERGE_L1, abs_x, abs_y, width, height, cur_cu->inter.mv[1][0], cur_cu->inter.mv[1][1]);
#endif
  } else {
    if (state->frame->slicetype == UVG_SLICE_B && cur_cu->type != CU_IBC) {
      // Code Inter Dir
      uint8_t inter_dir = cur_cu->inter.mv_dir;

      if (cu_loc->width + cu_loc->height > 12) { // ToDo: limit on 4x8/8x4
        uint32_t inter_dir_ctx = (7 - ((uvg_math_floor_log2(cu_loc->width) + uvg_math_floor_log2(cu_loc->height) + 1) >> 1));

        CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.inter_dir[inter_dir_ctx]), (inter_dir == 3), bits, "inter_pred_idc");
      }
      if (inter_dir < 3) {
        CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.inter_dir[5]), (inter_dir == 2), bits, "inter_pred_idc");
      }
   }

    for (uint32_t ref_list_idx = 0; ref_list_idx < 2; ref_list_idx++) {
      if (!(cur_cu->inter.mv_dir & (1 << ref_list_idx))) {
        continue;
      }
#ifdef UVG_DEBUG_PRINT_YUVIEW_CSV
      int abs_x = x + state->tile->offset_x;
      int abs_y = y + state->tile->offset_y;
      DBG_YUVIEW_MV(state->frame->poc, ref_list_idx ? DBG_YUVIEW_MVINTER_L1 : DBG_YUVIEW_MVINTER_L0, abs_x, abs_y, width, height, cur_cu->inter.mv[ref_list_idx][0], cur_cu->inter.mv[ref_list_idx][1]);
#endif
      // size of the current reference index list (L0/L1)
      uint8_t ref_LX_size = state->frame->ref_LX_size[ref_list_idx];

      if (ref_LX_size > 1 && cur_cu->type != CU_IBC) {
        // parseRefFrmIdx
        int32_t ref_frame = cur_cu->inter.mv_ref[ref_list_idx];
        
        CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_ref_pic_model[0]), (ref_frame != 0), bits, "ref_idx_lX");

        if (ref_frame > 0 && ref_LX_size > 2) {
          cabac->cur_ctx = &(cabac->ctx.cu_ref_pic_model[1]);
          CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_ref_pic_model[1]), (ref_frame > 1), bits, "ref_idx_lX");

          if (ref_frame > 1 && ref_LX_size > 3) {
            for (int idx = 3; idx < ref_LX_size; idx++)
            {
              uint8_t val = (ref_frame > idx - 1) ? 1 : 0;
              CABAC_BIN_EP(cabac, val, "ref_idx_lX");
              if (cabac->only_count) bits += 1;
              if (!val) break;

            }
          }
        }

      }

      if (state->frame->ref_list != REF_PIC_LIST_1 || cur_cu->inter.mv_dir != 3) {

        mv_t mv_cand[2][2];
        if (lcu) {
          uvg_inter_get_mv_cand(
            state, 
            mv_cand, cur_cu, lcu, ref_list_idx,
            cu_loc);
        }
        else {
          uvg_inter_get_mv_cand_cua(
            state,
            mv_cand, cur_cu, ref_list_idx, cu_loc
            );
        }

        uint8_t cu_mv_cand = CU_GET_MV_CAND(cur_cu, ref_list_idx);
        mv_t mvd_hor = cur_cu->inter.mv[ref_list_idx][0] - mv_cand[cu_mv_cand][0];
        mv_t mvd_ver = cur_cu->inter.mv[ref_list_idx][1] - mv_cand[cu_mv_cand][1];

        uvg_change_precision(INTERNAL_MV_PREC, uvg_g_imv_to_prec[(cur_cu->type == CU_IBC)?UVG_IMV_FPEL:UVG_IMV_OFF], &mvd_hor, &mvd_ver);
        uvg_encode_mvd(state, cabac, mvd_hor, mvd_ver, bits_out);

        non_zero_mvd |= (mvd_hor != 0) || (mvd_ver != 0);
      }

      // Signal which candidate MV to use
      CABAC_FBITS_UPDATE(cabac,&(cabac->ctx.mvp_idx_model), CU_GET_MV_CAND(cur_cu, ref_list_idx), bits, "mvp_flag");

    } // for ref_list
  } // if !merge
  if(bits_out) *bits_out += bits;
  return non_zero_mvd;
}

static void encode_chroma_intra_cu(
  cabac_data_t* const cabac,
  const cu_info_t* const cur_cu,
  const int cclm_enabled,
  int8_t luma_intra_dir,
  double* bits_out) {
  unsigned pred_mode = 0;
  unsigned chroma_pred_modes[8] = {0, 50, 18, 1, 67, 81, 82, 83};
  int8_t chroma_intra_dir = cur_cu->intra.mode_chroma;
  for(int i = 0; i < 4; i++) {
    if(chroma_pred_modes[i] == luma_intra_dir) {
      chroma_pred_modes[i] = 66;
    }
  }

  double bits = 0;
  bool derived_mode = chroma_intra_dir == luma_intra_dir;
  bool cclm_mode = chroma_intra_dir > 67;

  if (cclm_enabled) {
    CABAC_FBITS_UPDATE(cabac, &cabac->ctx.cclm_flag, cclm_mode, bits, "cclm_flag");
    if(cclm_mode) {
      CABAC_FBITS_UPDATE(cabac, &cabac->ctx.cclm_model, chroma_intra_dir != 81, bits, "cclm_model_1");
      if(chroma_intra_dir != 81) {
        CABAC_BIN_EP(cabac, chroma_intra_dir == 83, "cclm_model_2");
        bits += 1;
      }
      if (cabac->only_count && bits_out) *bits_out += bits;
      return;
    }

  }
  CABAC_FBITS_UPDATE(cabac, &cabac->ctx.chroma_pred_model, derived_mode ? 0 : 1, bits, "intra_chroma_pred_mode");


  if (!derived_mode) {
    /*for (int i = 0; i < 4; i++) {
        if (luma_intra_dir == chroma_pred_modes[i]) {
          chroma_pred_modes[i] = 66;
          break;
        }
      }*/
    for (; pred_mode < 5; pred_mode++) {
      if (chroma_intra_dir == chroma_pred_modes[pred_mode]) {
        break;
      }
    }
    assert(pred_mode != 5 && "Invalid chroma prediction mode");
    /*else if (intra_pred_mode_chroma == 66) {
        // Angular 66 mode is possible only if intra pred mode is one of the
        // possible chroma pred modes, in which case it is signaled with that
        // duplicate mode.
        for (int i = 0; i < 4; ++i) {
          if (intra_pred_mode_actual[0] == chroma_pred_modes[i]) pred_mode = i;
        }
      }
      else {
        for (int i = 0; i < 4; ++i) {
          if (intra_pred_mode_chroma == chroma_pred_modes[i]) pred_mode = i;
        }
      }

      // pred_mode == 67 mean intra_pred_mode_chroma is something that can't
      // be coded.
      assert(pred_mode != 67);
      */
    /**
       * Table 9-35 - Binarization for intra_chroma_pred_mode
       *   intra_chroma_pred_mode  bin_string
       *                        4           0
       *                        0         100
       *                        1         101
       *                        2         110
       *                        3         111
       * Table 9-37 - Assignment of ctxInc to syntax elements with context coded bins
       *   intra_chroma_pred_mode[][] = 0, bypass, bypass
       */
    /*cabac->cur_ctx = &(cabac->ctx.chroma_pred_model[0]);
      if (pred_mode == 68) {
        CABAC_BIN(cabac, 0, "intra_chroma_pred_mode");
      }
      else {
        CABAC_BIN(cabac, 1, "intra_chroma_pred_mode");*/
    CABAC_BINS_EP(cabac, pred_mode, 2, "intra_chroma_pred_mode");
    if (cabac->only_count && bits_out) *bits_out += 2 + bits;
    //}
  }
  else if (cabac->only_count && bits_out)*bits_out += bits;
}

void uvg_encode_intra_luma_coding_unit(
  const encoder_state_t * const state,
  cabac_data_t * const cabac,
  const cu_info_t * const cur_cu,
  const cu_loc_t* const cu_loc,
  const lcu_t* lcu,
  double* bits_out)
{
  const videoframe_t * const frame = state->tile->frame;
  uint8_t intra_pred_mode_actual;
  uint8_t *intra_pred_mode = &intra_pred_mode_actual;

  //uint8_t intra_pred_mode_chroma = cur_cu->intra.mode_chroma;
  int8_t intra_preds[INTRA_MPM_COUNT] = {-1, -1, -1, -1, -1, -1};
  int8_t mpm_preds = -1;
  uint32_t flag;
  double bits = 0;

  const int x = cu_loc->x;
  const int y = cu_loc->y;

  /*
  if ((cur_cu->type == CU_INTRA && (LCU_WIDTH >> cur_cu->depth <= 32))) {
    cabac->cur_ctx = &(cabac->ctx.bdpcm_mode[0]);
    CABAC_BIN(cabac, cur_cu->bdpcmMode > 0 ? 1 : 0, "bdpcm_mode");
    if (cur_cu->bdpcmMode) {
      cabac->cur_ctx = &(cabac->ctx.bdpcm_mode[1]);
      CABAC_BIN(cabac, cur_cu->bdpcmMode > 1 ? 1 : 0, "bdpcm_mode > 1");
    }
  }
  */

  #if ENABLE_PCM == 1
  // Code must start after variable initialization
  uvg_cabac_encode_bin_trm(cabac, 0); // IPCMFlag == 0
  #endif

  /*
  if (cur_cu->type == 1 && (LCU_WIDTH >> depth <= 32)) {
    cabac->cur_ctx = &(cabac->ctx.bdpcm_mode[0]);
    CABAC_BIN(cabac, 0, "bdpcm_mode");
  }
  */
  
  uint32_t width = cu_loc->width;
  uint32_t height = cu_loc->height; // TODO: height for non-square blocks

  // Code MIP related bits
  bool enable_mip = state->encoder_control->cfg.mip;
  int8_t mip_flag = enable_mip ? cur_cu->intra.mip_flag : false;
  bool mip_transpose = enable_mip ? cur_cu->intra.mip_is_transposed : false;
  int8_t mip_mode = enable_mip ? cur_cu->intra.mode : 0;
  uint8_t num_mip_modes;
  
  // Number of MIP modes for this block
  if (width == 4 && height == 4) {
    num_mip_modes = 16;
  }
  else if (width == 4 || height == 4 || (width == 8 && height == 8)) {
    num_mip_modes = 8;
  }
  else {
    num_mip_modes = 6;
  }

  if (mip_flag) {
    assert(mip_mode >= 0 && mip_mode < num_mip_modes && "Invalid MIP mode.");
  }

  if (cur_cu->type == CU_INTRA && !cur_cu->bdpcmMode && enable_mip) {
    uint8_t ctx_id = uvg_get_mip_flag_context(cu_loc, lcu, lcu ? NULL : frame->cu_array);

    // Write MIP flag
    CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.mip_flag[ctx_id]), mip_flag, bits, "mip_flag");
    if (mip_flag) {
      // Write MIP transpose flag & mode
      CABAC_BIN_EP(cabac, mip_transpose, "mip_transposed");
      if (cabac->only_count) bits += 1;
      uvg_cabac_encode_trunc_bin(cabac, mip_mode, num_mip_modes, bits_out);
      if (cabac->only_count && bits_out) *bits_out += bits;
      return;
    }
  }

  // Code MRL related bits
  bool enable_mrl = state->encoder_control->cfg.mrl;
  int multi_ref_idx = enable_mrl ? cur_cu->intra.multi_ref_idx : 0;
  
#ifdef UVG_DEBUG_PRINT_YUVIEW_CSV
  if(multi_ref_idx) DBG_YUVIEW_VALUE(state->frame->poc, DBG_YUVIEW_MRL, x, y, width, height, multi_ref_idx);
#endif

  if (cur_cu->type == CU_INTRA && (y % LCU_WIDTH) != 0 && !cur_cu->bdpcmMode && enable_mrl && !mip_flag) {
    if (MAX_REF_LINE_IDX > 1) {
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.multi_ref_line[0]), multi_ref_idx != 0, bits, "multi_ref_line");
      if (MAX_REF_LINE_IDX > 2 && multi_ref_idx != 0) {
        CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.multi_ref_line[1]), multi_ref_idx != 1, bits, "multi_ref_line");
      }
    }
  }

  bool enable_isp = state->encoder_control->cfg.isp;
  // Need at least 16 samples in sub blocks to use isp. If both dimensions are 4, not enough samples. Blocks of size 2 do not exist yet (not for luma at least)
  bool allow_isp = enable_isp ? uvg_can_use_isp(width, height) : false;
  uint8_t isp_mode = allow_isp ? cur_cu->intra.isp_mode : 0;

  if (allow_isp && !multi_ref_idx /*&& !bdpcm && !color_transform*/) {
    if (isp_mode == ISP_MODE_NO_ISP) {
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.intra_subpart_model[0]), 0, bits, "intra_subpartitions_mode");
    }
    else {
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.intra_subpart_model[0]), 1, bits, "intra_subpartitions_mode");
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.intra_subpart_model[1]), isp_mode - 1, bits, "intra_subpartitions_split_type"); // Vertical or horizontal split
    }
  }
  
    // PREDINFO CODING
    // If intra prediction mode is found from the predictors,
    // it can be signaled with two EP's. Otherwise we can send
    // 5 EP bins with the full predmode
    // ToDo: fix comments for VVC
    
  const cu_info_t* cur_pu = cur_cu; // uvg_cu_array_at_const(frame->cu_array, pu_x, pu_y);

  const cu_info_t* left_pu = NULL;
  const cu_info_t* above_pu = NULL;

  if (x > 0) {
    assert(x >> 2 > 0);
    const int x_scu = SUB_SCU(x) - 1;
    const int y_scu = SUB_SCU(y + height - 1);
    left_pu = lcu ?
                LCU_GET_CU_AT_PX(
                  lcu,
                  x_scu,
                  y_scu) :
                uvg_cu_array_at_const(
                  frame->cu_array,
                  x - 1,
                  y + height - 1);
  }
  // Don't take the above PU across the LCU boundary.
  if (y % LCU_WIDTH > 0 && y > 0) {
    assert(y >> 2 > 0);
    above_pu = lcu ?
                 LCU_GET_CU_AT_PX(
                   lcu,
                   SUB_SCU(x + width - 1),
                   SUB_SCU(y) - 1) :
                 uvg_cu_array_at_const(
                   frame->cu_array,
                   x + width - 1,
                   y - 1);
  }
  
  uvg_intra_get_dir_luma_predictor(x, y,
    intra_preds,
    cur_pu,
    left_pu, above_pu);

  intra_pred_mode_actual = cur_pu->intra.mode;

  for (int i = 0; i < INTRA_MPM_COUNT; i++) {
    if (intra_preds[i] == *intra_pred_mode) {
      mpm_preds = (int8_t)i;
      break;
    }
  }
  // Is the mode in the MPM array or not
  flag = (mpm_preds == -1) ? 0 : 1;
  if (cur_pu->intra.multi_ref_idx == 0) {
    CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.intra_luma_mpm_flag_model), flag, bits, "intra_luma_mpm_flag");
  }
    
  // Signal index of the prediction mode in the prediction list, if it is there
  if (flag) {
    
    const cu_info_t* cur_pu = cur_cu;
    if (cur_pu->intra.multi_ref_idx == 0) {
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.luma_planar_model[(isp_mode ? 0 : 1)]), (mpm_preds > 0 ? 1 : 0), bits, "mpm_idx_luma_planar");
    }

    if (mpm_preds > 0) {
      CABAC_BIN_EP(cabac, (mpm_preds > 1 ? 1 : 0), "mpm_idx");
      if (cabac->only_count) bits += 1;
    }
    if (mpm_preds > 1) {
      CABAC_BIN_EP(cabac, (mpm_preds > 2 ? 1 : 0), "mpm_idx");
      if (cabac->only_count) bits += 1;
    }
    if (mpm_preds > 2) {
      CABAC_BIN_EP(cabac, (mpm_preds > 3 ? 1 : 0), "mpm_idx");
      if (cabac->only_count) bits += 1;
    }
    if (mpm_preds > 3) {
      CABAC_BIN_EP(cabac, (mpm_preds > 4 ? 1 : 0), "mpm_idx");
      if (cabac->only_count) bits += 1;
    }
  }
  else {
    // Signal the actual prediction mode.
    int32_t tmp_pred = *intra_pred_mode;

    uint8_t intra_preds_temp[INTRA_MPM_COUNT + 2];
    memcpy(intra_preds_temp, intra_preds, sizeof(int8_t) * 3);
    memcpy(intra_preds_temp + 4, &intra_preds[3], sizeof(int8_t) * 3);
    intra_preds_temp[3] = 255;
    intra_preds_temp[7] = 255;

    // Improvised merge sort
    // Sort prediction list from lowest to highest.
    if (intra_preds_temp[0] > intra_preds_temp[1]) SWAP(intra_preds_temp[0], intra_preds_temp[1], uint8_t);
    if (intra_preds_temp[0] > intra_preds_temp[2]) SWAP(intra_preds_temp[0], intra_preds_temp[2], uint8_t);
    if (intra_preds_temp[1] > intra_preds_temp[2]) SWAP(intra_preds_temp[1], intra_preds_temp[2], uint8_t);

    if (intra_preds_temp[4] > intra_preds_temp[5]) SWAP(intra_preds_temp[4], intra_preds_temp[5], uint8_t);
    if (intra_preds_temp[4] > intra_preds_temp[6]) SWAP(intra_preds_temp[4], intra_preds_temp[6], uint8_t);
    if (intra_preds_temp[5] > intra_preds_temp[6]) SWAP(intra_preds_temp[5], intra_preds_temp[6], uint8_t);

    // Merge two subarrays
    int32_t array1 = 0;
    int32_t array2 = 4;
    for (int item = 0; item < INTRA_MPM_COUNT; item++) {
      if (intra_preds_temp[array1] < intra_preds_temp[array2]) {
        intra_preds[item] = intra_preds_temp[array1];
        array1++;
      }
      else {
        intra_preds[item] = intra_preds_temp[array2];
        array2++;
      }
    }

    // Reduce the index of the signaled prediction mode according to the
    // prediction list, as it has been already signaled that it's not one
    // of the prediction modes.
    for (int i = INTRA_MPM_COUNT - 1; i >= 0; i--) {
      if (tmp_pred > intra_preds[i]) {
        tmp_pred--;
      }
    }

    uvg_cabac_encode_trunc_bin(cabac, tmp_pred, 67 - INTRA_MPM_COUNT, bits_out);
  }
  if (cabac->only_count && bits_out) *bits_out += bits;
}


uint8_t uvg_write_split_flag(
  const encoder_state_t* const state,
  cabac_data_t* cabac,
  const cu_info_t* left_cu,
  const cu_info_t* above_cu,
  const cu_loc_t* const cu_loc,
  split_tree_t split_tree,
  enum uvg_tree_type tree_type,
  bool* is_implicit_out,
  double* bits_out)
{
  double bits = 0;
  // Implisit split flag when on border
  // Exception made in VVC with flag not being implicit if the BT can be used for
  // horizontal or vertical split, then this flag tells if QT or BT is used
  const int cu_width =  cu_loc->width;
  const int cu_height =  cu_loc->height;


  bool can_split[6];
  const bool is_implicit = uvg_get_possible_splits(state, cu_loc, split_tree, tree_type, can_split);


  bool allow_split = can_split[1] || can_split[2] || can_split[3] || can_split[4] || can_split[5];

  enum split_type split_flag = (split_tree.split_tree >> (split_tree.current_depth * 3)) & 7;

  assert(can_split[split_flag] && "Trying to write an illegal split");

  // split_flag = is_implicit ? (can_split[QT_SPLIT] ? QT_SPLIT : (can_split[BT_HOR_SPLIT] ? BT_HOR_SPLIT : BT_VER_SPLIT)) : split_flag;
  *is_implicit_out = is_implicit;

  int split_model = 0;
  if (can_split[NO_SPLIT] && allow_split) {
    // Get left and top block split_flags and if they are present and true, increase model number
    if (left_cu && (1 << left_cu->log2_height) < cu_height) {
      split_model++;
    }

    if (above_cu && (1 << above_cu->log2_width) < cu_width) {
      split_model++;
    }

    uint32_t split_num = 0;
    if (can_split[QT_SPLIT]) split_num += 2;
    if (can_split[BT_HOR_SPLIT]) split_num++;
    if (can_split[BT_VER_SPLIT]) split_num++;
    if (can_split[TT_HOR_SPLIT]) split_num++;
    if (can_split[TT_VER_SPLIT]) split_num++;

    if (split_num > 0) split_num--;

    split_model += 3 * (split_num >> 1);

    cabac->cur_ctx = &(cabac->ctx.split_flag_model[split_model]);

    CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.split_flag_model[split_model]), split_flag != NO_SPLIT, bits, "split_cu_flag");
  }


  if ((!is_implicit || (can_split[QT_SPLIT] && (can_split[BT_HOR_SPLIT] || can_split[BT_VER_SPLIT]))) 
    && (can_split[BT_HOR_SPLIT] || can_split[BT_VER_SPLIT] || can_split[TT_HOR_SPLIT] || can_split[TT_VER_SPLIT]) 
    && split_flag != NO_SPLIT) {
    bool qt_split = split_flag == QT_SPLIT;
    if((can_split[BT_VER_SPLIT] || can_split[BT_HOR_SPLIT] || can_split[TT_VER_SPLIT] || can_split[TT_HOR_SPLIT]) && can_split[QT_SPLIT]) {
      unsigned left_qt_depth = 0;
      unsigned top_qt_depth = 0;
      if(left_cu) {
        while (((left_cu->split_tree >> (left_qt_depth * 3)) & 7u) == QT_SPLIT) {
          left_qt_depth++;
        }
      }
      if(above_cu) {
        while (((above_cu->split_tree >> (top_qt_depth * 3)) & 7u) == QT_SPLIT) {
          top_qt_depth++;
        }
      }
      split_model = (left_cu && (left_qt_depth > split_tree.current_depth)) + (above_cu && (top_qt_depth > split_tree.current_depth)) + (split_tree.current_depth < 2 ? 0 : 3);
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.qt_split_flag_model[split_model]), qt_split, bits, "qt_split_flag");
    }
    if (!qt_split) {
      const bool is_vertical = split_flag == BT_VER_SPLIT || split_flag == TT_VER_SPLIT;
      if((can_split[BT_HOR_SPLIT] || can_split[TT_HOR_SPLIT]) && (can_split[BT_VER_SPLIT] || can_split[TT_VER_SPLIT])) {
        split_model = 0;
        if(can_split[BT_VER_SPLIT] + can_split[TT_VER_SPLIT] > can_split[BT_HOR_SPLIT] + can_split[TT_HOR_SPLIT]) {
          split_model = 4;
        } else if(can_split[BT_VER_SPLIT] + can_split[TT_VER_SPLIT] < can_split[BT_HOR_SPLIT] + can_split[TT_HOR_SPLIT]) {
          split_model = 3;
        } else {
          const int d_a = cu_width / (above_cu ? (1 << above_cu->log2_width) : 1);
          const int d_l = cu_height / (left_cu ? (1 << left_cu->log2_height) : 1);
          if(d_a != d_l && above_cu && left_cu) {
            split_model = d_a < d_l ? 1 : 2;
          }
        }
        CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.mtt_vertical_model[split_model]), is_vertical, bits, "mtt_vertical_flag");
      }
      if ((can_split[BT_VER_SPLIT] && can_split[TT_VER_SPLIT] && is_vertical) || (can_split[BT_HOR_SPLIT] && can_split[TT_HOR_SPLIT] && !is_vertical)) {
        split_model = (2 * is_vertical) + (split_tree.mtt_depth <= 1);
        CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.mtt_binary_model[split_model]), 
          split_flag == BT_VER_SPLIT || split_flag == BT_HOR_SPLIT, bits, "mtt_binary_flag");
      }
    }
  }

  if (bits_out) *bits_out += bits;
  return split_flag;
}

void uvg_encode_coding_tree(
  encoder_state_t * const state,
  lcu_coeff_t *coeff,
  enum uvg_tree_type tree_type,
  const cu_loc_t* const cu_loc,
  const cu_loc_t* const chroma_loc,
  split_tree_t split_tree,
  bool has_chroma)
{
  cabac_data_t * const cabac = &state->cabac;
  const encoder_control_t * const ctrl = state->encoder_control;
  const videoframe_t * const frame = state->tile->frame;
  const cu_array_t* used_array = tree_type != UVG_CHROMA_T ? frame->cu_array : frame->chroma_cu_array;
  
  const int cu_width  = cu_loc->width;
  const int cu_height = cu_loc->height;
 
  const int x = cu_loc->x;
  const int y = cu_loc->y;

  const cu_info_t* cur_cu = uvg_cu_array_at_const(used_array, x, y);

  const int depth = split_tree.current_depth;

  const cu_info_t *left_cu  = NULL;
  if (x > 0) {
    left_cu = uvg_cu_array_at_const(used_array, x - 1, y);
  }
  const cu_info_t *above_cu = NULL;
  if (y > 0) {
    above_cu = uvg_cu_array_at_const(used_array, x, y - 1);
  }


  // Absolute coordinates
  uint16_t abs_x = x + state->tile->offset_x;
  uint16_t abs_y = y + state->tile->offset_y ;

  int32_t frame_width =  ctrl->in.width;
  int32_t frame_height =  ctrl->in.height;

  // Stop if we are outside of the frame
  if (abs_x >= frame_width || abs_y >= frame_height) return;

  if (depth <= state->frame->max_qp_delta_depth) {
    state->must_code_qp_delta = true;
  }

  // When not in MAX_DEPTH, insert split flag and split the blocks if needed
  if (cu_width + cu_height > 8) {
    split_tree.split_tree = cur_cu->split_tree;
    bool is_implicit;
    const int split_flag = uvg_write_split_flag(
      state,
      cabac,
      left_cu,
      above_cu, 
      tree_type != UVG_CHROMA_T ? cu_loc : chroma_loc,
      split_tree,
      tree_type,
      &is_implicit,
      NULL
      );
    
    if (split_flag != NO_SPLIT) {
      split_tree_t new_split_tree = { cur_cu->split_tree,
        split_tree.current_depth + 1,
        split_tree.mtt_depth + (split_flag != QT_SPLIT),
        split_tree.implicit_mtt_depth + (split_flag != QT_SPLIT && is_implicit),
      0};

      cu_loc_t new_cu_loc[4];
      uint8_t separate_chroma = 0;
      const int splits = uvg_get_split_locs(cu_loc, split_flag, new_cu_loc, &separate_chroma);
      separate_chroma |= !has_chroma;
      for (int split = 0; split <splits; ++split) {
        new_split_tree.part_index = split;
        uvg_encode_coding_tree(state, coeff, tree_type,
          &new_cu_loc[split], 
          separate_chroma ? chroma_loc : &new_cu_loc[split],
          new_split_tree, !separate_chroma || (split == splits - 1 && has_chroma));
      }
      return;
    }
  }
  
  DBG_YUVIEW_VALUE(state->frame->poc, DBG_YUVIEW_CU_TYPE, abs_x, abs_y, cu_width, cu_height, cur_cu->type-1);

  // fprintf(stderr, "%4d %4d %2d %2d %d %d %d\n", x, y, cu_width, cu_height, has_chroma, tree_type, cur_cu->split_tree);

  if (ctrl->cfg.lossless) {
    cabac->cur_ctx = &cabac->ctx.cu_transquant_bypass;
    CABAC_BIN(cabac, 1, "cu_transquant_bypass_flag");
  }

  // Encode skip flag
  if ((state->frame->slicetype != UVG_SLICE_I || state->encoder_control->cfg.ibc)) {

    int8_t ctx_skip = 0;

    if (left_cu && left_cu->skipped) {
      ctx_skip++;
    }
    if (above_cu && above_cu->skipped) {
      ctx_skip++;
    }
    if (cu_width > 4 || state->encoder_control->cfg.ibc) {
      cabac->cur_ctx = &(cabac->ctx.cu_skip_flag_model[ctx_skip]);
      CABAC_BIN(cabac, cur_cu->skipped, "SkipFlag");
    }

    if (cur_cu->skipped) {

      if (state->encoder_control->cfg.ibc && state->frame->slicetype != UVG_SLICE_I)
      { // ToDo: Only for luma channel
        // ToDo: Disable for blocks over 64x64 pixels
        int8_t ctx_ibc = 0;
        if (left_cu && left_cu->type == CU_IBC) ctx_ibc++;
        if (above_cu && above_cu->type == CU_IBC) ctx_ibc++;
        cabac->cur_ctx = &(cabac->ctx.ibc_flag[ctx_ibc]);
        CABAC_BIN(cabac, (cur_cu->type == CU_IBC), "IBCFlag");
      }
      DBG_PRINT_MV(state, x, y, (uint32_t)cu_width, (uint32_t)cu_height, cur_cu);
      uvg_hmvp_add_mv(state, x, y, cu_width, cu_height, cur_cu);
      int16_t num_cand = state->encoder_control->cfg.max_merge;
      if (num_cand > 1) {
        for (int ui = 0; ui < num_cand - 1; ui++) {
          int32_t symbol = (ui != cur_cu->merge_idx);
          if (ui == 0) {
            cabac->cur_ctx = &(cabac->ctx.cu_merge_idx_ext_model);
            CABAC_BIN(cabac, symbol, "MergeIndex");
          } else {
            CABAC_BIN_EP(cabac,symbol,"MergeIndex");
          }
          if (symbol == 0) {
            break;
          }
        }
      }
#ifdef UVG_DEBUG_PRINT_YUVIEW_CSV
      if (cur_cu->inter.mv_dir & 1) DBG_YUVIEW_MV(state->frame->poc, DBG_YUVIEW_MVSKIP_L0, abs_x, abs_y, cu_width, cu_height, cur_cu->inter.mv[0][0], cur_cu->inter.mv[0][1]);
      if (cur_cu->inter.mv_dir & 2) DBG_YUVIEW_MV(state->frame->poc, DBG_YUVIEW_MVSKIP_L1, abs_x, abs_y, cu_width, cu_height, cur_cu->inter.mv[1][0], cur_cu->inter.mv[1][1]);
#endif

      goto end;
    }
  }

  // Prediction mode
  if ((state->frame->slicetype == UVG_SLICE_I || cu_width == 4) && state->encoder_control->cfg.ibc) { // ToDo: Only for luma channel
    // ToDo: Disable for blocks over 64x64 pixels
    int8_t ctx_ibc = 0;
    if (left_cu && left_cu->type == CU_IBC) ctx_ibc++;
    if (above_cu && above_cu->type == CU_IBC) ctx_ibc++;
    cabac->cur_ctx = &(cabac->ctx.ibc_flag[ctx_ibc]);
    CABAC_BIN(cabac, (cur_cu->type == CU_IBC), "IBCFlag");
  }

  if (state->frame->slicetype != UVG_SLICE_I && cu_width != 4 && cu_height != 4)  {

    int8_t ctx_predmode = 0;

    if ((left_cu && left_cu->type == CU_INTRA) || (above_cu && above_cu->type == CU_INTRA)) {
      ctx_predmode=1;
    }

    cabac->cur_ctx = &(cabac->ctx.cu_pred_mode_model[ctx_predmode]);
    CABAC_BIN(cabac, (cur_cu->type == CU_INTRA), "PredMode");

    // We need IBC flag if the mode is signalled as Inter
    if (state->encoder_control->cfg.ibc && cur_cu->type != CU_INTRA) {
      int8_t ctx_ibc = 0;
      if (left_cu && left_cu->type == CU_IBC) ctx_ibc++;
      if (above_cu && above_cu->type == CU_IBC) ctx_ibc++;
      cabac->cur_ctx = &(cabac->ctx.ibc_flag[ctx_ibc]);
      CABAC_BIN(cabac, (cur_cu->type == CU_IBC), "IBCFlag");
    }
  }
    

#if ENABLE_PCM
  // Code IPCM block
  if (FORCE_PCM || cur_cu->type == CU_PCM) {
    uvg_cabac_encode_bin_trm(cabac, 1); // IPCMFlag == 1
    uvg_cabac_finish(cabac);
    uvg_bitstream_add_rbsp_trailing_bits(cabac->stream);
    
    // PCM sample
    uvg_pixel *base_y = &frame->source->y[x + y * ctrl->in.width];
    uvg_pixel *base_u = &frame->source->u[x / 2 + y / 2 * ctrl->in.width / 2];
    uvg_pixel *base_v = &frame->source->v[x / 2 + y / 2 * ctrl->in.width / 2];

    uvg_pixel *rec_base_y = &frame->rec->y[x + y * ctrl->in.width];
    uvg_pixel *rec_base_u = &frame->rec->u[x / 2 + y / 2 * ctrl->in.width / 2];
    uvg_pixel *rec_base_v = &frame->rec->v[x / 2 + y / 2 * ctrl->in.width / 2];

    // Luma
    for (unsigned y_px = 0; y_px < cu_height; y_px++) {
      for (unsigned x_px = 0; x_px < cu_width; x_px++) {
        uvg_bitstream_put(cabac->stream, base_y[x_px + y_px * ctrl->in.width], 8);
        rec_base_y[x_px + y_px * ctrl->in.width] = base_y[x_px + y_px * ctrl->in.width];
      }
    }

    // Chroma
    if (ctrl->chroma_format != UVG_CSP_400) {
      for (unsigned y_px = 0; y_px < cu_loc->chroma_height; y_px++) {
        for (unsigned x_px = 0; x_px < cu_loc->chroma_width; x_px++) {
          uvg_bitstream_put(cabac->stream, base_u[x_px + y_px * (ctrl->in.width >> 1)], 8);
          rec_base_u[x_px + y_px * (ctrl->in.width >> 1)] = base_u[x_px + y_px * (ctrl->in.width >> 1)];
        }
      }
      for (unsigned y_px = 0; y_px < cu_loc->chroma_height; y_px++) {
        for (unsigned x_px = 0; x_px < cu_loc->chroma_width; x_px++) {
          uvg_bitstream_put(cabac->stream, base_v[x_px + y_px * (ctrl->in.width >> 1)], 8);
          rec_base_v[x_px + y_px * (ctrl->in.width >> 1)] = base_v[x_px + y_px * (ctrl->in.width >> 1)];
        }
      }
    }
    uvg_cabac_start(cabac);
  } else 
#endif

  if (cur_cu->type == CU_INTER || cur_cu->type == CU_IBC) {
    uint8_t imv_mode = UVG_IMV_OFF;
    bool non_zero_mvd = false;
  
    // TODO: height for non-square blocks
    const cu_info_t *cur_pu = uvg_cu_array_at_const(used_array, cu_loc->x, cu_loc->y);

    non_zero_mvd |= uvg_encode_inter_prediction_unit(state, cabac, cur_pu, NULL, NULL, cu_loc);
    DBG_PRINT_MV(state, x, y, cu_width, cu_height, cur_pu);
    uvg_hmvp_add_mv(state, x, y, cu_width, cu_height, cur_pu);
    

    // imv mode, select between fullpel, half-pel and quarter-pel resolutions
    // 0 = off, 1 = fullpel, 2 = 4-pel, 3 = half-pel
    if (ctrl->cfg.amvr && non_zero_mvd) {
      cabac->cur_ctx = &(cabac->ctx.imv_flag[0]);
      if(cur_cu->type != CU_IBC) CABAC_BIN(cabac, (imv_mode > UVG_IMV_OFF), "imv_flag");
      if (imv_mode > UVG_IMV_OFF) {
        cabac->cur_ctx = &(cabac->ctx.imv_flag[4]);
        if(cur_cu->type != CU_IBC) CABAC_BIN(cabac, (imv_mode < UVG_IMV_HPEL), "imv_flag");
        if (imv_mode < UVG_IMV_HPEL) {
          cabac->cur_ctx = &(cabac->ctx.imv_flag[1]);
          CABAC_BIN(cabac, (imv_mode > UVG_IMV_FPEL), "imv_flag"); // 1 indicates 4PEL, 0 FPEL
        }
      }
    }

    {
      // Only need to signal coded block flag if not skipped or merged
      // skip = no coded residual, merge = coded residual
      const bool has_coeffs = cur_pu->root_cbf || cur_pu->cbf;
      if (!cur_cu->merged) {
        cabac->cur_ctx = &(cabac->ctx.cu_qt_root_cbf_model);
        CABAC_BIN(cabac, has_coeffs, "rqt_root_cbf");
      }
      // Code (possible) coeffs to bitstream
      if (has_coeffs) {
        int luma_cbf_ctx = 0;
        encode_transform_coeff(state, cu_loc, 0, coeff, cur_cu, tree_type, true, false, &luma_cbf_ctx, cu_loc, cu_loc);
      }

      encode_mts_idx(state, cabac, cur_cu, cu_loc);

    }
  } else if (cur_cu->type == CU_INTRA) {
    if(tree_type != UVG_CHROMA_T) {
      uvg_encode_intra_luma_coding_unit(state, cabac, cur_cu, cu_loc, NULL, NULL);
    }
    
    const bool is_local_dual_tree = (chroma_loc->width != cu_loc->width || chroma_loc->height != cu_loc->height);

    // Code chroma prediction mode.
    if (state->encoder_control->chroma_format != UVG_CSP_400 
      && (chroma_loc->width == cu_loc->width && chroma_loc->height == cu_loc->height) 
      && tree_type == UVG_BOTH_T) {
      encode_chroma_intra_cu(cabac, cur_cu, state->encoder_control->cfg.cclm, !cur_cu->intra.mip_flag ? cur_cu->intra.mode : 0, NULL);
    }
    int luma_cbf_ctx = 0;

    if (tree_type != UVG_CHROMA_T) {
      // Cycle through sub partitions if ISP enabled.
      // ISP split is done horizontally or vertically depending on ISP mode, 2 or 4 times depending on block dimensions.
      // Small blocks are split only twice.
      int split_type = cur_cu->intra.isp_mode;
      int split_limit = split_type == ISP_MODE_NO_ISP ? 1 : uvg_get_isp_split_num(cu_width, cu_height, split_type, true);
      luma_cbf_ctx = split_limit != 1 ? 2 : 0;
      // If all first three splits have luma cbf 0, the last one must be one. Since the value ca be derived, no need to write it
      bool can_skip_last_cbf = true;
      for (int i = 0; i < split_limit; ++i) {
        cu_loc_t split_loc;
        uvg_get_isp_split_loc(&split_loc, x, y, cu_width, cu_height, i, split_type, true);

        // Check if last split to write chroma
        bool last_split = (i + 1) == split_limit;
        encode_transform_coeff(state, &split_loc,
          0, coeff, NULL, tree_type, last_split, can_skip_last_cbf, &luma_cbf_ctx, 
          cu_loc, is_local_dual_tree ? NULL : chroma_loc);
        can_skip_last_cbf &= luma_cbf_ctx == 2;
      }
    }

    if (tree_type != UVG_CHROMA_T) {
      encode_lfnst_idx(state, cabac, cur_cu, is_local_dual_tree && state->encoder_control->chroma_format != UVG_CSP_400 ? UVG_LUMA_T : tree_type, COLOR_Y, cu_loc);

      encode_mts_idx(state, cabac, cur_cu, cu_loc);
    }

    // For 4x4 the chroma PU/TU is coded after the last 
    if (state->encoder_control->chroma_format != UVG_CSP_400 &&
      ((is_local_dual_tree &&
      has_chroma) || tree_type == UVG_CHROMA_T) &&
      tree_type != UVG_LUMA_T)   {
      int8_t luma_dir = uvg_get_co_located_luma_mode(tree_type != UVG_CHROMA_T ? chroma_loc : cu_loc, cu_loc, cur_cu, NULL, frame->cu_array, UVG_CHROMA_T);
      encode_chroma_intra_cu(cabac, cur_cu, state->encoder_control->cfg.cclm && uvg_cclm_is_allowed(state, cu_loc, cur_cu, tree_type), luma_dir,NULL);
      // LFNST constraints must be reset here. Otherwise the left over values will interfere when calculating new constraints
      cu_info_t* tmp = uvg_cu_array_at((cu_array_t *)used_array, chroma_loc->x, chroma_loc->y);
      tmp->violates_lfnst_constrained_luma = false;
      tmp->violates_lfnst_constrained_chroma = false;
      tmp->lfnst_last_scan_pos = false;
      encode_transform_coeff(state, chroma_loc, 1, coeff, NULL, tree_type, true, false, &luma_cbf_ctx, chroma_loc, chroma_loc);
      // Write LFNST only once for single tree structure
      encode_lfnst_idx(state, cabac, tmp, is_local_dual_tree ? UVG_CHROMA_T : tree_type, COLOR_UV, chroma_loc);
    }
  }

  else {
    // CU type not set. Should not happen.
    assert(0);
    exit(1);
  }
  if (state->encoder_control->cabac_debug_file) {
    fprintf(state->encoder_control->cabac_debug_file, "E %4d %4d %9d %d", x, y, split_tree.split_tree, tree_type);
    fwrite(&cabac->ctx, 1, sizeof(cabac->ctx), state->encoder_control->cabac_debug_file);
  }

end:

  if (is_last_cu_in_qg(state, cu_loc)) {
    state->last_qp = cur_cu->qp;
  }

}

double uvg_mock_encode_coding_unit(
  encoder_state_t* const state,
  cabac_data_t* cabac,
  const cu_loc_t* const cu_loc,
  const cu_loc_t* const chroma_loc,
  lcu_t* lcu,
  cu_info_t* cur_cu,
  enum uvg_tree_type tree_type,
  const split_tree_t split_tree) {
  double bits = 0;
  const encoder_control_t* const ctrl = state->encoder_control;

  const int x = cu_loc->x;
  const int y = cu_loc->y;

  const uint8_t depth = 6 - uvg_g_convert_to_log2[cu_loc->width];

  int x_local = cu_loc->local_x;
  int y_local = cu_loc->local_y;
  const bool is_separate_tree = chroma_loc == NULL || cu_loc->height != chroma_loc->height || cu_loc->width != chroma_loc->width;
    
  const cu_info_t* left_cu = NULL, *above_cu = NULL;
  if (x) {
    if(x_local || tree_type != UVG_CHROMA_T) {
      left_cu = LCU_GET_CU_AT_PX(lcu, x_local - 1, y_local);
    }
    else {
      left_cu = uvg_cu_array_at_const(state->tile->frame->chroma_cu_array, x  - 1, y);
    }
  }
  if (y) {
    if(y_local || tree_type != UVG_CHROMA_T) {
      above_cu = LCU_GET_CU_AT_PX(lcu, x_local, y_local-1);
    }
    else {
      above_cu = uvg_cu_array_at_const(state->tile->frame->chroma_cu_array, x, y - 1);
    }
  }
  
  if (depth <= state->frame->max_qp_delta_depth) {
    state->must_code_qp_delta = true;
  }

  // When not in MAX_DEPTH, insert split flag and split the blocks if needed
  if (cur_cu->log2_height + cur_cu->log2_width > 4) {
    // We do not care about whether the split is implicit or not since there is never split here
    bool is_implicit;
    uvg_write_split_flag(
      state,
      cabac,
      left_cu,
      above_cu,
      cu_loc,
      split_tree,
      tree_type, &is_implicit,
      &bits
      );
  }

  // Encode skip flag
  if (state->frame->slicetype != UVG_SLICE_I && (cu_loc->width != 4 || cu_loc->height != 4)) {
    int8_t ctx_skip = 0;

    if (left_cu && left_cu->skipped) {
      ctx_skip++;
    }
    if (above_cu && above_cu->skipped) {
      ctx_skip++;
    }
    
    CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_skip_flag_model[ctx_skip]), cur_cu->skipped, bits, "SkipFlag");

    if (cur_cu->skipped) {
      int16_t num_cand = state->encoder_control->cfg.max_merge;
      if (num_cand > 1) {
        for (int ui = 0; ui < num_cand - 1; ui++) {
          int32_t symbol = (ui != cur_cu->merge_idx);
          if (ui == 0) {
            CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_merge_idx_ext_model), symbol, bits, "MergeIndex");
          }
          else {
            CABAC_BIN_EP(cabac, symbol, "MergeIndex");
            if(cabac->only_count) bits += 1;
          }
          if (symbol == 0) {
            break;
          }
        }
      }
      return bits;
    }
  }
  // Prediction mode
  if (state->frame->slicetype != UVG_SLICE_I && (cu_loc->width != 4 || cu_loc->height != 4)) {

    int8_t ctx_predmode = 0;

    if ((left_cu && left_cu->type == CU_INTRA) || (above_cu && above_cu->type == CU_INTRA)) {
      ctx_predmode = 1;
    }

    CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_pred_mode_model[ctx_predmode]), (cur_cu->type == CU_INTRA), bits, "PredMode");
  }
  
  if (cur_cu->type == CU_INTER || cur_cu->type == CU_IBC) {
    const uint8_t imv_mode = UVG_IMV_OFF;
    const int non_zero_mvd = uvg_encode_inter_prediction_unit(state, cabac, cur_cu, lcu, &bits, cu_loc);
    if (ctrl->cfg.amvr && non_zero_mvd) {
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.imv_flag[0]), imv_mode, bits, "imv_flag");
      if (imv_mode > UVG_IMV_OFF) {
        CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.imv_flag[4]), imv_mode, bits, "imv_flag");
        if (imv_mode < UVG_IMV_HPEL) {
          CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.imv_flag[1]), imv_mode, bits, "imv_flag"); // 1 indicates 4PEL, 0 FPEL
        }
      }
    }
  }
  else if (cur_cu->type == CU_INTRA) {
    if(tree_type != UVG_CHROMA_T) {
      uvg_encode_intra_luma_coding_unit(state, cabac, cur_cu, cu_loc, lcu, &bits);
    }
    if((chroma_loc || tree_type == UVG_CHROMA_T) && state->encoder_control->chroma_format != UVG_CSP_400 && tree_type != UVG_LUMA_T) {
      int8_t luma_dir = uvg_get_co_located_luma_mode(chroma_loc,cu_loc , cur_cu, tree_type != UVG_CHROMA_T ? lcu : NULL,
              tree_type == UVG_CHROMA_T ? state->tile->frame->cu_array : NULL,
              is_separate_tree ? UVG_CHROMA_T : tree_type);
      encode_chroma_intra_cu(cabac, cur_cu, state->encoder_control->cfg.cclm && uvg_cclm_is_allowed(state, chroma_loc, cur_cu, tree_type), luma_dir, &bits);
    }
  }
  else {
    assert(0 && "Unset cu type");
  }
  return bits;
}


void uvg_encode_mvd(encoder_state_t * const state,
                    cabac_data_t *cabac,
                    int32_t mvd_hor,
                    int32_t mvd_ver, double* bits_out)
{
  const int8_t hor_abs_gr0 = mvd_hor != 0;
  const int8_t ver_abs_gr0 = mvd_ver != 0;
  const uint32_t mvd_hor_abs = abs(mvd_hor);
  const uint32_t mvd_ver_abs = abs(mvd_ver);
  double         temp_bits_out = 0.0;

  cabac->cur_ctx = &cabac->ctx.cu_mvd_model[0];
  CABAC_FBITS_UPDATE(cabac, &cabac->ctx.cu_mvd_model[0], (mvd_hor != 0), temp_bits_out, "abs_mvd_greater0_flag_hor");
  CABAC_FBITS_UPDATE(cabac, &cabac->ctx.cu_mvd_model[0], (mvd_ver != 0), temp_bits_out, "abs_mvd_greater0_flag_ver");

  cabac->cur_ctx = &cabac->ctx.cu_mvd_model[1];
  if (hor_abs_gr0) {
    CABAC_FBITS_UPDATE(cabac, &cabac->ctx.cu_mvd_model[1], (mvd_hor_abs>1), temp_bits_out,"abs_mvd_greater1_flag_hor");
  }
  if (ver_abs_gr0) {
    CABAC_FBITS_UPDATE(cabac, &cabac->ctx.cu_mvd_model[1], (mvd_ver_abs>1), temp_bits_out, "abs_mvd_greater1_flag_ver");
  }

  if (hor_abs_gr0) {
    if (mvd_hor_abs > 1) {
      uint32_t bits = uvg_cabac_write_ep_ex_golomb(state, cabac, mvd_hor_abs - 2, 1);
      if(cabac->only_count) temp_bits_out += bits;
    }
    uint32_t mvd_hor_sign = (mvd_hor > 0) ? 0 : 1;
    CABAC_BIN_EP(cabac, mvd_hor_sign, "mvd_sign_flag_hor");
    if (cabac->only_count) temp_bits_out += 1;
  }
  if (ver_abs_gr0) {
    if (mvd_ver_abs > 1) {
      uint32_t bits = uvg_cabac_write_ep_ex_golomb(state, cabac, mvd_ver_abs - 2, 1);
      if (cabac->only_count) temp_bits_out += bits;
    }
    uint32_t mvd_ver_sign = mvd_ver > 0 ? 0 : 1;
    CABAC_BIN_EP(cabac, mvd_ver_sign, "mvd_sign_flag_ver");
    if (cabac->only_count) temp_bits_out += 1;
  }

  if(bits_out) *bits_out = temp_bits_out;
}


/**
 * \brief Get a subset of LCU coeff array.
 *
 * \param dst         Destination array. Should be coeff_t [32*32].
 * \param src         Coeff LCU array.
 * \param lcu_x       Local LCU x coordinate.
 * \param lcu_y       Local LCU y coordinate.
 * \param width       Block width.
 * \param height      Block height.
 * \param lcu_width   LCU_WIDTH for luma, LCU_WIDTH_C for chroma.
 *
 */
void uvg_get_sub_coeff(const coeff_t *dst, const coeff_t * const src, const int lcu_x, const int lcu_y, const int block_w, const int block_h, const int lcu_width)
{
  // Take subset of coeff array
  coeff_t* dst_ptr = (coeff_t*)dst;
  const coeff_t* coeff_ptr = &src[lcu_x + lcu_y * lcu_width];
  for (int j = 0; j < block_h; ++j) {
    //memcpy(dst_coeff + (j * lcu_width), &coeff[j * tr_width], tr_width * sizeof(coeff_t));
    memcpy(&dst_ptr[j * block_w], &coeff_ptr[j * lcu_width], block_w * sizeof(coeff_t));
  }
}
