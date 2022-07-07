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

bool uvg_is_mts_allowed(const encoder_state_t * const state, cu_info_t *const pred_cu)
{
  uint32_t ts_max_size = 1 << state->encoder_control->cfg.trskip_max_size; 
  const uint32_t max_size = 32; // CU::isIntra(cu) ? MTS_INTRA_MAX_CU_SIZE : MTS_INTER_MAX_CU_SIZE;
  const uint32_t cu_width    = LCU_WIDTH >> pred_cu->depth;
  const uint32_t cu_height   = LCU_WIDTH >> pred_cu->depth;
  //bool mts_allowed = cu.chType == CHANNEL_TYPE_LUMA && compID == COMPONENT_Y;

  uint8_t mts_type = state->encoder_control->cfg.mts;
  bool mts_allowed = mts_type == UVG_MTS_BOTH || (pred_cu->type == CU_INTRA ? mts_type == UVG_MTS_INTRA : pred_cu->type == CU_INTER && mts_type == UVG_MTS_INTER);
  mts_allowed &= cu_width <= max_size && cu_height <= max_size;
  //mts_allowed &= !cu.ispMode; // ISP_TODO: Uncomment this when ISP is implemented.
  //mts_allowed &= !cu.sbtInfo;
  mts_allowed &= !(pred_cu->bdpcmMode && cu_width <= ts_max_size && cu_height <= ts_max_size);
  mts_allowed &= pred_cu->tr_idx != MTS_SKIP && !pred_cu->violates_mts_coeff_constraint && pred_cu->mts_last_scan_pos ;
  mts_allowed &= pred_cu->lfnst_idx == 0;
  return mts_allowed;
}

static void encode_mts_idx(encoder_state_t * const state,
  cabac_data_t * const cabac,
  const cu_info_t *const pred_cu)
{
  //TransformUnit &tu = *cu.firstTU;
  int mts_idx = pred_cu->tr_idx;

  if (uvg_is_mts_allowed(state, (cu_info_t* const )pred_cu) && mts_idx != MTS_SKIP
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

// ISP_TODO: move these defines to a proper place when ISP is implemented
// As of now, these are only needed in lfnst checks
#define NOT_INTRA_SUBPARTITIONS 0
#define HOR_INTRA_SUBPARTITIONS 1
#define VER_INTRA_SUBPARTITIONS 2
#define NUM_INTRA_SUBPARTITIONS_MODES 3
#define INTRA_SUBPARTITIONS_RESERVED 4
#define TU_1D_HOR_SPLIT 8
#define TU_1D_VER_SPLIT 9

#define MIN_TB_SIZE_X 4
#define MIN_TB_SIZE_Y 4

static int get_isp_split_dim(const int width, const int height, const int isp_split_type)
{
  bool divide_tu_in_rows = isp_split_type == TU_1D_HOR_SPLIT;
  uint32_t split_dim_size, non_split_dim_size, partition_size, div_shift = 2;

  if (divide_tu_in_rows)
  {
    split_dim_size = height;
    non_split_dim_size = width;
  }
  else
  {
    split_dim_size = width;
    non_split_dim_size = height;
  }
  
  const unsigned min_num_samples_cu = 1 << ((uvg_math_floor_log2(MIN_TB_SIZE_Y) << 1));
  const unsigned factor_to_min_samples = non_split_dim_size < min_num_samples_cu ? min_num_samples_cu >> uvg_math_floor_log2(non_split_dim_size) : 1;
  partition_size = (split_dim_size >> div_shift) < factor_to_min_samples ? factor_to_min_samples : (split_dim_size >> div_shift);

  assert(!(uvg_math_floor_log2(partition_size) + uvg_math_floor_log2(non_split_dim_size) < uvg_math_floor_log2(min_num_samples_cu)) && "Partition has less than minimum amount of samples.");
  return partition_size;
}

static bool can_use_lfnst_with_isp(const int width, const int height, const int isp_split_type, const enum uvg_tree_type tree_type)
{
  if (tree_type == UVG_CHROMA_T) {
    return false;
  }
  if (isp_split_type == NOT_INTRA_SUBPARTITIONS) {
    return false;
  }

  const int tu_width = (isp_split_type == HOR_INTRA_SUBPARTITIONS) ? width : get_isp_split_dim(width, height, TU_1D_VER_SPLIT);
  const int tu_height = (isp_split_type == HOR_INTRA_SUBPARTITIONS) ? get_isp_split_dim(width, height, TU_1D_HOR_SPLIT) : height;

  if (!(tu_width >= MIN_TB_SIZE_Y && tu_height >= MIN_TB_SIZE_Y))
  {
    return false;
  }
  return true;
}

 bool uvg_is_lfnst_allowed(
  const encoder_state_t* const state,
  const cu_info_t* const pred_cu,
  const int width,
  const int height,
  const int x,
  const int y,
  enum uvg_tree_type tree_type,
  const color_t color,
  const lcu_t* lcu) 
{
  if (state->encoder_control->cfg.lfnst && pred_cu->type == CU_INTRA) {
    const int isp_mode = 0; // ISP_TODO: assign proper ISP mode when ISP is implemented
    const int isp_split_type = 0;
    const int depth = pred_cu->depth;
    const int chroma_width = width >> 1;
    const int chroma_height = height >> 1;
    const int cu_width = tree_type != UVG_LUMA_T || depth == 4 ? width : chroma_width;
    const int cu_height = tree_type != UVG_LUMA_T || depth == 4 ? height : chroma_height;
    bool can_use_lfnst_with_mip = (width >= 16 && height >= 16);
    bool is_sep_tree = depth == 4 || tree_type != UVG_BOTH_T;
    bool mip_flag = pred_cu->type == CU_INTRA && color == COLOR_Y ? pred_cu->intra.mip_flag : false;

    if ((isp_mode && !can_use_lfnst_with_isp(width, height, isp_split_type, tree_type)) ||
      (pred_cu->type == CU_INTRA && mip_flag && !can_use_lfnst_with_mip) || 
      (is_sep_tree && MIN(cu_width, cu_height) < 4) || 
      (cu_width > TR_MAX_WIDTH || cu_height > TR_MAX_WIDTH)) {
      return false;
    }
    bool luma_flag = (depth == 4 && color == COLOR_Y) || (tree_type != UVG_CHROMA_T && depth != 4);
    bool chroma_flag = (depth == 4 && color != COLOR_Y) || tree_type != UVG_LUMA_T;
    bool non_zero_coeff_non_ts_corner_8x8 = (luma_flag && pred_cu->violates_lfnst_constrained_luma) || (chroma_flag && pred_cu->violates_lfnst_constrained_chroma);
    bool is_tr_skip = false;

    const videoframe_t* const frame = state->tile->frame;
    //const int num_pred_units = kvz_part_mode_num_parts[pred_cu->part_size];
    const int tr_depth = pred_cu->tr_depth;
    assert(depth <= tr_depth && "Depth greater than transform depth. This should never trigger.");
    const int num_transform_units = 1 << (2 * (tr_depth - depth));
    const int tu_row_length = 1 << (tr_depth - depth);
    const int tu_width = cu_width >> (tr_depth - depth);
    const int tu_height = tu_width; // TODO: height for non-square blocks

    // TODO: chroma transform skip
    if (color == COLOR_Y) {
      for (int i = 0; i < num_transform_units; i++) {
        // TODO: this works only for square blocks
        const int tu_x = x + ((i % tu_row_length) * tu_width);
        const int tu_y = y + ((i / tu_row_length) * tu_height);
        const cu_info_t* cur_tu = lcu ? LCU_GET_CU_AT_PX(lcu, tu_x, tu_y) : uvg_cu_array_at_const(frame->cu_array, tu_x, tu_y);
        assert(cur_tu != NULL && "NULL transform unit.");
        bool cbf_set = cbf_is_set(cur_tu->cbf, tr_depth, COLOR_Y);

        if (cur_tu != NULL && cbf_set && cur_tu->tr_idx == MTS_SKIP) {
          is_tr_skip = true;
        }
      }
    }

    if ((!pred_cu->lfnst_last_scan_pos && !isp_mode) || non_zero_coeff_non_ts_corner_8x8 || is_tr_skip) {
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
  const int x,
  const int y,
  const int depth,
  const int width,
  const int height,
  enum uvg_tree_type tree_type,
  const color_t color)
{
  
  if (uvg_is_lfnst_allowed(state, pred_cu, width, height, x, y, tree_type, color, NULL)) {
    // Getting separate tree bool from block size is a temporary fix until a proper dual tree check is possible (there is no dual tree structure at time of writing this).
    // VTM seems to force explicit dual tree structure for small 4x4 blocks
    bool is_separate_tree = depth == 4 || tree_type != UVG_BOTH_T;

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
    return false;
  }
}

void uvg_encode_ts_residual(encoder_state_t* const state,
  cabac_data_t* const cabac,
  const coeff_t* coeff,
  uint32_t width,
  uint8_t type,
  int8_t scan_mode,
  double* bits_out) {
  //const encoder_control_t * const encoder = state->encoder_control;
  //int c1 = 1;
  uint32_t i;
  int32_t blk_pos;
  // ToDo: large block support in VVC?
  uint32_t sig_coeffgroup_flag[32 * 32] = { 0 };
  

  // CONSTANTS

  const uint32_t log2_block_size = uvg_g_convert_to_bit[width] + 2;
  const uint32_t log2_cg_size = uvg_g_log2_sbb_size[log2_block_size][log2_block_size][0] + uvg_g_log2_sbb_size[log2_block_size][log2_block_size][1];
  const uint32_t* scan =    uvg_g_sig_last_scan[scan_mode][log2_block_size - 1];
  const uint32_t* scan_cg = g_sig_last_scan_cg[log2_block_size - 1][scan_mode];
  double bits = 0;

  // Init base contexts according to block type
  cabac_ctx_t* base_coeff_group_ctx = &(cabac->ctx.transform_skip_sig_coeff_group[0]);

  cabac->cur_ctx = base_coeff_group_ctx;
  
  int maxCtxBins = (width * width * 7) >> 2;
  unsigned scan_cg_last = (unsigned )-1;
  //unsigned scan_pos_last = (unsigned )-1;

  for (i = 0; i < width * width; i++) {
    if (coeff[scan[i]]) {
      //scan_pos_last = i;
      sig_coeffgroup_flag[scan_cg[i >> log2_cg_size]] = 1;
    }
  }
  scan_cg_last = (width * width - 1) >> log2_cg_size;
  const uint32_t cg_width = (MIN((uint8_t)32, width) >> (log2_cg_size / 2));

  bool no_sig_group_before_last = true;

  for (i = 0; i <= scan_cg_last; i++) {
    if (!(width == 4 || (i ==scan_cg_last && no_sig_group_before_last))) {
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
  const int index_y = uvg_math_floor_log2(width);
  const int prefix_ctx[8] = { 0, 0, 0, 3, 6, 10, 15, 21 };
  //ToDo: own ctx_offset and shift for X and Y 
  uint8_t ctx_offset_x = type ? 0 : prefix_ctx[index_x];
  uint8_t ctx_offset_y = type ? 0 : prefix_ctx[index_y];
  uint8_t shift_x = type ? CLIP(0, 2, width>>3) : (index_x+1)>>2;
  uint8_t shift_y = type ? CLIP(0, 2, width >> 3) : (index_y + 1) >> 2;
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
  int x,
  int y,
  int depth,
  const uint8_t width_c,
  cu_info_t* cur_pu,
  int8_t* scan_idx,
  lcu_coeff_t* coeff,
  uint8_t joint_chroma,
  enum
  uvg_tree_type tree_type)
{
  int x_local = ((x >> (tree_type != UVG_CHROMA_T)) & ~3) % LCU_WIDTH_C;
  int y_local = ((y >> (tree_type != UVG_CHROMA_T)) & ~3) % LCU_WIDTH_C;
  cabac_data_t* const cabac = &state->cabac;
  *scan_idx = uvg_get_scan_order(cur_pu->type, cur_pu->intra.mode_chroma, depth);
  if(!joint_chroma){
    const coeff_t *coeff_u = &coeff->u[xy_to_zorder(LCU_WIDTH_C, x_local, y_local)];
    const coeff_t *coeff_v = &coeff->v[xy_to_zorder(LCU_WIDTH_C, x_local, y_local)];

    if (cbf_is_set(cur_pu->cbf, depth, COLOR_U)) {
      if(state->encoder_control->cfg.trskip_enable && width_c <= (1 << state->encoder_control->cfg.trskip_max_size)){
        cabac->cur_ctx = &cabac->ctx.transform_skip_model_chroma;
        // HEVC only supports transform_skip for Luma
        // TODO: transform skip for chroma blocks
        CABAC_BIN(cabac, (cur_pu->tr_skip >> COLOR_U) & 1, "transform_skip_flag");
      }
      uvg_encode_coeff_nxn(state, &state->cabac, coeff_u, width_c, COLOR_U, *scan_idx, cur_pu, NULL);
    }

    if (cbf_is_set(cur_pu->cbf, depth, COLOR_V)) {
      if (state->encoder_control->cfg.trskip_enable && width_c <= (1 << state->encoder_control->cfg.trskip_max_size)) {
        cabac->cur_ctx = &cabac->ctx.transform_skip_model_chroma;
        CABAC_BIN(cabac, (cur_pu->tr_skip >> COLOR_V) & 1, "transform_skip_flag");
      }
      uvg_encode_coeff_nxn(state, &state->cabac, coeff_v, width_c, COLOR_V, *scan_idx, cur_pu, NULL);
    }
  }
  else {
    const coeff_t *coeff_uv = &coeff->joint_uv[xy_to_zorder(LCU_WIDTH_C, x_local, y_local)];
    if (state->encoder_control->cfg.trskip_enable && width_c <= (1 << state->encoder_control->cfg.trskip_max_size)) {
      cabac->cur_ctx = &cabac->ctx.transform_skip_model_chroma;
      CABAC_BIN(cabac, 0, "transform_skip_flag");
    }
    uvg_encode_coeff_nxn(state, &state->cabac, coeff_uv, width_c, COLOR_V, *scan_idx, cur_pu, NULL);
    
  }
}

static void encode_transform_unit(
  encoder_state_t * const state,
  int x,
  int y,
  int depth,
  bool only_chroma,
  lcu_coeff_t* coeff,
  enum uvg_tree_type tree_type)
{
  assert(depth >= 1 && depth <= MAX_PU_DEPTH);

  const videoframe_t * const frame = state->tile->frame;
  cabac_data_t* const cabac = &state->cabac;
  const uint8_t width = LCU_WIDTH >> depth;
  const uint8_t width_c = (depth == MAX_PU_DEPTH ? width : width / 2);

  cu_array_t* used_cu_array = tree_type != UVG_CHROMA_T ? frame->cu_array : frame->chroma_cu_array;
  const cu_info_t *cur_pu = uvg_cu_array_at_const(used_cu_array, x, y);

  int8_t scan_idx = uvg_get_scan_order(cur_pu->type, cur_pu->intra.mode, depth);

  int cbf_y = cbf_is_set(cur_pu->cbf, depth, COLOR_Y);

  if (cbf_y && !only_chroma) {
    int x_local = x % LCU_WIDTH;
    int y_local = y % LCU_WIDTH;
    const coeff_t *coeff_y = &coeff->y[xy_to_zorder(LCU_WIDTH, x_local, y_local)];

    // CoeffNxN
    // Residual Coding

    if(state->encoder_control->cfg.trskip_enable && width <= (1 << state->encoder_control->cfg.trskip_max_size)) {
      cabac->cur_ctx = &cabac->ctx.transform_skip_model_luma;
      CABAC_BIN(cabac, cur_pu->tr_idx == MTS_SKIP, "transform_skip_flag");
      DBG_YUVIEW_VALUE(state->frame->poc, DBG_YUVIEW_TR_SKIP, x, y, width, width, (cur_pu->tr_idx == MTS_SKIP) ? 1 : 0);
    }
    if(cur_pu->tr_idx == MTS_SKIP) {
      uvg_encode_ts_residual(state, cabac, coeff_y, width, 0, scan_idx, NULL);      
    }
    else {
      uvg_encode_coeff_nxn(state,
                           cabac,
                           coeff_y,
                           width,
                           0,
                           scan_idx,
                           (cu_info_t * )cur_pu,
                           NULL);
    }
  }

  bool joint_chroma = cur_pu->joint_cb_cr != 0;
  if (depth == MAX_DEPTH) {
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

  bool chroma_cbf_set = cbf_is_set(cur_pu->cbf, depth, COLOR_U) ||
                        cbf_is_set(cur_pu->cbf, depth, COLOR_V);
  if (chroma_cbf_set || joint_chroma) {
    //Need to drop const to get lfnst constraints
    encode_chroma_tu(state, x, y, depth, width_c, (cu_info_t*)cur_pu, &scan_idx, coeff, joint_chroma, tree_type);
  }
}

/**
 * \param encoder
 * \param x_pu            Prediction units' x coordinate.
 * \param y_pu            Prediction units' y coordinate.
 * \param depth           Depth from LCU.
 * \param tr_depth        Depth from last CU.
 * \param parent_coeff_u  What was signaled at previous level for cbf_cb.
 * \param parent_coeff_v  What was signlaed at previous level for cbf_cr.
 */
static void encode_transform_coeff(
  encoder_state_t * const state,
  int32_t x,
  int32_t y,
  int8_t depth,
  int8_t tr_depth,
  uint8_t parent_coeff_u,
  uint8_t parent_coeff_v,
  bool only_chroma,
  lcu_coeff_t* coeff,
  enum uvg_tree_type tree_type)
{
  cabac_data_t * const cabac = &state->cabac;
  //const encoder_control_t *const ctrl = state->encoder_control;
  const videoframe_t * const frame = state->tile->frame;
  const cu_array_t* used_array = tree_type != UVG_CHROMA_T ? frame->cu_array : frame->chroma_cu_array;

  const cu_info_t *cur_pu = uvg_cu_array_at_const(used_array, x, y);
  // Round coordinates down to a multiple of 8 to get the location of the
  // containing CU.
  const int x_cu = 8 * (x / 8);
  const int y_cu = 8 * (y / 8);
  const cu_info_t *cur_cu = uvg_cu_array_at_const(used_array, x, y);

  // NxN signifies implicit transform split at the first transform level.
  // There is a similar implicit split for inter, but it is only used when
  // transform hierarchy is not in use.
  //int intra_split_flag = (cur_cu->type == CU_INTRA && cur_cu->part_size == SIZE_NxN);

  // The implicit split by intra NxN is not counted towards max_tr_depth.
  /*
  int max_tr_depth;
  if (cur_cu->type == CU_INTRA) {
    max_tr_depth = ctrl->cfg.tr_depth_intra + intra_split_flag;
  } else {
    max_tr_depth = ctrl->tr_depth_inter;
  }
  */

  int8_t split = (LCU_WIDTH >> depth > TR_MAX_WIDTH);

 

  const int cb_flag_y = tree_type != UVG_CHROMA_T ? cbf_is_set(cur_pu->cbf, depth, COLOR_Y) : 0;
  const int cb_flag_u = tree_type != UVG_LUMA_T ?( cur_pu->joint_cb_cr ? (cur_pu->joint_cb_cr >> 1) & 1 : cbf_is_set(cur_cu->cbf, depth, COLOR_U)) : 0;
  const int cb_flag_v = tree_type != UVG_LUMA_T ? (cur_pu->joint_cb_cr ? cur_pu->joint_cb_cr & 1 : cbf_is_set(cur_cu->cbf, depth, COLOR_V)) : 0;

  // The split_transform_flag is not signaled when:
  // - transform size is greater than 32 (depth == 0)
  // - transform size is 4 (depth == MAX_PU_DEPTH)
  // - transform depth is max
  // - cu is intra NxN and it's the first split
  
  //ToDo: check BMS transform split in QTBT
  /*
  if (depth > 0 &&
      depth < MAX_PU_DEPTH &&
      tr_depth < max_tr_depth &&
      !(intra_split_flag && tr_depth == 0))
  {
    cabac->cur_ctx = &(cabac->ctx.trans_subdiv_model[5 - ((uvg_g_convert_to_bit[LCU_WIDTH] + 2) - depth)]);
    CABAC_BIN(cabac, split, "split_transform_flag");
  }
  */

  // Chroma cb flags are not signaled when one of the following:
  // - transform size is 4 (2x2 chroma transform doesn't exist)
  // - they have already been signaled to 0 previously
  // When they are not present they are inferred to be 0, except for size 4
  // when the flags from previous level are used.
  if (state->encoder_control->chroma_format != UVG_CSP_400 && (depth != 4 || only_chroma) && tree_type != UVG_LUMA_T) {
    
    if (!split) {
      if (true) {
        assert(tr_depth < 5);
        cabac->cur_ctx = &(cabac->ctx.qt_cbf_model_cb[0]);
        CABAC_BIN(cabac, cb_flag_u, "cbf_cb");
      }
      if (true) {
        cabac->cur_ctx = &(cabac->ctx.qt_cbf_model_cr[cb_flag_u ? 1 : 0]);
        CABAC_BIN(cabac,  cb_flag_v, "cbf_cr");
      }
    }
  }

  if (split) {
    uint8_t offset = LCU_WIDTH >> (depth + 1);
    int x2 = x + offset;
    int y2 = y + offset;
    encode_transform_coeff(state, x,  y,  depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v, only_chroma, coeff, tree_type);
    encode_transform_coeff(state, x2, y,  depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v, only_chroma, coeff, tree_type);
    encode_transform_coeff(state, x,  y2, depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v, only_chroma, coeff, tree_type);
    encode_transform_coeff(state, x2, y2, depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v, only_chroma, coeff, tree_type);
    return;
  }

  // Luma coded block flag is signaled when one of the following:
  // - prediction mode is intra
  // - transform depth > 0
  // - we have chroma coefficients at this level
  // When it is not present, it is inferred to be 1.
  if ((cur_cu->type == CU_INTRA || tr_depth > 0 || cb_flag_u || cb_flag_v) && !only_chroma && tree_type != UVG_CHROMA_T) {
      cabac->cur_ctx = &(cabac->ctx.qt_cbf_model_luma[0]);
      CABAC_BIN(cabac, cb_flag_y, "cbf_luma");
  }

  if (cb_flag_y | cb_flag_u | cb_flag_v) {
    if (state->must_code_qp_delta && (only_chroma || cb_flag_y || depth != 4) ) {
      const int qp_pred      = uvg_get_cu_ref_qp(state, x_cu, y_cu, state->last_qp);
      const int qp_delta     = cur_cu->qp - qp_pred;
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
          && cur_cu->type == CU_INTRA)
        || (cb_flag_u && cb_flag_v)) 
      && (depth != 4 || only_chroma || tree_type == UVG_CHROMA_T) 
      && state->encoder_control->cfg.jccr
      ) {
      assert(cur_pu->joint_cb_cr < 4 && "JointCbCr is in search state.");
      cabac->cur_ctx = &cabac->ctx.joint_cb_cr[cb_flag_u * 2 + cb_flag_v - 1];
      CABAC_BIN(cabac, cur_pu->joint_cb_cr != 0, "tu_joint_cbcr_residual_flag");
    }
    encode_transform_unit(state, x, y, depth, only_chroma, coeff, tree_type);
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
int uvg_encode_inter_prediction_unit(encoder_state_t * const state,
                                      cabac_data_t * const cabac,
                                      const cu_info_t * const cur_cu,
                                      int x, int y, int width, int height,
                                      int depth, lcu_t* lcu, double* bits_out)
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
    if (state->frame->slicetype == UVG_SLICE_B) {
      // Code Inter Dir
      uint8_t inter_dir = cur_cu->inter.mv_dir;

      if (cur_cu->part_size == SIZE_2Nx2N || (LCU_WIDTH >> depth) != 4) { // ToDo: limit on 4x8/8x4
        uint32_t inter_dir_ctx = (7 - ((uvg_math_floor_log2(width) + uvg_math_floor_log2(height) + 1) >> 1));

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

      if (ref_LX_size > 1) {
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
            x, y, width, height,
            mv_cand, cur_cu, 
            lcu, ref_list_idx);
        }
        else {
          uvg_inter_get_mv_cand_cua(
            state,
            x, y, width, height,
            mv_cand, cur_cu, ref_list_idx
          );
        }

        uint8_t cu_mv_cand = CU_GET_MV_CAND(cur_cu, ref_list_idx);
        mv_t mvd_hor = cur_cu->inter.mv[ref_list_idx][0] - mv_cand[cu_mv_cand][0];
        mv_t mvd_ver = cur_cu->inter.mv[ref_list_idx][1] - mv_cand[cu_mv_cand][1];

        uvg_change_precision(INTERNAL_MV_PREC, uvg_g_imv_to_prec[UVG_IMV_OFF], &mvd_hor, &mvd_ver);
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
  double* bits_out) {
  unsigned pred_mode = 0;
  unsigned chroma_pred_modes[8] = {0, 50, 18, 1, 67, 81, 82, 83};
  int8_t chroma_intra_dir = cur_cu->intra.mode_chroma;
  int8_t luma_intra_dir = !cur_cu->intra.mip_flag ? cur_cu->intra.mode : 0;
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

void uvg_encode_intra_luma_coding_unit(const encoder_state_t * const state,
                                     cabac_data_t * const cabac,
                                     const cu_info_t * const cur_cu,
                                     int x, int y, int depth, const lcu_t* lcu, double* bits_out)
{
  const videoframe_t * const frame = state->tile->frame;
  uint8_t intra_pred_mode_actual;
  uint8_t *intra_pred_mode = &intra_pred_mode_actual;

  //uint8_t intra_pred_mode_chroma = cur_cu->intra.mode_chroma;
  int8_t intra_preds[INTRA_MPM_COUNT] = {-1, -1, -1, -1, -1, -1};
  int8_t mpm_preds = -1;
  uint32_t flag;
  double bits = 0;

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
  
  // Intra Subpartition mode
  uint32_t width = (LCU_WIDTH >> depth);
  uint32_t height = (LCU_WIDTH >> depth);

  bool enough_samples = uvg_g_convert_to_bit[width] + uvg_g_convert_to_bit[height] > (uvg_g_convert_to_bit[4 /* MIN_TB_SIZEY*/] << 1);
  uint8_t isp_mode = 0;
  // ToDo: add height comparison
  //isp_mode += ((width > TR_MAX_WIDTH) || !enough_samples) ? 1 : 0;
  //isp_mode += ((height > TR_MAX_WIDTH) || !enough_samples) ? 2 : 0;
  bool allow_isp = enough_samples;

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
    const int cu_width = LCU_WIDTH >> depth;
    const int cu_height = cu_width; // TODO: height for non-square blocks
    uint8_t ctx_id = uvg_get_mip_flag_context(x, y, cu_width, cu_height, lcu, lcu ? NULL : frame->cu_array);

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
  if(multi_ref_idx) DBG_YUVIEW_VALUE(state->frame->poc, DBG_YUVIEW_MRL, x, y, width, width, multi_ref_idx);
#endif

  if (cur_cu->type == CU_INTRA && (y % LCU_WIDTH) != 0 && !cur_cu->bdpcmMode && enable_mrl && !mip_flag) {
    if (MAX_REF_LINE_IDX > 1) {
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.multi_ref_line[0]), multi_ref_idx != 0, bits, "multi_ref_line");
      if (MAX_REF_LINE_IDX > 2 && multi_ref_idx != 0) {
        CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.multi_ref_line[1]), multi_ref_idx != 1, bits, "multi_ref_line");
      }
    }
  }


  // ToDo: update real usage, these if clauses as such don't make any sense
  if (isp_mode != 0 && multi_ref_idx == 0) {
    if (isp_mode) {
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.intra_subpart_model[0]),  0, bits, "intra_subPartitions");
    } else {
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.intra_subpart_model[0]), 1, bits, "intra_subPartitions");
      // ToDo: complete this if-clause
      if (isp_mode == 3) {
        CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.intra_subpart_model[0]), allow_isp - 1, bits, "intra_subPart_ver_hor");
      }
    }
  }

  const int cu_width = LCU_WIDTH >> depth;
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
    const int y_scu = SUB_SCU(y + cu_width - 1);
    left_pu = lcu ?
                LCU_GET_CU_AT_PX(
                  lcu,
                  x_scu,
                  y_scu) :
                uvg_cu_array_at_const(
                  frame->cu_array,
                  x - 1,
                  y + cu_width - 1);
  }
  // Don't take the above PU across the LCU boundary.
  if (y % LCU_WIDTH > 0 && y > 0) {
    assert(y >> 2 > 0);
    above_pu = lcu ?
                 LCU_GET_CU_AT_PX(
                   lcu,
                   SUB_SCU(x + cu_width - 1),
                   SUB_SCU(y) - 1) :
                 uvg_cu_array_at_const(
                   frame->cu_array,
                   x + cu_width - 1,
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
  if (!(cur_pu->intra.multi_ref_idx || (isp_mode))) {
    CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.intra_luma_mpm_flag_model), flag, bits, "prev_intra_luma_pred_flag");
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

/**
static void encode_part_mode(encoder_state_t * const state,
                             cabac_data_t * const cabac,
                             const cu_info_t * const cur_cu,
                             int depth)
{
  // Binarization from Table 9-34 of the HEVC spec:
  //
  //                |   log2CbSize >     |    log2CbSize ==
  //                |   MinCbLog2SizeY   |    MinCbLog2SizeY
  // -------+-------+----------+---------+-----------+----------
  //  pred  | part  | AMP      | AMP     |           |
  //  mode  | mode  | disabled | enabled | size == 8 | size > 8
  // -------+-------+----------+---------+-----------+----------
  //  intra | 2Nx2N |        -         - |         1          1
  //        |   NxN |        -         - |         0          0
  // -------+-------+--------------------+----------------------
  //  inter | 2Nx2N |        1         1 |         1          1
  //        |  2NxN |       01       011 |        01         01
  //        |  Nx2N |       00       001 |        00        001
  //        |   NxN |        -         - |         -        000
  //        | 2NxnU |        -      0100 |         -          -
  //        | 2NxnD |        -      0101 |         -          -
  //        | nLx2N |        -      0000 |         -          -
  //        | nRx2N |        -      0001 |         -          -
  // -------+-------+--------------------+----------------------
  //
  //
  // Context indices from Table 9-37 of the HEVC spec:
  //
  //                                      binIdx
  //                               |  0  1  2       3
  // ------------------------------+------------------
  //  log2CbSize == MinCbLog2SizeY |  0  1  2  bypass
  //  log2CbSize >  MinCbLog2SizeY |  0  1  3  bypass
  // ------------------------------+------------------
  double bits = 0;
  if (cur_cu->type == CU_INTRA) {
    if (depth == MAX_DEPTH) {
      cabac->cur_ctx = &(cabac->ctx.part_size_model[0]);
      if (cur_cu->part_size == SIZE_2Nx2N) {
        CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.part_size_model[0]), 1, bits, "part_mode 2Nx2N");
      } else {
        CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.part_size_model[0]), 0, bits, "part_mode NxN");
      }
    }
  } else {

    cabac->cur_ctx = &(cabac->ctx.part_size_model[0]);
    if (cur_cu->part_size == SIZE_2Nx2N) {
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.part_size_model[0]), 1, bits, "part_mode 2Nx2N");
      return bits;
    }
    CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.part_size_model[0]), 0, bits, "part_mode split");

    cabac->cur_ctx = &(cabac->ctx.part_size_model[1]);
    if (cur_cu->part_size == SIZE_2NxN ||
        cur_cu->part_size == SIZE_2NxnU ||
        cur_cu->part_size == SIZE_2NxnD) {
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.part_size_model[1]), 1, bits, "part_mode vertical");
    } else {
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.part_size_model[1]), 0, bits, "part_mode horizontal");
    }

    if (state->encoder_control->cfg.amp_enable && depth < MAX_DEPTH) {
      cabac->cur_ctx = &(cabac->ctx.part_size_model[3]);

      if (cur_cu->part_size == SIZE_2NxN ||
          cur_cu->part_size == SIZE_Nx2N) {
        CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.part_size_model[3]), 1, bits, "part_mode SMP");
        return bits;
      }
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.part_size_model[3]), 0, bits, "part_mode AMP");

      if (cur_cu->part_size == SIZE_2NxnU ||
          cur_cu->part_size == SIZE_nLx2N) {
        CABAC_BINS_EP(cabac, 0, 1, "part_mode AMP");
        if(cabac->only_count) bits += 1;
      } else {
        CABAC_BINS_EP(cabac, 1, 1, "part_mode AMP");
        if(cabac->only_count) bits += 1;
      }
    }
  }
  return bits;
}
**/


bool uvg_write_split_flag(
  const encoder_state_t * const state,
  cabac_data_t* cabac,
  const cu_info_t * left_cu,
  const cu_info_t * above_cu,
  uint8_t split_flag,
  int depth,
  int cu_width,
  int x,
  int y,
  enum uvg_tree_type tree_type,
  double* bits_out)
{
  uint16_t abs_x = x + (state->tile->offset_x >> (tree_type == UVG_CHROMA_T));
  uint16_t abs_y = y + (state->tile->offset_y >> (tree_type == UVG_CHROMA_T));
  double bits = 0;
  const encoder_control_t* const ctrl = state->encoder_control;
  // Implisit split flag when on border
  // Exception made in VVC with flag not being implicit if the BT can be used for
  // horizontal or vertical split, then this flag tells if QT or BT is used

  bool no_split, allow_qt, bh_split, bv_split, th_split, tv_split;
  no_split = allow_qt = bh_split = bv_split = th_split = tv_split = true;
  if (depth > MAX_DEPTH) allow_qt = false;
  // ToDo: update this when btt is actually used
  bool allow_btt = false;// when mt_depth < MAX_BT_DEPTH
  

  uint8_t implicit_split_mode = UVG_NO_SPLIT;
  //bool implicit_split = border;
  bool bottom_left_available = ((abs_y + cu_width - 1) < (ctrl->in.height >> (tree_type == UVG_CHROMA_T)));
  bool top_right_available = ((abs_x + cu_width - 1) < (ctrl->in.width >> (tree_type == UVG_CHROMA_T)));

  if (!bottom_left_available && !top_right_available && allow_qt) {
    implicit_split_mode = UVG_QUAD_SPLIT;
  }
  else if (!bottom_left_available && allow_btt) {
    implicit_split_mode = UVG_HORZ_SPLIT;
  }
  else if (!top_right_available && allow_btt) {
    implicit_split_mode = UVG_VERT_SPLIT;
  }
  else if (!bottom_left_available || !top_right_available) {
    implicit_split_mode = UVG_QUAD_SPLIT;
  }
  
  // Check split conditions
  if (implicit_split_mode != UVG_NO_SPLIT) {
    no_split = th_split = tv_split = false;
    bh_split = (implicit_split_mode == UVG_HORZ_SPLIT);
    bv_split = (implicit_split_mode == UVG_VERT_SPLIT);
  }

  if (!allow_btt) {
    bh_split = bv_split = th_split = tv_split = false;
  }

  bool allow_split = allow_qt | bh_split | bv_split | th_split | tv_split;

  split_flag |= implicit_split_mode != UVG_NO_SPLIT;

  int split_model = 0;
  if (no_split && allow_split) {
    // Get left and top block split_flags and if they are present and true, increase model number
    // ToDo: should use height and width to increase model, PU_GET_W() ?
    if (left_cu && PU_GET_H(left_cu->part_size, LCU_WIDTH >> left_cu->depth, 0) < LCU_WIDTH >> depth) {
      split_model++;
    }

    if (above_cu && PU_GET_W(above_cu->part_size, LCU_WIDTH >> above_cu->depth, 0) < LCU_WIDTH >> depth) {
      split_model++;
    }

    uint32_t split_num = 0;
    if (allow_qt) split_num += 2;
    if (bh_split) split_num++;
    if (bv_split) split_num++;
    if (th_split) split_num++;
    if (tv_split) split_num++;

    if (split_num > 0) split_num--;

    split_model += 3 * (split_num >> 1);

    cabac->cur_ctx = &(cabac->ctx.split_flag_model[split_model]);
    if(cabac->only_count && !split_flag) {

      //printf("%hu %hu %d %d %d\n", state->search_cabac.ctx.split_flag_model[split_model].state[0], state->search_cabac.ctx.split_flag_model[split_model].state[1],
      //  split_model, x, y);
    }
    CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.split_flag_model[split_model]), split_flag, bits, "split_flag");
  }

  bool qt_split = split_flag || implicit_split_mode == UVG_QUAD_SPLIT;

  if (!(implicit_split_mode == UVG_NO_SPLIT) && (allow_qt && allow_btt)) {
    split_model = (left_cu && GET_SPLITDATA(left_cu, depth)) + (above_cu && GET_SPLITDATA(above_cu, depth)) + (depth < 2 ? 0 : 3);
    CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.qt_split_flag_model[split_model]), qt_split, bits, "QT_split_flag");
  }

  // Only signal split when it is not implicit, currently only Qt split supported
  if (!(implicit_split_mode == UVG_NO_SPLIT) && !qt_split && (bh_split | bv_split | th_split | tv_split)) {

    split_model = 0;

    // Get left and top block split_flags and if they are present and true, increase model number
    if (left_cu && GET_SPLITDATA(left_cu, depth) == 1) {
      split_model++;
    }

    if (above_cu && GET_SPLITDATA(above_cu, depth) == 1) {
      split_model++;
    }

    split_model += (depth > 2 ? 0 : 3);
    CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.qt_split_flag_model[split_model]), split_flag, bits, "split_cu_mode");
  }
  if (bits_out) *bits_out += bits;
  return split_flag;
}

void uvg_encode_coding_tree(
  encoder_state_t * const state,
  uint16_t x,
  uint16_t y,
  uint8_t depth,
  lcu_coeff_t *coeff,
  enum uvg_tree_type tree_type)
{
  cabac_data_t * const cabac = &state->cabac;
  const encoder_control_t * const ctrl = state->encoder_control;
  const videoframe_t * const frame = state->tile->frame;
  const cu_array_t* used_array = tree_type != UVG_CHROMA_T ? frame->cu_array : frame->chroma_cu_array;
  const cu_info_t *cur_cu   = uvg_cu_array_at_const(used_array, x, y);

  const int cu_width = tree_type != UVG_CHROMA_T ? LCU_WIDTH >> depth : LCU_WIDTH_C >> depth;
  const int cu_height = cu_width; // TODO: height for non-square blocks
  const int half_cu  = cu_width >> 1;

  const cu_info_t *left_cu  = NULL;
  if (x > 0) {
    left_cu = uvg_cu_array_at_const(used_array, x - 1, y);
  }
  const cu_info_t *above_cu = NULL;
  if (y > 0) {
    above_cu = uvg_cu_array_at_const(used_array, x, y - 1);
  }


  // Absolute coordinates
  uint16_t abs_x = x + (state->tile->offset_x >> (tree_type == UVG_CHROMA_T));
  uint16_t abs_y = y + (state->tile->offset_y >> (tree_type == UVG_CHROMA_T));

  int32_t frame_width = tree_type !=  UVG_CHROMA_T ? ctrl->in.width : ctrl->in.width / 2;
  int32_t frame_height = tree_type != UVG_CHROMA_T ? ctrl->in.height : ctrl->in.height / 2;
  // Check for slice border
  bool border_x = frame_width  < abs_x + cu_width;
  bool border_y = frame_height < abs_y + cu_width;
  bool border_split_x = frame_width  >= abs_x + (LCU_WIDTH >> MAX_DEPTH) + half_cu;
  bool border_split_y = frame_height >= abs_y + (LCU_WIDTH >> MAX_DEPTH) + half_cu;
  bool border = border_x || border_y; /*!< are we in any border CU */

  if (depth <= state->frame->max_qp_delta_depth) {
    state->must_code_qp_delta = true;
  }

  // When not in MAX_DEPTH, insert split flag and split the blocks if needed
  if (depth != MAX_DEPTH && !(tree_type == UVG_CHROMA_T && depth == MAX_DEPTH -1)) {

    const int split_flag = uvg_write_split_flag(state, cabac, left_cu, above_cu, GET_SPLITDATA(cur_cu, depth), depth, cu_width, x, y, tree_type,NULL);
    
    if (split_flag || border) {
      // Split blocks and remember to change x and y block positions
      uvg_encode_coding_tree(state, x, y, depth + 1, coeff, tree_type);

      if (!border_x || border_split_x) {
        uvg_encode_coding_tree(state, x + half_cu, y, depth + 1, coeff, tree_type);
      }
      if (!border_y || border_split_y) {
        uvg_encode_coding_tree(state, x, y + half_cu, depth + 1, coeff, tree_type);
      }
      if (!border || (border_split_x && border_split_y)) {
        uvg_encode_coding_tree(state, x + half_cu, y + half_cu, depth + 1, coeff, tree_type);
      }
      return;
    }
  }

  //ToDo: check if we can actually split
  //ToDo: Implement MT split
  if (depth < MAX_PU_DEPTH)
  {
   // cabac->cur_ctx = &(cabac->ctx.trans_subdiv_model[5 - ((uvg_g_convert_to_bit[LCU_WIDTH] + 2) - depth)]);
   // CABAC_BIN(cabac, 0, "split_transform_flag");
  }

  DBG_YUVIEW_VALUE(state->frame->poc, DBG_YUVIEW_CU_TYPE, abs_x, abs_y, cu_width, cu_width, (cur_cu->type == CU_INTRA)?0:1);

  if (ctrl->cfg.lossless) {
    cabac->cur_ctx = &cabac->ctx.cu_transquant_bypass;
    CABAC_BIN(cabac, 1, "cu_transquant_bypass_flag");
  }

  // Encode skip flag
  if (state->frame->slicetype != UVG_SLICE_I && cu_width != 4) {

    int8_t ctx_skip = 0;

    if (left_cu && left_cu->skipped) {
      ctx_skip++;
    }
    if (above_cu && above_cu->skipped) {
      ctx_skip++;
    }

    cabac->cur_ctx = &(cabac->ctx.cu_skip_flag_model[ctx_skip]);
    CABAC_BIN(cabac, cur_cu->skipped, "SkipFlag");

    if (cur_cu->skipped) {
      DBG_PRINT_MV(state, x, y, (uint32_t)cu_width, (uint32_t)cu_width, cur_cu);
      uvg_hmvp_add_mv(state, x, y, (uint32_t)cu_width, (uint32_t)cu_width, cur_cu);
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
      if (cur_cu->inter.mv_dir & 1) DBG_YUVIEW_MV(state->frame->poc, DBG_YUVIEW_MVSKIP_L0, abs_x, abs_y, cu_width, cu_width, cur_cu->inter.mv[0][0], cur_cu->inter.mv[0][1]);
      if (cur_cu->inter.mv_dir & 2) DBG_YUVIEW_MV(state->frame->poc, DBG_YUVIEW_MVSKIP_L1, abs_x, abs_y, cu_width, cu_width, cur_cu->inter.mv[1][0], cur_cu->inter.mv[1][1]);
#endif

      goto end;
    }
  }

  // Prediction mode
  if (state->frame->slicetype != UVG_SLICE_I && cu_width != 4) {

    int8_t ctx_predmode = 0;

    if ((left_cu && left_cu->type == CU_INTRA) || (above_cu && above_cu->type == CU_INTRA)) {
      ctx_predmode=1;
    }

    cabac->cur_ctx = &(cabac->ctx.cu_pred_mode_model[ctx_predmode]);
    CABAC_BIN(cabac, (cur_cu->type == CU_INTRA), "PredMode");
  }

  // part_mode
  //encode_part_mode(state, cabac, cur_cu, depth);

  

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
    for (unsigned y_px = 0; y_px < LCU_WIDTH >> depth; y_px++) {
      for (unsigned x_px = 0; x_px < LCU_WIDTH >> depth; x_px++) {
        uvg_bitstream_put(cabac->stream, base_y[x_px + y_px * ctrl->in.width], 8);
        rec_base_y[x_px + y_px * ctrl->in.width] = base_y[x_px + y_px * ctrl->in.width];
      }
    }

    // Chroma
    if (ctrl->chroma_format != UVG_CSP_400) {
      for (unsigned y_px = 0; y_px < LCU_WIDTH >> (depth + 1); y_px++) {
        for (unsigned x_px = 0; x_px < LCU_WIDTH >> (depth + 1); x_px++) {
          uvg_bitstream_put(cabac->stream, base_u[x_px + y_px * (ctrl->in.width >> 1)], 8);
          rec_base_u[x_px + y_px * (ctrl->in.width >> 1)] = base_u[x_px + y_px * (ctrl->in.width >> 1)];
        }
      }
      for (unsigned y_px = 0; y_px < LCU_WIDTH >> (depth + 1); y_px++) {
        for (unsigned x_px = 0; x_px < LCU_WIDTH >> (depth + 1); x_px++) {
          uvg_bitstream_put(cabac->stream, base_v[x_px + y_px * (ctrl->in.width >> 1)], 8);
          rec_base_v[x_px + y_px * (ctrl->in.width >> 1)] = base_v[x_px + y_px * (ctrl->in.width >> 1)];
        }
      }
    }
    uvg_cabac_start(cabac);
  } else 
#endif

  if (cur_cu->type == CU_INTER) {
    uint8_t imv_mode = UVG_IMV_OFF;
    
    const int num_pu = uvg_part_mode_num_parts[cur_cu->part_size];
    bool non_zero_mvd = false;

    for (int i = 0; i < num_pu; ++i) {
      const int pu_x = PU_GET_X(cur_cu->part_size, cu_width, x, i);
      const int pu_y = PU_GET_Y(cur_cu->part_size, cu_width, y, i);
      const int pu_w = PU_GET_W(cur_cu->part_size, cu_width, i);
      const int pu_h = PU_GET_H(cur_cu->part_size, cu_width, i);
      const cu_info_t *cur_pu = uvg_cu_array_at_const(used_array, pu_x, pu_y);

      non_zero_mvd |= uvg_encode_inter_prediction_unit(state, cabac, cur_pu, pu_x, pu_y, pu_w, pu_h, depth, NULL, NULL);
      DBG_PRINT_MV(state, pu_x, pu_y, pu_w, pu_h, cur_pu);
      uvg_hmvp_add_mv(state, x, y, pu_w, pu_h, cur_pu);
    }

    // imv mode, select between fullpel, half-pel and quarter-pel resolutions
    // 0 = off, 1 = fullpel, 2 = 4-pel, 3 = half-pel
    if (ctrl->cfg.amvr && non_zero_mvd) {
      cabac->cur_ctx = &(cabac->ctx.imv_flag[0]);
      CABAC_BIN(cabac, (imv_mode > UVG_IMV_OFF), "imv_flag");
      if (imv_mode > UVG_IMV_OFF) {
        cabac->cur_ctx = &(cabac->ctx.imv_flag[4]);
        CABAC_BIN(cabac, (imv_mode < UVG_IMV_HPEL), "imv_flag");
        if (imv_mode < UVG_IMV_HPEL) {
          cabac->cur_ctx = &(cabac->ctx.imv_flag[1]);
          CABAC_BIN(cabac, (imv_mode > UVG_IMV_FPEL), "imv_flag"); // 1 indicates 4PEL, 0 FPEL
        }
      }
    }

    {
      int cbf = cbf_is_set_any(cur_cu->cbf, depth);
      // Only need to signal coded block flag if not skipped or merged
      // skip = no coded residual, merge = coded residual
      if (cur_cu->part_size != SIZE_2Nx2N || !cur_cu->merged) {
        cabac->cur_ctx = &(cabac->ctx.cu_qt_root_cbf_model);
        CABAC_BIN(cabac, cbf, "rqt_root_cbf");
      }
      // Code (possible) coeffs to bitstream

      if (cbf) {
        encode_transform_coeff(state, x, y, depth, 0, 0, 0, 0, coeff, tree_type);
      }

      encode_mts_idx(state, cabac, cur_cu);

    }
  } else if (cur_cu->type == CU_INTRA) {
    if(tree_type != UVG_CHROMA_T) {
      uvg_encode_intra_luma_coding_unit(state, cabac, cur_cu, x, y, depth, NULL, NULL);
    }

    // Code chroma prediction mode.
    if (state->encoder_control->chroma_format != UVG_CSP_400 && depth != 4 && tree_type == UVG_BOTH_T) {
      encode_chroma_intra_cu(cabac, cur_cu, state->encoder_control->cfg.cclm, NULL);
    }

    if (tree_type != UVG_CHROMA_T) {
      encode_transform_coeff(state, x, y, depth, 0, 0, 0, 0, coeff, tree_type);
    }

    if (tree_type != UVG_CHROMA_T) {
      bool lfnst_written = encode_lfnst_idx(state, cabac, cur_cu, x, y, depth, cu_width, cu_height, tree_type, COLOR_Y);
    }
    encode_mts_idx(state, cabac, cur_cu);

    // For 4x4 the chroma PU/TU is coded after the last 
    if (state->encoder_control->chroma_format != UVG_CSP_400 && 
      ((depth == 4 && x % 8 && y % 8) || tree_type == UVG_CHROMA_T) &&
      tree_type != UVG_LUMA_T)   {
      encode_chroma_intra_cu(cabac, cur_cu, state->encoder_control->cfg.cclm, NULL);
      // LFNST constraints must be reset here. Otherwise the left over values will interfere when calculating new constraints
      cu_info_t* tmp = uvg_cu_array_at((cu_array_t*)used_array, x, y);
      tmp->violates_lfnst_constrained_luma = false;
      tmp->violates_lfnst_constrained_chroma = false;
      tmp->lfnst_last_scan_pos = false;
      encode_transform_coeff(state, x, y, depth, 0, 0, 0, 1, coeff, tree_type);
      // Write LFNST only once for single tree structure
      encode_lfnst_idx(state, cabac, tmp, x, y, depth, cu_width, cu_height, tree_type, COLOR_UV);
    }
  }

  else {
    // CU type not set. Should not happen.
    assert(0);
    exit(1);
  }
  if (state->encoder_control->cabac_debug_file) {
    fprintf(state->encoder_control->cabac_debug_file, "E %4d %4d %d %d", x << (tree_type == UVG_CHROMA_T), y << (tree_type == UVG_CHROMA_T), depth, tree_type);
    fwrite(&cabac->ctx, 1, sizeof(cabac->ctx), state->encoder_control->cabac_debug_file);
  }

end:

  if (is_last_cu_in_qg(state, x, y, depth)) {
    state->last_qp = cur_cu->qp;
  }

}

double uvg_mock_encode_coding_unit(
  encoder_state_t* const state,
  cabac_data_t* cabac,
  int x,
  int y,
  int depth,
  lcu_t* lcu,
  cu_info_t* cur_cu,
  enum uvg_tree_type tree_type) {
  double bits = 0;
  const encoder_control_t* const ctrl = state->encoder_control;

  int x_local = SUB_SCU(x) >> (tree_type == UVG_CHROMA_T);
  int y_local = SUB_SCU(y) >> (tree_type == UVG_CHROMA_T);

  const int cu_width = LCU_WIDTH >> depth;
  
  const cu_info_t* left_cu = NULL, *above_cu = NULL;
  if (x) {
    if(x_local || tree_type != UVG_CHROMA_T) {
      left_cu = LCU_GET_CU_AT_PX(lcu, x_local - 1, y_local);
    }
    else {
      left_cu = uvg_cu_array_at_const(state->tile->frame->chroma_cu_array, (x >> 1) - 1, y >> 1);
    }
  }
  if (y) {
    if(y_local || tree_type != UVG_CHROMA_T) {
      above_cu = LCU_GET_CU_AT_PX(lcu, x_local, y_local-1);
    }
    else {
      above_cu = uvg_cu_array_at_const(state->tile->frame->chroma_cu_array, x >> 1, (y >> 1) - 1);
    }
  }
  
  if (depth <= state->frame->max_qp_delta_depth) {
    state->must_code_qp_delta = true;
  }

  // When not in MAX_DEPTH, insert split flag and split the blocks if needed
  if (tree_type != UVG_CHROMA_T ? depth != MAX_DEPTH : depth != MAX_DEPTH - 1) {
    uvg_write_split_flag(
      state,
      cabac,
      left_cu,
      above_cu,
      0,
      depth,
      cu_width >> (tree_type == UVG_CHROMA_T),
      x >> (tree_type == UVG_CHROMA_T),
      y >> (tree_type == UVG_CHROMA_T),
      tree_type,
      &bits);
  }

  // Encode skip flag
  if (state->frame->slicetype != UVG_SLICE_I && cu_width != 4) {
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
  if (state->frame->slicetype != UVG_SLICE_I && cu_width != 4) {

    int8_t ctx_predmode = 0;

    if ((left_cu && left_cu->type == CU_INTRA) || (above_cu && above_cu->type == CU_INTRA)) {
      ctx_predmode = 1;
    }

    CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_pred_mode_model[ctx_predmode]), (cur_cu->type == CU_INTRA), bits, "PredMode");
  }
  
  if (cur_cu->type == CU_INTER) {
    const uint8_t imv_mode = UVG_IMV_OFF;
    const int non_zero_mvd = uvg_encode_inter_prediction_unit(state, cabac, cur_cu, x, y, cu_width, cu_width, depth, lcu, &bits);
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
      uvg_encode_intra_luma_coding_unit(state, cabac, cur_cu, x, y, depth, lcu, &bits);
    }
    if((depth != 4 || (x % 8 != 0 && y % 8 != 0)) && state->encoder_control->chroma_format != UVG_CSP_400 && tree_type != UVG_LUMA_T) {
      encode_chroma_intra_cu(cabac, cur_cu, state->encoder_control->cfg.cclm, &bits);
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

  cabac->cur_ctx = &cabac->ctx.cu_mvd_model[0];
  CABAC_FBITS_UPDATE(cabac, &cabac->ctx.cu_mvd_model[0], (mvd_hor != 0), *bits_out, "abs_mvd_greater0_flag_hor");
  CABAC_FBITS_UPDATE(cabac, &cabac->ctx.cu_mvd_model[0], (mvd_ver != 0), *bits_out, "abs_mvd_greater0_flag_ver");

  cabac->cur_ctx = &cabac->ctx.cu_mvd_model[1];
  if (hor_abs_gr0) {
    CABAC_FBITS_UPDATE(cabac, &cabac->ctx.cu_mvd_model[1], (mvd_hor_abs>1), *bits_out,"abs_mvd_greater1_flag_hor");
  }
  if (ver_abs_gr0) {
    CABAC_FBITS_UPDATE(cabac, &cabac->ctx.cu_mvd_model[1], (mvd_ver_abs>1), *bits_out, "abs_mvd_greater1_flag_ver");
  }

  if (hor_abs_gr0) {
    if (mvd_hor_abs > 1) {
      uint32_t bits = uvg_cabac_write_ep_ex_golomb(state, cabac, mvd_hor_abs - 2, 1);
      if(cabac->only_count) *bits_out += bits;
    }
    uint32_t mvd_hor_sign = (mvd_hor > 0) ? 0 : 1;
    CABAC_BIN_EP(cabac, mvd_hor_sign, "mvd_sign_flag_hor");
    if (cabac->only_count) *bits_out += 1;
  }
  if (ver_abs_gr0) {
    if (mvd_ver_abs > 1) {
      uint32_t bits = uvg_cabac_write_ep_ex_golomb(state, cabac, mvd_ver_abs - 2, 1);
      if (cabac->only_count) *bits_out += bits;
    }
    uint32_t mvd_ver_sign = mvd_ver > 0 ? 0 : 1;
    CABAC_BIN_EP(cabac, mvd_ver_sign, "mvd_sign_flag_ver");
    if (cabac->only_count) *bits_out += 1;
  }
}
