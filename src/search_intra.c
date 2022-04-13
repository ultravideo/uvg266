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

#include "search_intra.h"

#include <limits.h>

#include "cabac.h"
#include "encoder.h"
#include "encoderstate.h"
#include "encode_coding_tree.h"
#include "image.h"
#include "intra.h"
#include "kvazaar.h"
#include "rdo.h"
#include "search.h"
#include "strategies/strategies-picture.h"
#include "videoframe.h"


// Normalize SAD for comparison against SATD to estimate transform skip
// for 4x4 blocks.
#ifndef TRSKIP_RATIO
# define TRSKIP_RATIO 1.7
#endif

/**
* \brief Select mode with the smallest cost.
*/
static INLINE uint8_t select_best_mode_index(const int8_t *modes, const double *costs, uint8_t length)
{
  uint8_t best_index = 0;
  double best_cost = costs[0];
  
  for (uint8_t i = 1; i < length; ++i) {
    if (costs[i] < best_cost) {
      best_cost = costs[i];
      best_index = i;
    }
  }

  return best_index;
}


/**
 * \brief Calculate quality of the reconstruction.
 *
 * \param pred  Predicted pixels in continous memory.
 * \param orig_block  Orignal (target) pixels in continous memory.
 * \param satd_func  SATD function for this block size.
 * \param sad_func  SAD function this block size.
 * \param width  Pixel width of the block.
 *
 * \return  Estimated RD cost of the reconstruction and signaling the
 *     coefficients of the residual.
 */
static double get_cost(encoder_state_t * const state, 
                       kvz_pixel *pred, kvz_pixel *orig_block,
                       cost_pixel_nxn_func *satd_func,
                       cost_pixel_nxn_func *sad_func,
                       int width)
{
  double satd_cost = satd_func(pred, orig_block);
  if (TRSKIP_RATIO != 0 && width <= (1 << state->encoder_control->cfg.trskip_max_size) && state->encoder_control->cfg.trskip_enable) {
    // If the mode looks better with SAD than SATD it might be a good
    // candidate for transform skip. How much better SAD has to be is
    // controlled by TRSKIP_RATIO.

    // Add the offset bit costs of signaling 'luma and chroma use trskip',
    // versus signaling 'luma and chroma don't use trskip' to the SAD cost.
    const cabac_ctx_t *ctx = &state->search_cabac.ctx.transform_skip_model_luma;
    double trskip_bits = CTX_ENTROPY_FBITS(ctx, 1) - CTX_ENTROPY_FBITS(ctx, 0);

    
    // ToDo: Check cost
    if (state->encoder_control->chroma_format != KVZ_CSP_400) {
      ctx = &state->search_cabac.ctx.transform_skip_model_chroma;
      trskip_bits += 2.0 * (CTX_ENTROPY_FBITS(ctx, 1) - CTX_ENTROPY_FBITS(ctx, 0));
    }
    

    double sad_cost = TRSKIP_RATIO * sad_func(pred, orig_block) + state->lambda_sqrt * trskip_bits;
    if (sad_cost < satd_cost) {
      return sad_cost;
    }
  }
  return satd_cost;
}


/**
 * \brief Calculate quality of the reconstruction.
 *
 * \param a  bc
 *
 * \return  
 */
static void get_cost_dual(encoder_state_t * const state, 
                       const pred_buffer preds, const kvz_pixel *orig_block,
                       cost_pixel_nxn_multi_func *satd_twin_func,
                       cost_pixel_nxn_multi_func *sad_twin_func,
                       int width, double *costs_out)
{
  #define PARALLEL_BLKS 2
  unsigned satd_costs[PARALLEL_BLKS] = { 0 };
  satd_twin_func(preds, orig_block, PARALLEL_BLKS, satd_costs);
  costs_out[0] = (double)satd_costs[0];
  costs_out[1] = (double)satd_costs[1];

  // TODO: width and height
  if (TRSKIP_RATIO != 0 && width <= (1 << state->encoder_control->cfg.trskip_max_size) && state->encoder_control->cfg.trskip_enable) {
    // If the mode looks better with SAD than SATD it might be a good
    // candidate for transform skip. How much better SAD has to be is
    // controlled by TRSKIP_RATIO.

    // Add the offset bit costs of signaling 'luma and chroma use trskip',
    // versus signaling 'luma and chroma don't use trskip' to the SAD cost.
    const cabac_ctx_t *ctx = &state->cabac.ctx.transform_skip_model_luma;
    double trskip_bits = CTX_ENTROPY_FBITS(ctx, 1) - CTX_ENTROPY_FBITS(ctx, 0);

    
    // ToDo: check cost
    if (state->encoder_control->chroma_format != KVZ_CSP_400) {
      ctx = &state->cabac.ctx.transform_skip_model_chroma;
      trskip_bits += 2.0 * (CTX_ENTROPY_FBITS(ctx, 1) - CTX_ENTROPY_FBITS(ctx, 0));
    }
    

    unsigned unsigned_sad_costs[PARALLEL_BLKS] = { 0 };
    double sad_costs[PARALLEL_BLKS] = { 0 };
    sad_twin_func(preds, orig_block, PARALLEL_BLKS, unsigned_sad_costs);
    for (int i = 0; i < PARALLEL_BLKS; ++i) {
      sad_costs[i] = TRSKIP_RATIO * (double)unsigned_sad_costs[i] + state->lambda_sqrt * trskip_bits;
      if (sad_costs[i] < (double)satd_costs[i]) {
        costs_out[i] = sad_costs[i];
      }
    }
  }

  #undef PARALLEL_BLKS
}


/**
* \brief Derives mts_last_scan_pos and violates_mts_coeff_constraint for pred_cu.
*
* Deriving mts_last_scan_pos and violates_mts_coeff_constraint is
* based on checking coefficients of current luma block.
*
* \param pred_cu  Current prediction coding unit.
* \param lcu      Current lcu.
* \param depth    Current transform depth.
* \param lcu_px   Position of the top left pixel of current CU within current LCU.
*/
static void derive_mts_constraints(cu_info_t *const pred_cu,
                                   lcu_t *const lcu, const int depth,
                                   const vector2d_t lcu_px)
{
  const int width = LCU_WIDTH >> depth;
  int8_t scan_idx = kvz_get_scan_order(pred_cu->type, pred_cu->intra.mode, depth);
  int32_t i;
  // ToDo: large block support in VVC?
  uint32_t sig_coeffgroup_flag[32 * 32] = { 0 };

  const uint32_t log2_block_size = kvz_g_convert_to_bit[width] + 2;
  const uint32_t log2_cg_size = kvz_g_log2_sbb_size[log2_block_size][log2_block_size][0]
    + kvz_g_log2_sbb_size[log2_block_size][log2_block_size][1];
  const uint32_t *scan = kvz_g_sig_last_scan[scan_idx][log2_block_size - 1];
  const uint32_t *scan_cg = g_sig_last_scan_cg[log2_block_size - 1][scan_idx];
  const coeff_t* coeff = &lcu->coeff.y[xy_to_zorder(LCU_WIDTH, lcu_px.x, lcu_px.y)];

  signed scan_cg_last = -1;
  signed scan_pos_last = -1;

  for (int i = 0; i < width * width; i++) {
    if (coeff[scan[i]]) {
      scan_pos_last = i;
      sig_coeffgroup_flag[scan_cg[i >> log2_cg_size]] = 1;
    }
  }
  if (scan_pos_last < 0) return;

  scan_cg_last = scan_pos_last >> log2_cg_size;

  // significant_coeff_flag
  for (i = scan_cg_last; i >= 0; i--) {
    int32_t cg_blk_pos = scan_cg[i];
    int32_t cg_pos_y = cg_blk_pos / (MIN((uint8_t)32, width) >> (log2_cg_size / 2));
    int32_t cg_pos_x = cg_blk_pos - (cg_pos_y * (MIN((uint8_t)32, width) >> (log2_cg_size / 2)));

    // Encode significant coeff group flag when not the last or the first
    if (i == scan_cg_last || i == 0) {
      sig_coeffgroup_flag[cg_blk_pos] = 1;
    }

    if (sig_coeffgroup_flag[cg_blk_pos]) {
      int32_t min_sub_pos = i << log2_cg_size; // LOG2_SCAN_SET_SIZE;
      int32_t first_sig_pos = (i == scan_cg_last) ? scan_pos_last : (min_sub_pos + (1 << log2_cg_size) - 1);

      //if (is_luma /*&& cur_cu->tr_idx != MTS_SKIP*/)
      {
        pred_cu->mts_last_scan_pos |= first_sig_pos > 0;
      }
    }

    if ((cg_pos_y > 3 || cg_pos_x > 3) && sig_coeffgroup_flag[cg_blk_pos] != 0)
    {
      pred_cu->violates_mts_coeff_constraint = true;
    }
  }
}


/**
* \brief Perform search for best intra transform split configuration.
*
* This function does a recursive search for the best intra transform split
* configuration for a given intra prediction mode.
*
* \return RD cost of best transform split configuration. Splits in lcu->cu.
* \param depth          Current transform depth.
* \param max_depth      Depth to which TR split will be tried.
* \param intra_mode     Intra prediction mode.
* \param cost_treshold  RD cost at which search can be stopped.
* \param mts_mode       Selected MTS mode for current intra mode.
*/
static double search_intra_trdepth(
  encoder_state_t * const state,
  int x_px,
  int y_px,
  int depth,
  int max_depth,
  int cost_treshold,
  intra_search_data_t *const search_data,
  lcu_t *const lcu)
{
  assert(depth >= 0 && depth <= MAX_PU_DEPTH);

  const int width = LCU_WIDTH >> depth;
  const int width_c = width > TR_MIN_WIDTH ? width / 2 : width;

  const int offset = width / 2;
  const vector2d_t lcu_px = { SUB_SCU(x_px), SUB_SCU(y_px) };

  const bool reconstruct_chroma =  (depth != 4 || (depth == 4 && (x_px & 4 && y_px & 4))) && state->encoder_control->chroma_format != KVZ_CSP_400;
  cu_info_t* pred_cu = &search_data->pred_cu;
  cu_info_t* const tr_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x, lcu_px.y);

  struct {
    kvz_pixel y[TR_MAX_WIDTH*TR_MAX_WIDTH];
    kvz_pixel u[TR_MAX_WIDTH*TR_MAX_WIDTH];
    kvz_pixel v[TR_MAX_WIDTH*TR_MAX_WIDTH];
  } nosplit_pixels;
  uint16_t nosplit_cbf = 0;

  double split_cost = INT32_MAX;
  double nosplit_cost = INT32_MAX;

  if (depth > 0) {
    const bool mts_enabled = state->encoder_control->cfg.mts == KVZ_MTS_INTRA || state->encoder_control->cfg.mts == KVZ_MTS_BOTH;
    tr_cu->tr_depth = depth;
    pred_cu->tr_depth = depth;

    nosplit_cost = 0.0;

    cbf_clear(&pred_cu->cbf, depth, COLOR_Y);
    if (reconstruct_chroma) {
      cbf_clear(&pred_cu->cbf, depth, COLOR_U);
      cbf_clear(&pred_cu->cbf, depth, COLOR_V);
    }

    const int8_t chroma_mode = reconstruct_chroma ? pred_cu->intra.mode : -1;
    double best_rd_cost = MAX_INT;
    int best_tr_idx = 0;

    int trafo;
    int num_transforms = 1;
    if (pred_cu->tr_idx != MTS_TR_NUM)
    {
      trafo = pred_cu->tr_idx;
      num_transforms = pred_cu->tr_idx + 1;
    }
    else
    {
      trafo = 0;
      num_transforms = (mts_enabled ? MTS_TR_NUM : 1);
    }
    //TODO: height
    if(state->encoder_control->cfg.trskip_enable && width <= (1 << state->encoder_control->cfg.trskip_max_size) /*&& height == 4*/) {
      num_transforms = MAX(num_transforms, 2);
    }

    for (; trafo < num_transforms; trafo++) {
      pred_cu->tr_idx = trafo;
      if (mts_enabled)
      {
        pred_cu->mts_last_scan_pos = 0;
        pred_cu->violates_mts_coeff_constraint = 0;

        if (trafo == MTS_SKIP && width > (1 << state->encoder_control->cfg.trskip_max_size)) {
          //TODO: parametrize that this is not hardcoded
          // TODO: this probably should currently trip for chroma?
          continue;
        }
      }
     
      kvz_intra_recon_cu(state,
                         x_px, y_px,
                         depth, search_data,
                         pred_cu,
                         lcu);

      // TODO: Not sure if this should be 0 or 1 but at least seems to work with 1
      if (pred_cu->tr_idx > 1)
      {
        derive_mts_constraints(pred_cu, lcu, depth, lcu_px);
        if (pred_cu->violates_mts_coeff_constraint || !pred_cu->mts_last_scan_pos)
        {
          assert(pred_cu->tr_idx == MTS_TR_NUM); //mts mode should not be decided and then not allowed to be used. (might be some exception here)
          continue;
        }
      }

      double rd_cost = kvz_cu_rd_cost_luma(state, lcu_px.x, lcu_px.y, depth, pred_cu, lcu);
      //if (reconstruct_chroma) {
      //  rd_cost += kvz_cu_rd_cost_chroma(state, lcu_px.x, lcu_px.y, depth, pred_cu, lcu);
      //}

      if (rd_cost < best_rd_cost) {
        best_rd_cost = rd_cost;
        best_tr_idx = pred_cu->tr_idx;
      }
    }
    if(reconstruct_chroma) {
      int8_t luma_mode = pred_cu->intra.mode;
      pred_cu->intra.mode = -1;
      pred_cu->intra.mode_chroma = chroma_mode;
      pred_cu->joint_cb_cr= 4; // TODO: Maybe check the jccr mode here also but holy shit is the interface of search_intra_rdo bad currently
      kvz_intra_recon_cu(state,
                         x_px & ~7, y_px & ~7,
                         depth, search_data,
                         pred_cu, 
                         lcu);
      best_rd_cost += kvz_cu_rd_cost_chroma(state, lcu_px.x, lcu_px.y, depth, pred_cu, lcu);
      pred_cu->intra.mode = luma_mode;
    }
    pred_cu->tr_skip = best_tr_idx == MTS_SKIP;
    pred_cu->tr_idx = best_tr_idx;
    nosplit_cost += best_rd_cost;

    // Early stop codition for the recursive search.
    // If the cost of any 1/4th of the transform is already larger than the
    // whole transform, assume that splitting further is a bad idea.
    if (nosplit_cost >= cost_treshold) {
      return nosplit_cost;
    }

    nosplit_cbf = pred_cu->cbf;

    kvz_pixels_blit(lcu->rec.y, nosplit_pixels.y, width, width, LCU_WIDTH, width);
    if (reconstruct_chroma) {
      kvz_pixels_blit(lcu->rec.u, nosplit_pixels.u, width_c, width_c, LCU_WIDTH_C, width_c);
      kvz_pixels_blit(lcu->rec.v, nosplit_pixels.v, width_c, width_c, LCU_WIDTH_C, width_c);
    }
  }

  // Recurse further if all of the following:
  // - Current depth is less than maximum depth of the search (max_depth).
  //   - Maximum transform hierarchy depth is constrained by clipping
  //     max_depth.
  // - Min transform size hasn't been reached (MAX_PU_DEPTH).
  if (depth < max_depth && depth < MAX_PU_DEPTH) {
    split_cost = 0;

    split_cost += search_intra_trdepth(state, x_px, y_px, depth + 1, max_depth, nosplit_cost, search_data, lcu);
    if (split_cost < nosplit_cost) {
      split_cost += search_intra_trdepth(state, x_px + offset, y_px, depth + 1, max_depth, nosplit_cost, search_data, lcu);
    }
    if (split_cost < nosplit_cost) {
      split_cost += search_intra_trdepth(state, x_px, y_px + offset, depth + 1, max_depth, nosplit_cost, search_data, lcu);
    }
    if (split_cost < nosplit_cost) {
      split_cost += search_intra_trdepth(state, x_px + offset, y_px + offset, depth + 1, max_depth, nosplit_cost, search_data, lcu);
    }

    double cbf_bits = 0.0;

    // Add cost of cbf chroma bits on transform tree.
    // All cbf bits are accumulated to pred_cu.cbf and cbf_is_set returns true
    // if cbf is set at any level >= depth, so cbf chroma is assumed to be 0
    // if this and any previous transform block has no chroma coefficients.
    // When searching the first block we don't actually know the real values,
    // so this will code cbf as 0 and not code the cbf at all for descendants.
    if (state->encoder_control->chroma_format != KVZ_CSP_400) {
      const uint8_t tr_depth = depth - pred_cu->depth;
      cabac_data_t* cabac = (cabac_data_t *)&state->search_cabac;

      cabac_ctx_t* ctx = &(cabac->ctx.qt_cbf_model_cb[0]);
      if (tr_depth == 0 || cbf_is_set(pred_cu->cbf, depth - 1, COLOR_U)) {
        CABAC_FBITS_UPDATE(cabac, ctx, cbf_is_set(pred_cu->cbf, depth, COLOR_U), cbf_bits, "cbf_cb");
      }
      ctx = &(state->cabac.ctx.qt_cbf_model_cr[cbf_is_set(pred_cu->cbf, depth, COLOR_U)]);
      if (tr_depth == 0 || cbf_is_set(pred_cu->cbf, depth - 1, COLOR_V)) {
        CABAC_FBITS_UPDATE(cabac, ctx, cbf_is_set(pred_cu->cbf, depth, COLOR_V), cbf_bits, "cbf_cr");
      }
    }

    double bits = cbf_bits;
    split_cost += bits * state->lambda;
  } else {
    assert(width <= TR_MAX_WIDTH);
  }

  if (depth == 0 || split_cost < nosplit_cost) {
    return split_cost;
  } else {
    kvz_lcu_fill_trdepth(lcu, x_px, y_px, depth, depth);

    pred_cu->cbf = nosplit_cbf;

    // We only restore the pixel data and not coefficients or cbf data.
    // The only thing we really need are the border pixels.kvz_intra_get_dir_luma_predictor
    kvz_pixels_blit(nosplit_pixels.y, lcu->rec.y, width, width, width, LCU_WIDTH);
    if (reconstruct_chroma) {
      kvz_pixels_blit(nosplit_pixels.u, lcu->rec.u, width_c, width_c, width_c, LCU_WIDTH_C);
      kvz_pixels_blit(nosplit_pixels.v, lcu->rec.v, width_c, width_c, width_c, LCU_WIDTH_C);
    }

    return nosplit_cost;
  }
}
static void sort_modes(intra_search_data_t* __restrict modes, uint8_t length)
{
  // Length for intra is always between 5 and 23, and is either 21, 17, 9 or 8 about
  // 60% of the time, so there should be no need for anything more complex
  // than insertion sort.
  // Length for merge is 5 or less.
  for (uint8_t i = 1; i < length; ++i) {
    const intra_search_data_t cur_cost = modes[i];
    uint8_t j = i;
    while (j > 0 && cur_cost.cost < modes[j - 1].cost) {
      modes[j] = modes[j - 1];
      --j;
    }
    modes[j] = cur_cost;
  }
}

static void search_intra_chroma_rough(
  encoder_state_t * const state,
  int x_px,
  int y_px,
  int depth,
  const kvz_pixel *orig_u,
  const kvz_pixel *orig_v,
  int16_t origstride,
  kvz_intra_references *refs_u,
  kvz_intra_references *refs_v,
  intra_search_data_t* chroma_data,
  lcu_t* lcu)
{
  assert(!(x_px & 4 || y_px & 4));

  const unsigned width = MAX(LCU_WIDTH_C >> depth, TR_MIN_WIDTH);

  cost_pixel_nxn_func *const satd_func = kvz_pixels_get_satd_func(width);
  //cost_pixel_nxn_func *const sad_func = kvz_pixels_get_sad_func(width);
  cu_loc_t loc = { x_px, y_px, width, width, width, width };
    
  kvz_pixel _pred[32 * 32 + SIMD_ALIGNMENT];
  kvz_pixel *pred = ALIGNED_POINTER(_pred, SIMD_ALIGNMENT);

  kvz_pixel _orig_block[32 * 32 + SIMD_ALIGNMENT];
  kvz_pixel *orig_block = ALIGNED_POINTER(_orig_block, SIMD_ALIGNMENT);

  kvz_pixels_blit(orig_u, orig_block, width, width, origstride, width);
  int modes_count = (state->encoder_control->cfg.cclm ? 8 : 5);
  for (int i = 0; i < modes_count; ++i) {
    if (chroma_data[i].pred_cu.intra.mode_chroma == -1) continue;
    kvz_intra_predict(state, refs_u, &loc, COLOR_U, pred, &chroma_data[i], lcu);
    //costs[i] += get_cost(encoder_state, pred, orig_block, satd_func, sad_func, width);
    chroma_data[i].cost += satd_func(pred, orig_block);
  }

  kvz_pixels_blit(orig_v, orig_block, width, width, origstride, width);
  for (int i = 0; i < modes_count; ++i) {
    if (chroma_data[i].pred_cu.intra.mode_chroma == -1) continue;
    kvz_intra_predict(state, refs_v, &loc, COLOR_V, pred, &chroma_data[i], lcu);
    //costs[i] += get_cost(encoder_state, pred, orig_block, satd_func, sad_func, width);
    chroma_data[i].cost += satd_func(pred, orig_block);
  }

  for (int i = 0; i < modes_count; ++i) {
    const double bits = kvz_chroma_mode_bits(state, chroma_data[i].pred_cu.intra.mode_chroma, chroma_data[i].pred_cu.intra.mode);
    chroma_data[i].bits = bits;
    chroma_data[i].cost = bits * state->lambda_sqrt;
  }
}


/**
 * \brief  Order the intra prediction modes according to a fast criteria. 
 *
 * This function uses SATD to order the intra prediction modes. For 4x4 modes
 * SAD might be used instead, if the cost given by SAD is much better than the
 * one given by SATD, to take into account that 4x4 modes can be coded with
 * transform skip. This version of the function calculates two costs
 * simultaneously to better utilize large SIMD registers with AVX and newer
 * extensions.
 *
 * The modes are searched using halving search and the total number of modes
 * that are tried is dependent on size of the predicted block. More modes
 * are tried for smaller blocks.
 *
 * \param orig  Pointer to the top-left corner of current CU in the picture
 *     being encoded.
 * \param orig_stride  Stride of param orig..
 * \param rec  Pointer to the top-left corner of current CU in the picture
 *     being encoded.
 * \param rec_stride  Stride of param rec.
 * \param width  Width of the prediction block.
 * \param intra_preds  Array of the 3 predicted intra modes.
 *
 * \param[out] modes  The modes ordered according to their RD costs, from best
 *     to worst. The number of modes and costs output is given by parameter
 *     modes_to_check.
 * \param[out] costs  The RD costs of corresponding modes in param modes.
 *
 * \return  Number of prediction modes in param modes.
 */
static int16_t search_intra_rough(
  encoder_state_t * const state,
  kvz_pixel *orig,
  int32_t origstride,
  kvz_intra_references *refs,
  int log2_width,
  int8_t *intra_preds,
  intra_search_data_t* modes_out,
  cu_info_t* const pred_cu)
{
  #define PARALLEL_BLKS 2 // TODO: use 4 for AVX-512 in the future?
  assert(log2_width >= 2 && log2_width <= 5);
  int_fast8_t width = 1 << log2_width;
  cost_pixel_nxn_func *satd_func = kvz_pixels_get_satd_func(width);
  cost_pixel_nxn_func *sad_func = kvz_pixels_get_sad_func(width);
  cost_pixel_nxn_multi_func *satd_dual_func = kvz_pixels_get_satd_dual_func(width);
  cost_pixel_nxn_multi_func *sad_dual_func = kvz_pixels_get_sad_dual_func(width);
  int8_t modes[KVZ_NUM_INTRA_MODES];
  double costs[KVZ_NUM_INTRA_MODES];

  // const kvz_config *cfg = &state->encoder_control->cfg;
  // const bool filter_boundary = !(cfg->lossless && cfg->implicit_rdpcm);

  // Temporary block arrays
  kvz_pixel _preds[PARALLEL_BLKS * 32 * 32 + SIMD_ALIGNMENT];
  pred_buffer preds = ALIGNED_POINTER(_preds, SIMD_ALIGNMENT);
  
  kvz_pixel _orig_block[32 * 32 + SIMD_ALIGNMENT];
  kvz_pixel *orig_block = ALIGNED_POINTER(_orig_block, SIMD_ALIGNMENT);

  // Store original block for SAD computation
  kvz_pixels_blit(orig, orig_block, width, width, origstride, width);

  int8_t modes_selected = 0;
  // Note: get_cost and get_cost_dual may return negative costs.
  int32_t min_cost = INT_MAX;
  int32_t max_cost = INT_MIN;
  
  // Initial offset decides how many modes are tried before moving on to the
  // recursive search.
  int offset;
  if (state->encoder_control->cfg.full_intra_search) {
    offset = 1;
  } else {
    static const int8_t offsets[4] = { 2, 4, 8, 8 };
    offset = offsets[log2_width - 2];
  }

  // Calculate SAD for evenly spaced modes to select the starting point for 
  // the recursive search.
  cu_loc_t loc = { 0, 0, width, width, width, width };
  intra_search_data_t search_proxy;
  FILL(search_proxy, 0);
  search_proxy.pred_cu = *pred_cu;

  for (int mode = 2; mode <= 66; mode += PARALLEL_BLKS * offset) {
    
    double costs_out[PARALLEL_BLKS] = { 0 };
    for (int i = 0; i < PARALLEL_BLKS; ++i) {
      if (mode + i * offset <= 66) {
        search_proxy.pred_cu.intra.mode = mode + i*offset;
        kvz_intra_predict(state, refs, &loc, COLOR_Y, preds[i], &search_proxy, NULL);
      }
    }
    
    //TODO: add generic version of get cost  multi
    get_cost_dual(state, preds, orig_block, satd_dual_func, sad_dual_func, width, costs_out);

    for (int i = 0; i < PARALLEL_BLKS; ++i) {
      if (mode + i * offset <= 66) {
        costs[modes_selected] = costs_out[i];
        modes[modes_selected] = mode + i * offset;
        min_cost = MIN(min_cost, costs[modes_selected]);
        max_cost = MAX(max_cost, costs[modes_selected]);
        ++modes_selected;
      }
    }
  }

  int8_t best_mode = modes[select_best_mode_index(modes, costs, modes_selected)];
  double best_cost = min_cost;
  
  // Skip recursive search if all modes have the same cost.
  if (min_cost != max_cost) {
    // Do a recursive search to find the best mode, always centering on the
    // current best mode.
    while (offset > 1) {
      offset >>= 1;

      int8_t center_node = best_mode;
      int8_t test_modes[] = { center_node - offset, center_node + offset };

      double costs_out[PARALLEL_BLKS] = { 0 };
      char mode_in_range = 0;

      for (int i = 0; i < PARALLEL_BLKS; ++i) mode_in_range |= (test_modes[i] >= 2 && test_modes[i] <= 66);

      if (mode_in_range) {
        for (int i = 0; i < PARALLEL_BLKS; ++i) {
          if (test_modes[i] >= 2 && test_modes[i] <= 66) {
            search_proxy.pred_cu.intra.mode = test_modes[i];
            kvz_intra_predict(state, refs, &loc, COLOR_Y, preds[i], &search_proxy, NULL);
          }
        }

        //TODO: add generic version of get cost multi
        get_cost_dual(state, preds, orig_block, satd_dual_func, sad_dual_func, width, costs_out);

        for (int i = 0; i < PARALLEL_BLKS; ++i) {
          if (test_modes[i] >= 2 && test_modes[i] <= 66) {
            costs[modes_selected] = costs_out[i];
            modes[modes_selected] = test_modes[i];
            if (costs[modes_selected] < best_cost) {
              best_cost = costs[modes_selected];
              best_mode = modes[modes_selected];
            }
            ++modes_selected;
          }
        }
      }
    }
  }

  int8_t add_modes[INTRA_MPM_COUNT + 2] = {intra_preds[0], intra_preds[1], intra_preds[2], intra_preds[3], intra_preds[4], intra_preds[5], 0, 1};

  // Add DC, planar and missing predicted modes.
  for (int8_t pred_i = 0; pred_i < (INTRA_MPM_COUNT + 2); ++pred_i) {
    bool has_mode = false;
    int8_t mode = add_modes[pred_i];

    for (int mode_i = 0; mode_i < modes_selected; ++mode_i) {
      if (modes[mode_i] == add_modes[pred_i]) {
        has_mode = true;
        break;
      }
    }

    if (!has_mode) {
      search_proxy.pred_cu.intra.mode = mode;
      kvz_intra_predict(state, refs, &loc, COLOR_Y, preds[0], &search_proxy, NULL);
      costs[modes_selected] = get_cost(state, preds[0], orig_block, satd_func, sad_func, width);
      modes[modes_selected] = mode;
      ++modes_selected;
    }
  }

  // Add prediction mode coding cost as the last thing. We don't want this
  // affecting the halving search.
  const double mpm_mode_bit = CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.intra_luma_mpm_flag_model), 1);
  const double not_mpm_mode_bit = CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.intra_luma_mpm_flag_model), 0);
  const double planar_mode_flag = CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.luma_planar_model[1]), 1);
  const double not_planar_mode_flag = CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.luma_planar_model[1]), 0);
  for (int mode_i = 0; mode_i < modes_selected; ++mode_i) {
    int i = 0;
    int smaller_than_pred = 0;
    double bits;
    for(; i < INTRA_MPM_COUNT; i++) {
      if (intra_preds[i] == modes[mode_i]) {
        break;
      }
      if(modes[mode_i] > intra_preds[i]) {
        smaller_than_pred += 1;
      }
    }
    if (i == 0) {
      bits = planar_mode_flag + mpm_mode_bit;
    }
    else if (i < INTRA_MPM_COUNT) {
      bits = not_planar_mode_flag + mpm_mode_bit + MAX(i, 3);
    }
    else {
      bits = not_mpm_mode_bit + 5 + (mode_i - smaller_than_pred > 3);
    }
    costs[mode_i] += state->lambda_sqrt * bits;
    modes_out[mode_i].cost = costs[mode_i];
    modes_out[mode_i].pred_cu = *pred_cu;
    modes_out[mode_i].pred_cu.intra.mode = modes[mode_i];
    modes_out[mode_i].pred_cu.intra.mode_chroma = modes[mode_i];
  }

  #undef PARALLEL_BLKS
  return modes_selected;
}


static void get_rough_cost_for_n_modes(
  encoder_state_t* const state,
  kvz_intra_references* refs,
  const cu_loc_t* const cu_loc,
  kvz_pixel *orig,
  int orig_stride,
  intra_search_data_t *search_data,
  int num_modes)
{
#define PARALLEL_BLKS 2
  assert(num_modes % 2 == 0 && "passing odd number of modes to get_rough_cost_for_n_modes");
  const int width = cu_loc->width;
  cost_pixel_nxn_multi_func* satd_dual_func = kvz_pixels_get_satd_dual_func(width);
  cost_pixel_nxn_multi_func* sad_dual_func = kvz_pixels_get_sad_dual_func(width);

  kvz_pixel _preds[PARALLEL_BLKS * MIN(LCU_WIDTH, 64)* MIN(LCU_WIDTH, 64)+ SIMD_ALIGNMENT];
  pred_buffer preds = ALIGNED_POINTER(_preds, SIMD_ALIGNMENT);

  kvz_pixel _orig_block[MIN(LCU_WIDTH, 64) * MIN(LCU_WIDTH, 64) + SIMD_ALIGNMENT];
  kvz_pixel* orig_block = ALIGNED_POINTER(_orig_block, SIMD_ALIGNMENT);

  kvz_pixels_blit(orig, orig_block, width, width, orig_stride, width);
  
  double costs_out[PARALLEL_BLKS] = { 0 };
  for(int mode = 0; mode < num_modes; mode += PARALLEL_BLKS) {
    for (int i = 0; i < PARALLEL_BLKS; ++i) {
      kvz_intra_predict(state, refs, cu_loc, COLOR_Y, preds[i], &search_data[mode + i], NULL);
    }
    get_cost_dual(state, preds, orig_block, satd_dual_func, sad_dual_func, width, costs_out);
    search_data[mode].cost = costs_out[0];
    search_data[mode + 1].cost = costs_out[1];
  }
#undef PARALLEL_BLKS
}


/**
 * \brief  Find best intra mode out of the ones listed in parameter modes.
 *
 * This function perform intra search by doing full quantization,
 * reconstruction and CABAC coding of coefficients. It is very slow
 * but results in better RD quality than using just the rough search.
 *
 * \param x_px  Luma picture coordinate.
 * \param y_px  Luma picture coordinate.
 * \param orig  Pointer to the top-left corner of current CU in the picture
 *     being encoded.
 * \param orig_stride  Stride of param orig.
 * \param rec  Pointer to the top-left corner of current CU in the picture
 *     being encoded.
 * \param rec_stride  Stride of param rec.
 * \param intra_preds  Array of the 3 predicted intra modes.
 * \param modes_to_check  How many of the modes in param modes are checked.
 * \param[in] modes  The intra prediction modes that are to be checked.
 * 
 * \param[out] modes  The modes ordered according to their RD costs, from best
 *     to worst. The number of modes and costs output is given by parameter
 *     modes_to_check.
 * \param[out] costs  The RD costs of corresponding modes in param modes.
 * \param[out] lcu  If transform split searching is used, the transform split
 *     information for the best mode is saved in lcu.cu structure.
 */
static int8_t search_intra_rdo(
  encoder_state_t * const state,
  int x_px,
  int y_px,
  int depth,
  int modes_to_check,
  intra_search_data_t *search_data,
  lcu_t *lcu)
{
  const int tr_depth = CLIP(1, MAX_PU_DEPTH, depth + state->encoder_control->cfg.tr_depth_intra);
  
  for (int mode = 0; mode < modes_to_check; mode++) {
    double rdo_bitcost = kvz_luma_mode_bits(state, &search_data[mode].pred_cu, x_px, y_px, depth, lcu);
    search_data[mode].bits = rdo_bitcost;
    search_data[mode].cost = rdo_bitcost * state->lambda;

    double mode_cost = search_intra_trdepth(state, x_px, y_px, depth, tr_depth, MAX_INT, &search_data[mode], lcu);
    search_data[mode].cost += mode_cost;
    if (state->encoder_control->cfg.intra_rdo_et && !cbf_is_set_any(search_data[mode].pred_cu.cbf, depth)) {
      modes_to_check = mode + 1;
      break;
    }
  }

  // Update order according to new costs
  double best_cost = MAX_INT;
  int best_mode = 0;
  for (int mode = 0; mode < modes_to_check; mode++) {
    if(search_data[mode].cost < best_cost) {
      best_cost = search_data[mode].cost;
      best_mode = mode;
    }
  }
  search_data[0] = search_data[best_mode];

  return modes_to_check;
}


double kvz_luma_mode_bits(const encoder_state_t *state, const cu_info_t* const cur_cu, int x, int y, int8_t depth, const lcu_t* lcu)
{
  cabac_data_t* cabac = (cabac_data_t *)&state->search_cabac;
  double mode_bits = 0;
  cabac_data_t cabac_copy;
  memcpy(&cabac_copy, cabac, sizeof cabac_copy);
  kvz_encode_intra_luma_coding_unit(
    state,
    &cabac_copy, cur_cu,
    x, y, depth, lcu, &mode_bits
  );

  return mode_bits;
}


double kvz_chroma_mode_bits(const encoder_state_t *state, int8_t chroma_mode, int8_t luma_mode)
{
  cabac_data_t* cabac = (cabac_data_t*)&state->search_cabac;
  const cabac_ctx_t *ctx = &(cabac->ctx.chroma_pred_model);
  double mode_bits;
  if (chroma_mode == luma_mode) {
    mode_bits = CTX_ENTROPY_FBITS(ctx, 0);
  } else {
    if(chroma_mode > 67) {
      mode_bits = 2.0 + CTX_ENTROPY_FBITS(ctx, 1);
    }
    else {
      ctx = &(state->cabac.ctx.cclm_model);
      mode_bits = CTX_ENTROPY_FBITS(ctx, chroma_mode != 81);
      if (chroma_mode != 81) mode_bits += 1;
    }
  }
  // Technically this is encoded first but for this method of counting bits it does not matter
  if(state->encoder_control->cfg.cclm) {
    ctx = &(state->cabac.ctx.cclm_flag);
    mode_bits += CTX_ENTROPY_FBITS(ctx, chroma_mode > 67);
  }

  if(cabac->update) {
    if(chroma_mode != luma_mode) {
      // Again it does not matter what we actually write here
      CABAC_BINS_EP(cabac, 0, 2, "intra_chroma_pred_mode");      
    }
  }

  return mode_bits;
}


int8_t kvz_search_intra_chroma_rdo(
  encoder_state_t * const state,
  int x_px,
  int y_px,
  int depth,
  int8_t num_modes,
  lcu_t *const lcu,
  intra_search_data_t* chroma_data)
{
  const bool reconstruct_chroma = (depth != 4) || (x_px & 4 && y_px & 4);


  kvz_intra_references refs[2];
  const vector2d_t luma_px = { x_px & ~7, y_px & ~7 };
  const vector2d_t pic_px = {
    state->tile->frame->width,
    state->tile->frame->height,
  };


  if (reconstruct_chroma) {
    kvz_intra_build_reference(MAX(LOG2_LCU_WIDTH - depth - 1, 2), COLOR_U, &luma_px, &pic_px, lcu, &refs[0], state->encoder_control->cfg.wpp, NULL, 0);
    kvz_intra_build_reference(MAX(LOG2_LCU_WIDTH - depth - 1, 2), COLOR_V, &luma_px, &pic_px, lcu, &refs[1], state->encoder_control->cfg.wpp, NULL, 0);
    
    const vector2d_t lcu_px = { SUB_SCU(x_px), SUB_SCU(y_px) };
    cu_info_t *const tr_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x, lcu_px.y);
    
    for (int8_t i = 0; i < num_modes; ++i) {
      const uint8_t mode = chroma_data[i].pred_cu.intra.mode_chroma;
      if(mode < 67 || depth == 0) {
        kvz_intra_recon_cu(state,
                           x_px, y_px,
                           depth, &chroma_data[i],
                           NULL,
                           lcu);
      }
      
      if(tr_cu->depth != tr_cu->tr_depth || !state->encoder_control->cfg.jccr) {
        chroma_data[i].cost = kvz_cu_rd_cost_chroma(state, lcu_px.x, lcu_px.y, depth, tr_cu, lcu);
      } else {
        kvz_select_jccr_mode(state, lcu_px.x, lcu_px.y, depth, tr_cu, lcu, &chroma_data[i].cost);
      }

      double mode_bits = kvz_chroma_mode_bits(state, mode, chroma_data[i].pred_cu.intra.mode);
      chroma_data[i].cost += mode_bits * state->lambda;
    }
    sort_modes(chroma_data, num_modes);

    return chroma_data[0].pred_cu.intra.mode_chroma;
  }

  return 100;
}


int8_t kvz_search_cu_intra_chroma(encoder_state_t * const state,
                              const int x_px, const int y_px,
                              const int depth, lcu_t *lcu, cclm_parameters_t *best_cclm)
{
  const vector2d_t lcu_px = { SUB_SCU(x_px), SUB_SCU(y_px) };

  cu_info_t *cur_pu = LCU_GET_CU_AT_PX(lcu, lcu_px.x, lcu_px.y);
  int8_t intra_mode = cur_pu->intra.mode;
  
  int8_t modes[8] = { 0, 50, 18, 1, -1, 81, 82, 83 };
  uint8_t total_modes = (state->encoder_control->cfg.cclm ? 8 : 5);
  if (intra_mode != 0 && intra_mode != 50 && intra_mode != 18 && intra_mode != 1) {
    modes[4] = intra_mode;
  }
  else {
    total_modes -= 1;
    modes[4] = modes[5];
    modes[5] = modes[6];
    modes[6] = modes[7];
  }


  // The number of modes to select for slower chroma search. Luma mode
  // is always one of the modes, so 2 means the final decision is made
  // between luma mode and one other mode that looks the best
  // according to search_intra_chroma_rough.
  const int8_t modes_in_depth[5] = { 1, 1, 1, 1, 2 };
  int num_modes = modes_in_depth[depth];

  if (state->encoder_control->cfg.rdo >= 3) {
    num_modes = total_modes;
  }

  intra_search_data_t chroma_data[8];
  FILL(chroma_data, 0);
  for (int i = 0; i < num_modes; i++) {
    chroma_data[i].pred_cu = *cur_pu;
    chroma_data[i].pred_cu.intra.mode_chroma = modes[i];
  }
  // Don't do rough mode search if all modes are selected.
  // FIXME: It might make more sense to only disable rough search if
  // num_modes is 0.is 0.

  if (total_modes != num_modes) {
    const int_fast8_t log2_width_c = MAX(LOG2_LCU_WIDTH - depth - 1, 2);
    const vector2d_t pic_px = { state->tile->frame->width, state->tile->frame->height };
    const vector2d_t luma_px = { x_px, y_px };

    kvz_intra_references refs_u;
    kvz_intra_build_reference(log2_width_c, COLOR_U, &luma_px, &pic_px, lcu, &refs_u, state->encoder_control->cfg.wpp, NULL, 0);

    kvz_intra_references refs_v;
    kvz_intra_build_reference(log2_width_c, COLOR_V, &luma_px, &pic_px, lcu, &refs_v, state->encoder_control->cfg.wpp, NULL, 0);

    vector2d_t lcu_cpx = { lcu_px.x / 2, lcu_px.y / 2 };
    kvz_pixel *ref_u = &lcu->ref.u[lcu_cpx.x + lcu_cpx.y * LCU_WIDTH_C];
    kvz_pixel *ref_v = &lcu->ref.v[lcu_cpx.x + lcu_cpx.y * LCU_WIDTH_C];

    search_intra_chroma_rough(state, x_px, y_px, depth,
                              ref_u, ref_v,
                              LCU_WIDTH_C,
                              &refs_u, &refs_v,
      chroma_data, lcu);
    sort_modes(chroma_data, total_modes);
  }

  int8_t intra_mode_chroma = intra_mode;
  if (num_modes > 1) {
    intra_mode_chroma = kvz_search_intra_chroma_rdo(state, x_px, y_px, depth, num_modes, lcu, chroma_data);
  }

  return intra_mode_chroma;
}


/**
 * Update lcu to have best modes at this depth.
 * \return Cost of best mode.
 */
void kvz_search_cu_intra(
  encoder_state_t * const state,
  const int x_px,
  const int y_px,
  const int depth,
  intra_search_data_t* mode_out,
  lcu_t *lcu)
{
  const vector2d_t lcu_px = { SUB_SCU(x_px), SUB_SCU(y_px) };
  const int8_t cu_width = LCU_WIDTH >> depth;
  const cu_loc_t cu_loc = { x_px, y_px, cu_width, cu_width,
    MAX(cu_width >> 1, TR_MIN_WIDTH), MAX(cu_width >> 1, TR_MIN_WIDTH) };
  const int_fast8_t log2_width = LOG2_LCU_WIDTH - depth;

  cu_info_t *cur_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x, lcu_px.y);

  kvz_intra_references refs;

  int8_t candidate_modes[INTRA_MPM_COUNT];
  // Normal intra modes + mrl modes + mip modes
  intra_search_data_t search_data[KVZ_NUM_INTRA_MODES +(MAX_REF_LINE_IDX - 1) * (INTRA_MPM_COUNT - 1) + 32];

  cu_info_t *left_cu = 0;
  cu_info_t *above_cu = 0;

  // Select left and top CUs if they are available.
  // Top CU is not available across LCU boundary.
  if (x_px >= SCU_WIDTH) {
    left_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x - 1, lcu_px.y+ cu_width-1);
  }
  if (y_px >= SCU_WIDTH && lcu_px.y > 0) {
    above_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x+ cu_width-1, lcu_px.y - 1);
  }
  kvz_intra_get_dir_luma_predictor(x_px, y_px, candidate_modes, cur_cu, left_cu, above_cu);

  if (depth > 0) {
    const vector2d_t luma_px = { x_px, y_px };
    const vector2d_t pic_px = { state->tile->frame->width, state->tile->frame->height };

    // These references will only be used with rough search. No need for MRL stuff here.
    kvz_intra_build_reference(log2_width, COLOR_Y, &luma_px, &pic_px, lcu, &refs, state->encoder_control->cfg.wpp, NULL, 0);
  }

  // The maximum number of possible MIP modes depend on block size & shape
  int width = LCU_WIDTH >> depth;
  int height = width; // TODO: proper height for non-square blocks.

  // Find best intra mode for 2Nx2N.
  kvz_pixel *ref_pixels = &lcu->ref.y[lcu_px.x + lcu_px.y * LCU_WIDTH];

  // Need to set some data for all cus
  cu_info_t temp_pred_cu;
  temp_pred_cu = *cur_cu;
  temp_pred_cu.type = CU_INTRA;
  FILL(temp_pred_cu.intra, 0);

  int16_t number_of_modes;
  bool skip_rough_search = (depth == 0 || state->encoder_control->cfg.rdo >= 4);
  if (!skip_rough_search) {
    number_of_modes = search_intra_rough(state,
                                         ref_pixels,
                                         LCU_WIDTH,
                                         &refs,
                                         log2_width, candidate_modes,
                                         search_data, &temp_pred_cu);

  } else {
    for (int8_t i = 0; i < KVZ_NUM_INTRA_MODES; i++) {
      search_data[i].pred_cu = temp_pred_cu;
      search_data[i].pred_cu.intra.mode = i;
      search_data[i].pred_cu.intra.mode_chroma = i;
      search_data[i].cost = MAX_INT;
    }
    number_of_modes = KVZ_NUM_INTRA_MODES;
  }

  int num_mip_modes = 0;
  if (state->encoder_control->cfg.mip) {
    // MIP is not allowed for 64 x 4 or 4 x 64 blocks
    if (!((width == 64 && height == 4) || (width == 4 && height == 64))) {
      num_mip_modes = NUM_MIP_MODES_FULL(width, height);

      for (int transpose = 0; transpose < 2; transpose++) {
        const int half_mip_modes = NUM_MIP_MODES_HALF(width, height);
        for (int i = 0; i < half_mip_modes; ++i) {
          const int index = i + number_of_modes + transpose * half_mip_modes;
          search_data[index].pred_cu = temp_pred_cu;
          search_data[index].pred_cu.intra.mip_flag = 1;
          search_data[index].pred_cu.intra.mode = i;
          search_data[index].pred_cu.intra.mip_is_transposed = transpose;
          search_data[index].pred_cu.intra.mode_chroma = i;
          search_data[index].cost = MAX_INT;
        }
      }
      if(!skip_rough_search) {
        get_rough_cost_for_n_modes(state, &refs, &cu_loc,
          ref_pixels, LCU_WIDTH, search_data + number_of_modes, num_mip_modes);
      }
    }
    number_of_modes += num_mip_modes;
  }

  int num_mrl_modes = 0;
  // Find modes with multiple reference lines if in use. Do not use if CU in first row.
  uint8_t lines = state->encoder_control->cfg.mrl && (y_px % LCU_WIDTH) != 0 ? MAX_REF_LINE_IDX : 1;

  for(int line = 1; line < lines; ++line) {
    for(int i = 1; i < INTRA_MPM_COUNT; i++) {
      num_mrl_modes++;
      const int index = (i - 1) + (INTRA_MPM_COUNT -1)*(line-1) + number_of_modes;
      search_data[index].pred_cu = temp_pred_cu;
      search_data[index].pred_cu.intra.mode = candidate_modes[i];
      search_data[index].pred_cu.intra.multi_ref_idx = line;
      search_data[index].pred_cu.intra.mode_chroma = candidate_modes[i];
      search_data[index].cost = MAX_INT;
    }
  }
  if (!skip_rough_search && lines != 1) {
    get_rough_cost_for_n_modes(state, &refs, &cu_loc,
      ref_pixels, LCU_WIDTH, search_data + number_of_modes, num_mrl_modes);
  }
  number_of_modes += num_mrl_modes;

  // Set transform depth to current depth, meaning no transform splits.
  kvz_lcu_fill_trdepth(lcu, x_px, y_px, depth, depth);
  // Refine results with slower search or get some results if rough search was skipped.
  const int32_t rdo_level = state->encoder_control->cfg.rdo;
  if (rdo_level >= 2 || skip_rough_search) {
    int number_of_modes_to_search;
    if (rdo_level == 4) {
      number_of_modes_to_search = number_of_modes;
    } else if (rdo_level == 2 || rdo_level == 3) {
      number_of_modes_to_search = (cu_width == 4) ? 3 : 2;
    } else {
      // Check only the predicted modes.
      number_of_modes_to_search = 0;
    }
    if(!skip_rough_search) {
      sort_modes(search_data, number_of_modes);
    }

    for(int pred_mode = 0; pred_mode < INTRA_MPM_COUNT; ++pred_mode) {
      bool mode_found = false;
      for(int i = 0; i < number_of_modes_to_search; i++) {
        if(search_data[i].pred_cu.intra.mode == candidate_modes[pred_mode]) {
          mode_found = true;
          break;
        }
      }
      if(!mode_found) {
        search_data[number_of_modes_to_search].pred_cu = temp_pred_cu;
        search_data[number_of_modes_to_search].pred_cu.intra.mode = candidate_modes[pred_mode];
        search_data[number_of_modes_to_search].pred_cu.intra.mode_chroma = candidate_modes[pred_mode];
        number_of_modes_to_search++;
      }
    }

    // TODO: if rough search is implemented for MIP, sort mip_modes here.
    search_intra_rdo(
      state,
      x_px,
      y_px,
      depth,
      number_of_modes_to_search,
      search_data,
      lcu);
    
  }
  
  *mode_out = search_data[0];
}
