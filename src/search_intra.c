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
#include "image.h"
#include "intra.h"
#include "uvg266.h"
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
                       uvg_pixel *pred, uvg_pixel *orig_block,
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
    const cabac_ctx_t *ctx = &state->cabac.ctx.transform_skip_model_luma;
    double trskip_bits = CTX_ENTROPY_FBITS(ctx, 1) - CTX_ENTROPY_FBITS(ctx, 0);

    
    // ToDo: Check cost
    if (state->encoder_control->chroma_format != UVG_CSP_400) {
      ctx = &state->cabac.ctx.transform_skip_model_chroma;
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
                       const pred_buffer preds, const uvg_pixel *orig_block,
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
    if (state->encoder_control->chroma_format != UVG_CSP_400) {
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
  int8_t scan_idx = uvg_get_scan_order(pred_cu->type, pred_cu->intra.mode, depth);
  int32_t i;
  // ToDo: large block support in VVC?
  uint32_t sig_coeffgroup_flag[32 * 32] = { 0 };

  const uint32_t log2_block_size = uvg_g_convert_to_bit[width] + 2;
  const uint32_t log2_cg_size = uvg_g_log2_sbb_size[log2_block_size][log2_block_size][0]
    + uvg_g_log2_sbb_size[log2_block_size][log2_block_size][1];
  const uint32_t *scan = uvg_g_sig_last_scan[scan_idx][log2_block_size - 1];
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
static double search_intra_trdepth(encoder_state_t * const state,
                                   int x_px, int y_px, int depth, int max_depth,
                                   int intra_mode, int cost_treshold,
                                   cu_info_t *const pred_cu,
                                   lcu_t *const lcu,
                                   cclm_parameters_t *cclm_params,
                                   const int mts_mode)
{
  assert(depth >= 0 && depth <= MAX_PU_DEPTH);

  const int width = LCU_WIDTH >> depth;
  const int width_c = width > TR_MIN_WIDTH ? width / 2 : width;

  const int offset = width / 2;
  const vector2d_t lcu_px = { SUB_SCU(x_px), SUB_SCU(y_px) };
  cu_info_t *const tr_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x, lcu_px.y);

  const bool reconstruct_chroma = (depth != 4 || (depth == 4 && (x_px & 4 && y_px & 4))) && state->encoder_control->chroma_format != UVG_CSP_400;

  struct {
    uvg_pixel y[TR_MAX_WIDTH*TR_MAX_WIDTH];
    uvg_pixel u[TR_MAX_WIDTH*TR_MAX_WIDTH];
    uvg_pixel v[TR_MAX_WIDTH*TR_MAX_WIDTH];
  } nosplit_pixels;
  uint16_t nosplit_cbf = 0;

  double split_cost = INT32_MAX;
  double nosplit_cost = INT32_MAX;

  if (depth > 0) {
    const bool mts_enabled = state->encoder_control->cfg.mts == UVG_MTS_INTRA || state->encoder_control->cfg.mts == UVG_MTS_BOTH;
    tr_cu->tr_depth = depth;
    pred_cu->tr_depth = depth;

    nosplit_cost = 0.0;

    cbf_clear(&pred_cu->cbf, depth, COLOR_Y);
    if (reconstruct_chroma) {
      cbf_clear(&pred_cu->cbf, depth, COLOR_U);
      cbf_clear(&pred_cu->cbf, depth, COLOR_V);
    }

    const int8_t chroma_mode = reconstruct_chroma ? intra_mode : -1;
    double best_rd_cost = MAX_INT;
    int best_tr_idx = 0;

    int trafo;
    int num_transforms = 1;
    if (mts_mode != -1)
    {
      trafo = mts_mode;
      num_transforms = mts_mode + 1;
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
     
      uvg_intra_recon_cu(state,
        x_px, y_px,
        depth,
        intra_mode, -1,
        pred_cu, cclm_params, pred_cu->intra.multi_ref_idx, 
        pred_cu->intra.mip_flag, pred_cu->intra.mip_is_transposed,
        lcu);

      // TODO: Not sure if this should be 0 or 1 but at least seems to work with 1
      if (pred_cu->tr_idx > 1)
      {
        derive_mts_constraints(pred_cu, lcu, depth, lcu_px);
        if (pred_cu->violates_mts_coeff_constraint || !pred_cu->mts_last_scan_pos)
        {
          assert(mts_mode == -1); //mts mode should not be decided and then not allowed to be used. (might be some exception here)
          continue;
        }
      }

      double rd_cost = uvg_cu_rd_cost_luma(state, lcu_px.x, lcu_px.y, depth, pred_cu, lcu);
      //if (reconstruct_chroma) {
      //  rd_cost += uvg_cu_rd_cost_chroma(state, lcu_px.x, lcu_px.y, depth, pred_cu, lcu);
      //}

      if (rd_cost < best_rd_cost) {
        best_rd_cost = rd_cost;
        best_tr_idx = pred_cu->tr_idx;
      }
    }
    if(reconstruct_chroma) {
      uvg_intra_recon_cu(state,
        x_px, y_px,
        depth,
        -1, chroma_mode,
        pred_cu, cclm_params, 0, 
        pred_cu->intra.mip_flag, pred_cu->intra.mip_is_transposed,
        lcu);
      best_rd_cost += uvg_cu_rd_cost_chroma(state, lcu_px.x, lcu_px.y, depth, pred_cu, lcu);
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

    uvg_pixels_blit(lcu->rec.y, nosplit_pixels.y, width, width, LCU_WIDTH, width);
    if (reconstruct_chroma) {
      uvg_pixels_blit(lcu->rec.u, nosplit_pixels.u, width_c, width_c, LCU_WIDTH_C, width_c);
      uvg_pixels_blit(lcu->rec.v, nosplit_pixels.v, width_c, width_c, LCU_WIDTH_C, width_c);
    }
  }

  // Recurse further if all of the following:
  // - Current depth is less than maximum depth of the search (max_depth).
  //   - Maximum transform hierarchy depth is constrained by clipping
  //     max_depth.
  // - Min transform size hasn't been reached (MAX_PU_DEPTH).
  if (depth < max_depth && depth < MAX_PU_DEPTH) {
    split_cost = 3 * state->lambda;

    split_cost += search_intra_trdepth(state, x_px, y_px, depth + 1, max_depth, intra_mode, nosplit_cost, pred_cu, lcu, cclm_params, -1);
    if (split_cost < nosplit_cost) {
      split_cost += search_intra_trdepth(state, x_px + offset, y_px, depth + 1, max_depth, intra_mode, nosplit_cost, pred_cu, lcu, cclm_params, -1);
    }
    if (split_cost < nosplit_cost) {
      split_cost += search_intra_trdepth(state, x_px, y_px + offset, depth + 1, max_depth, intra_mode, nosplit_cost, pred_cu, lcu, cclm_params, -1);
    }
    if (split_cost < nosplit_cost) {
      split_cost += search_intra_trdepth(state, x_px + offset, y_px + offset, depth + 1, max_depth, intra_mode, nosplit_cost, pred_cu, lcu, cclm_params, -1);
    }

    double cbf_bits = 0.0;

    // Add cost of cbf chroma bits on transform tree.
    // All cbf bits are accumulated to pred_cu.cbf and cbf_is_set returns true
    // if cbf is set at any level >= depth, so cbf chroma is assumed to be 0
    // if this and any previous transform block has no chroma coefficients.
    // When searching the first block we don't actually know the real values,
    // so this will code cbf as 0 and not code the cbf at all for descendants.
    if (state->encoder_control->chroma_format != UVG_CSP_400) {
      const uint8_t tr_depth = depth - pred_cu->depth;

      const cabac_ctx_t* ctx = &(state->cabac.ctx.qt_cbf_model_cb[0]);
      if (tr_depth == 0 || cbf_is_set(pred_cu->cbf, depth - 1, COLOR_U)) {
        cbf_bits += CTX_ENTROPY_FBITS(ctx, cbf_is_set(pred_cu->cbf, depth, COLOR_U));
      }
      ctx = &(state->cabac.ctx.qt_cbf_model_cr[cbf_is_set(pred_cu->cbf, depth, COLOR_U)]);
      if (tr_depth == 0 || cbf_is_set(pred_cu->cbf, depth - 1, COLOR_V)) {
        cbf_bits += CTX_ENTROPY_FBITS(ctx, cbf_is_set(pred_cu->cbf, depth, COLOR_V));
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
    uvg_lcu_fill_trdepth(lcu, x_px, y_px, depth, depth);

    pred_cu->cbf = nosplit_cbf;

    // We only restore the pixel data and not coefficients or cbf data.
    // The only thing we really need are the border pixels.uvg_intra_get_dir_luma_predictor
    uvg_pixels_blit(nosplit_pixels.y, lcu->rec.y, width, width, width, LCU_WIDTH);
    if (reconstruct_chroma) {
      uvg_pixels_blit(nosplit_pixels.u, lcu->rec.u, width_c, width_c, width_c, LCU_WIDTH_C);
      uvg_pixels_blit(nosplit_pixels.v, lcu->rec.v, width_c, width_c, width_c, LCU_WIDTH_C);
    }

    return nosplit_cost;
  }
}


static void search_intra_chroma_rough(encoder_state_t * const state,
                                      int x_px, int y_px, int depth,
                                      const uvg_pixel *orig_u, const uvg_pixel *orig_v, int16_t origstride,
                                      uvg_intra_references *refs_u, uvg_intra_references *refs_v,
                                      int8_t luma_mode,
                                      int8_t modes[8], double costs[8], lcu_t* lcu)
{
  assert(!(x_px & 4 || y_px & 4));

  const unsigned width = MAX(LCU_WIDTH_C >> depth, TR_MIN_WIDTH);
  const int_fast8_t log2_width_c = MAX(LOG2_LCU_WIDTH - (depth + 1), 2);

  for (int i = 0; i < 8; ++i) {
    costs[i] = 0;
  }

  cost_pixel_nxn_func *const satd_func = uvg_pixels_get_satd_func(width);
  //cost_pixel_nxn_func *const sad_func = uvg_pixels_get_sad_func(width);

  cclm_parameters_t cclm_params;
  
  uvg_pixel _pred[32 * 32 + SIMD_ALIGNMENT];
  uvg_pixel *pred = ALIGNED_POINTER(_pred, SIMD_ALIGNMENT);

  uvg_pixel _orig_block[32 * 32 + SIMD_ALIGNMENT];
  uvg_pixel *orig_block = ALIGNED_POINTER(_orig_block, SIMD_ALIGNMENT);

  uvg_pixels_blit(orig_u, orig_block, width, width, origstride, width);
  for (int i = 0; i < 5; ++i) {
    if (modes[i] == -1) continue;
    uvg_intra_predict(state, refs_u, log2_width_c, modes[i], COLOR_U, pred, false, 0);
    //costs[i] += get_cost(encoder_state, pred, orig_block, satd_func, sad_func, width);
    costs[i] += satd_func(pred, orig_block);
  }
  for (int i = 5; i < 8; i++) {
    assert(state->encoder_control->cfg.cclm);
    uvg_predict_cclm(
      state,
      COLOR_U, width, width, x_px, y_px, state->tile->frame->source->stride, modes[i], lcu, refs_u,  pred, &cclm_params);
  }

  uvg_pixels_blit(orig_v, orig_block, width, width, origstride, width);
  for (int i = 0; i < 5; ++i) {
    if (modes[i] == -1) continue;
    uvg_intra_predict(state, refs_v, log2_width_c, modes[i], COLOR_V, pred, false, 0);
    //costs[i] += get_cost(encoder_state, pred, orig_block, satd_func, sad_func, width);
    costs[i] += satd_func(pred, orig_block);
  }
  for (int i = 5; i < 8; i++) {
    assert(state->encoder_control->cfg.cclm);
    uvg_predict_cclm(
      state,
      COLOR_V, width, width, x_px, y_px, state->tile->frame->source->stride, modes[i], lcu, refs_u, pred, &cclm_params);
  }

  uvg_sort_modes(modes, costs, 5);
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
static int8_t search_intra_rough(encoder_state_t * const state, 
                                 uvg_pixel *orig, int32_t origstride,
                                 uvg_intra_references *refs,
                                 int log2_width, int8_t *intra_preds,
                                 int8_t modes[67], double costs[67])
{
  #define PARALLEL_BLKS 2 // TODO: use 4 for AVX-512 in the future?
  assert(log2_width >= 2 && log2_width <= 5);
  int_fast8_t width = 1 << log2_width;
  cost_pixel_nxn_func *satd_func = uvg_pixels_get_satd_func(width);
  cost_pixel_nxn_func *sad_func = uvg_pixels_get_sad_func(width);
  cost_pixel_nxn_multi_func *satd_dual_func = uvg_pixels_get_satd_dual_func(width);
  cost_pixel_nxn_multi_func *sad_dual_func = uvg_pixels_get_sad_dual_func(width);

  const uvg_config *cfg = &state->encoder_control->cfg;
  const bool filter_boundary = !(cfg->lossless && cfg->implicit_rdpcm);

  // Temporary block arrays
  uvg_pixel _preds[PARALLEL_BLKS * 32 * 32 + SIMD_ALIGNMENT];
  pred_buffer preds = ALIGNED_POINTER(_preds, SIMD_ALIGNMENT);
  
  uvg_pixel _orig_block[32 * 32 + SIMD_ALIGNMENT];
  uvg_pixel *orig_block = ALIGNED_POINTER(_orig_block, SIMD_ALIGNMENT);

  // Store original block for SAD computation
  uvg_pixels_blit(orig, orig_block, width, width, origstride, width);

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
  for (int mode = 2; mode <= 66; mode += PARALLEL_BLKS * offset) {
    
    double costs_out[PARALLEL_BLKS] = { 0 };
    for (int i = 0; i < PARALLEL_BLKS; ++i) {
      if (mode + i * offset <= 66) {
        uvg_intra_predict(state, refs, log2_width, mode + i * offset, COLOR_Y, preds[i], filter_boundary, 0);
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
            uvg_intra_predict(state, refs, log2_width, test_modes[i], COLOR_Y, preds[i], filter_boundary, 0);
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

  int8_t add_modes[5] = {intra_preds[0], intra_preds[1], intra_preds[2], 0, 1};

  // Add DC, planar and missing predicted modes.
  for (int8_t pred_i = 0; pred_i < 5; ++pred_i) {
    bool has_mode = false;
    int8_t mode = add_modes[pred_i];

    for (int mode_i = 0; mode_i < modes_selected; ++mode_i) {
      if (modes[mode_i] == add_modes[pred_i]) {
        has_mode = true;
        break;
      }
    }

    if (!has_mode) {
      uvg_intra_predict(state, refs, log2_width, mode, COLOR_Y, preds[0], filter_boundary, 0);
      costs[modes_selected] = get_cost(state, preds[0], orig_block, satd_func, sad_func, width);
      modes[modes_selected] = mode;
      ++modes_selected;
    }
  }

  // Add prediction mode coding cost as the last thing. We don't want this
  // affecting the halving search.
  int lambda_cost = (int)(state->lambda_sqrt + 0.5);
  for (int mode_i = 0; mode_i < modes_selected; ++mode_i) {
    costs[mode_i] += lambda_cost * uvg_luma_mode_bits(state, modes[mode_i], intra_preds, 0, 0, 0);
  }

  #undef PARALLEL_BLKS

  return modes_selected;
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
static int8_t search_intra_rdo(encoder_state_t * const state, 
                             int x_px, int y_px, int depth,
                             uvg_pixel *orig, int32_t origstride,
                             int8_t *intra_preds,
                             int modes_to_check,
                             int8_t modes[67], int8_t trafo[67], double costs[67],
                             int num_mip_modes_full,
                             int8_t mip_modes[32], int8_t mip_trafo[32], double mip_costs[32],
                             lcu_t *lcu,
                             uint8_t multi_ref_idx)
{
  const int tr_depth = CLIP(1, MAX_PU_DEPTH, depth + state->encoder_control->cfg.tr_depth_intra);
  const int width = LCU_WIDTH >> depth;
  const int height = width; // TODO: proper height for non-square blocks

  uvg_pixel orig_block[LCU_WIDTH * LCU_WIDTH + 1];

  uvg_pixels_blit(orig, orig_block, width, height, origstride, width);

  // Check that the predicted modes are in the RDO mode list
  if (modes_to_check < 67) {
    int pred_mode = 0;
    // Skip planar if searching modes for MRL
    if (multi_ref_idx != 0) {
      pred_mode = 1;
    }
    for (; pred_mode < 6; pred_mode++) {
      int mode_found = 0;
      for (int rdo_mode = 0; rdo_mode < modes_to_check; rdo_mode++) {
        if (intra_preds[pred_mode] == modes[rdo_mode]) {
          mode_found = 1;
          break;
        }
      }
      // Add this prediction mode to RDO checking
      if (!mode_found) {
        modes[modes_to_check] = intra_preds[pred_mode];
        modes_to_check++;
      }
    }
  }

  // MIP_TODO: implement this inside the standard intra for loop. Code duplication is bad.
  // MIP_TODO: loop through normal intra modes first
  
  for (int mip = 0; mip <= 1; mip++) {
    const int transp_off = mip ? num_mip_modes_full >> 1 : 0;
    uint8_t ctx_id = mip ? uvg_get_mip_flag_context(x_px, y_px, width, height, lcu, NULL) : 0;
    uint8_t multi_ref_index = mip ? 0 : multi_ref_idx;
    int *num_modes = mip ? &num_mip_modes_full : &modes_to_check;

    for (uint8_t i = 0; i < *num_modes; i++) {
      int8_t mode = mip ? mip_modes[i] : modes[i];
      double *mode_cost_p = mip ? &mip_costs[i] : &costs[i];
      int8_t *mode_trafo_p = mip ? &mip_trafo[i] : &trafo[i];
      int rdo_bitcost = uvg_luma_mode_bits(state, mode, intra_preds, multi_ref_index, transp_off, ctx_id);

      *mode_cost_p = rdo_bitcost * (int)(state->lambda + 0.5);

      // Mip related stuff
      // There can be 32 MIP modes, but only mode numbers [0, 15] are ever written to bitstream.
      // Half of the modes [16, 31] are indicated with the separate transpose flag.
      // Number of possible modes is less for larger blocks.
      const bool is_transposed = mip ? (mode >= transp_off ? true : false) : 0;
      int8_t pred_mode = (is_transposed ? mode - transp_off : mode);

      // Perform transform split search and save mode RD cost for the best one.
      cu_info_t pred_cu;
      pred_cu.depth = depth;
      pred_cu.type = CU_INTRA;
      pred_cu.part_size = ((depth == MAX_PU_DEPTH) ? SIZE_NxN : SIZE_2Nx2N); // TODO: non-square blocks
      pred_cu.intra.mode = pred_mode;
      pred_cu.intra.mode_chroma = pred_mode;
      pred_cu.intra.multi_ref_idx = multi_ref_index;
      pred_cu.intra.mip_is_transposed = is_transposed;
      pred_cu.intra.mip_flag = mip ? true : false;
      pred_cu.joint_cb_cr = 0;
      FILL(pred_cu.cbf, 0);

      // Reset transform split data in lcu.cu for this area.
      uvg_lcu_fill_trdepth(lcu, x_px, y_px, depth, depth);

      double mode_cost = search_intra_trdepth(state, x_px, y_px, depth, tr_depth, pred_mode, MAX_INT, &pred_cu, lcu, NULL, -1);
      *mode_cost_p += mode_cost;
      *mode_trafo_p = pred_cu.tr_idx;

      // Early termination if no coefficients has to be coded
      if (state->encoder_control->cfg.intra_rdo_et && !cbf_is_set_any(pred_cu.cbf, depth)) {
        *num_modes = i + 1;
        break;
      }
    }
  }

  // Update order according to new costs
  uvg_sort_modes_intra_luma(modes, trafo, costs, modes_to_check);
  bool use_mip = false;
  if (num_mip_modes_full) {
    uvg_sort_modes_intra_luma(mip_modes, mip_trafo, mip_costs, num_mip_modes_full);
    if (costs[0] > mip_costs[0]) {
      use_mip = true;
    }
  }
  

  // The best transform split hierarchy is not saved anywhere, so to get the
  // transform split hierarchy the search has to be performed again with the
  // best mode.
  if (tr_depth != depth) {
    cu_info_t pred_cu;
    pred_cu.depth = depth;
    pred_cu.type = CU_INTRA;
    pred_cu.part_size = ((depth == MAX_PU_DEPTH) ? SIZE_NxN : SIZE_2Nx2N);
    if (use_mip) {
      int transp_off = num_mip_modes_full >> 1;
      bool is_transposed = (mip_modes[0] >= transp_off ? true : false);
      int8_t pred_mode = (is_transposed ? mip_modes[0] - transp_off : mip_modes[0]);
      pred_cu.intra.mode = pred_mode;
      pred_cu.intra.mode_chroma = pred_mode;
      pred_cu.intra.multi_ref_idx = 0;
      pred_cu.intra.mip_flag = true;
      pred_cu.intra.mip_is_transposed = is_transposed;
    }
    else {
      pred_cu.intra.mode = modes[0];
      pred_cu.intra.mode_chroma = modes[0];
      pred_cu.intra.multi_ref_idx = multi_ref_idx;
      pred_cu.intra.mip_flag = false;
      pred_cu.intra.mip_is_transposed = false;
    }
    FILL(pred_cu.cbf, 0);
    search_intra_trdepth(state, x_px, y_px, depth, tr_depth, pred_cu.intra.mode, MAX_INT, &pred_cu, lcu, NULL, trafo[0]);
  }

  // TODO: modes to check does not consider mip modes. Maybe replace with array when mip search is optimized?
  return modes_to_check;
}


double uvg_luma_mode_bits(const encoder_state_t *state, int8_t luma_mode, const int8_t *intra_preds, const uint8_t multi_ref_idx, const uint8_t num_mip_modes_half, int mip_flag_ctx_id)
{
  double mode_bits = 0.0;

  bool enable_mip = state->encoder_control->cfg.mip;
  bool mip_flag = enable_mip ? (num_mip_modes_half > 0 ? true : false) : false;

  // Mip flag cost must be calculated even if mip is not used in this block
  if (enable_mip) {
    // Make a copy of state->cabac for bit cost estimation.
    cabac_data_t state_cabac_copy;
    cabac_data_t* cabac;
    memcpy(&state_cabac_copy, &state->cabac, sizeof(cabac_data_t));
    // Clear data and set mode to count only
    state_cabac_copy.only_count = 1;
    state_cabac_copy.num_buffered_bytes = 0;
    state_cabac_copy.bits_left = 23;

    cabac = &state_cabac_copy;

    // Do cabac writes as normal
    const int transp_off = num_mip_modes_half;
    const bool is_transposed = luma_mode >= transp_off ? true : false;
    int8_t mip_mode = is_transposed ? luma_mode - transp_off : luma_mode;

    // Write MIP flag
    cabac->cur_ctx = &(cabac->ctx.mip_flag[mip_flag_ctx_id]);
    CABAC_BIN(cabac, mip_flag, "mip_flag");
    
    if (mip_flag) {
      // Write MIP transpose flag & mode
      CABAC_BIN_EP(cabac, is_transposed, "mip_transposed");
      uvg_cabac_encode_trunc_bin(cabac, mip_mode, transp_off);
    }
    
    // Write is done. Get bit cost out of cabac
    mode_bits += (23 - state_cabac_copy.bits_left) + (state_cabac_copy.num_buffered_bytes << 3);
  }

  if (!mip_flag) {
    int8_t mode_in_preds = -1;
    for (int i = 0; i < INTRA_MPM_COUNT; ++i) {
      if (luma_mode == intra_preds[i]) {
        mode_in_preds = i;
        break;
      }
    }

    bool enable_mrl = state->encoder_control->cfg.mrl;
    uint8_t multi_ref_index = enable_mrl ? multi_ref_idx : 0;

    const cabac_ctx_t* ctx = &(state->cabac.ctx.intra_luma_mpm_flag_model);

    if (multi_ref_index == 0) {
      mode_bits += CTX_ENTROPY_FBITS(ctx, mode_in_preds != -1);
    }

    // Add MRL bits.
    if (enable_mrl && MAX_REF_LINE_IDX > 1) {
      ctx = &(state->cabac.ctx.multi_ref_line[0]);
      mode_bits += CTX_ENTROPY_FBITS(ctx, multi_ref_index != 0);

      if (multi_ref_index != 0 && MAX_REF_LINE_IDX > 2) {
        ctx = &(state->cabac.ctx.multi_ref_line[1]);
        mode_bits += CTX_ENTROPY_FBITS(ctx, multi_ref_index != 1);
      }
    }

    if (mode_in_preds != -1 || multi_ref_index != 0) {
      ctx = &(state->cabac.ctx.luma_planar_model[0]);
      if (multi_ref_index == 0) {
        mode_bits += CTX_ENTROPY_FBITS(ctx, mode_in_preds > 0);
      }
      mode_bits += MIN(4.0, mode_in_preds);
    }
    else {
      mode_bits += 6.0;
    }
  }

  return mode_bits;
}


double uvg_chroma_mode_bits(const encoder_state_t *state, int8_t chroma_mode, int8_t luma_mode)
{
  const cabac_ctx_t *ctx = &(state->cabac.ctx.chroma_pred_model);
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

  return mode_bits;
}


int8_t uvg_search_intra_chroma_rdo(encoder_state_t * const state,
                                  int x_px, int y_px, int depth,
                                  int8_t intra_mode,
                                  int8_t modes[8], int8_t num_modes,
                                  lcu_t *const lcu, cclm_parameters_t *best_cclm)
{
  const bool reconstruct_chroma = (depth != 4) || (x_px & 4 && y_px & 4);


  uvg_intra_references refs[2];
  const vector2d_t luma_px = { x_px & ~7, y_px & ~7 };
  const vector2d_t pic_px = {
    state->tile->frame->width,
    state->tile->frame->height,
  };


  if (reconstruct_chroma) {

    int c_width = MAX(32 >> (depth), 4);

    uvg_intra_build_reference(MAX(LOG2_LCU_WIDTH - depth - 1, 2), COLOR_U, &luma_px, &pic_px, lcu, &refs[0], state->encoder_control->cfg.wpp, NULL, 0);
    uvg_intra_build_reference(MAX(LOG2_LCU_WIDTH - depth - 1, 2), COLOR_V, &luma_px, &pic_px, lcu, &refs[1], state->encoder_control->cfg.wpp, NULL, 0);

    cclm_parameters_t cclm_params[2] = { 0 };

    const vector2d_t lcu_px = { SUB_SCU(x_px), SUB_SCU(y_px) };
    cu_info_t *const tr_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x, lcu_px.y);

    struct {
      double cost;
      int8_t mode;
      cclm_parameters_t cclm[2];
    } chroma, best_chroma;

    // chroma.cclm = cclm_params;

    best_chroma.mode = 0;
    best_chroma.cost = MAX_INT;

    for (int8_t chroma_mode_i = 0; chroma_mode_i < num_modes; ++chroma_mode_i) {
      chroma.mode = modes[chroma_mode_i];
      if (chroma.mode == -1) continue;
      if(chroma.mode < 67 || depth == 0) {
        uvg_intra_recon_cu(state,
          x_px, y_px,
          depth,
          -1, chroma.mode, // skip luma
          NULL, NULL, 0, false, false, lcu);
      }
      else {

        uvg_predict_cclm(
          state, COLOR_U,
          c_width, c_width,
          x_px & ~7, y_px & ~7,
          state->tile->frame->source->stride,
          chroma.mode, 
          lcu,
          &refs[0], NULL,
          &cclm_params[0]);

        chroma.cclm[0] = cclm_params[0];

        uvg_predict_cclm(
          state, COLOR_V,
          c_width, c_width,
          x_px & ~7, y_px & ~7,
          state->tile->frame->source->stride, 
          chroma.mode, 
          lcu, 
          &refs[1], NULL,
          &cclm_params[1]);

        chroma.cclm[1] = cclm_params[1];

        uvg_intra_recon_cu(
          state,
          x_px, y_px,
          depth,
          -1, chroma.mode, // skip luma
          NULL, cclm_params, 0, false, false, lcu);
      }
      chroma.cost = uvg_cu_rd_cost_chroma(state, lcu_px.x, lcu_px.y, depth, tr_cu, lcu);

      double mode_bits = uvg_chroma_mode_bits(state, chroma.mode, intra_mode);
      chroma.cost += mode_bits * state->lambda;

      if (chroma.cost < best_chroma.cost) {
        best_chroma = chroma;
      }
    }
    best_cclm[0] = best_chroma.cclm[0];
    best_cclm[1] = best_chroma.cclm[1];

    return best_chroma.mode;
  }

  return 100;
}


int8_t uvg_search_cu_intra_chroma(encoder_state_t * const state,
                              const int x_px, const int y_px,
                              const int depth, lcu_t *lcu, cclm_parameters_t *best_cclm)
{
  const vector2d_t lcu_px = { SUB_SCU(x_px), SUB_SCU(y_px) };

  cu_info_t *cur_pu = LCU_GET_CU_AT_PX(lcu, lcu_px.x, lcu_px.y);
  int8_t intra_mode = cur_pu->intra.mode;

  double costs[8];
  int8_t modes[8] = { 0, 50, 18, 1, -1, 81, 82, 83 };
  if (intra_mode != 0 && intra_mode != 50 && intra_mode != 18 && intra_mode != 1) {
    modes[4] = intra_mode;
  }

  // The number of modes to select for slower chroma search. Luma mode
  // is always one of the modes, so 2 means the final decision is made
  // between luma mode and one other mode that looks the best
  // according to search_intra_chroma_rough.
  const int8_t modes_in_depth[5] = { 1, 1, 1, 1, 2 };
  int num_modes = modes_in_depth[depth];

  if (state->encoder_control->cfg.rdo >= 3) {
    num_modes = state->encoder_control->cfg.cclm ? 8 : 5;
  }

  // Don't do rough mode search if all modes are selected.
  // FIXME: It might make more sense to only disable rough search if
  // num_modes is 0.is 0.
  if (num_modes != 1 && num_modes != 5 && num_modes != 4 && num_modes != 8) {
    const int_fast8_t log2_width_c = MAX(LOG2_LCU_WIDTH - depth - 1, 2);
    const vector2d_t pic_px = { state->tile->frame->width, state->tile->frame->height };
    const vector2d_t luma_px = { x_px, y_px };

    uvg_intra_references refs_u;
    uvg_intra_build_reference(log2_width_c, COLOR_U, &luma_px, &pic_px, lcu, &refs_u, state->encoder_control->cfg.wpp, NULL, 0);

    uvg_intra_references refs_v;
    uvg_intra_build_reference(log2_width_c, COLOR_V, &luma_px, &pic_px, lcu, &refs_v, state->encoder_control->cfg.wpp, NULL, 0);

    vector2d_t lcu_cpx = { lcu_px.x / 2, lcu_px.y / 2 };
    uvg_pixel *ref_u = &lcu->ref.u[lcu_cpx.x + lcu_cpx.y * LCU_WIDTH_C];
    uvg_pixel *ref_v = &lcu->ref.v[lcu_cpx.x + lcu_cpx.y * LCU_WIDTH_C];

    search_intra_chroma_rough(state, x_px, y_px, depth,
                              ref_u, ref_v, LCU_WIDTH_C,
                              &refs_u, &refs_v,
                              intra_mode, modes, costs, lcu);
  }

  int8_t intra_mode_chroma = intra_mode;
  if (num_modes > 1) {
    intra_mode_chroma = uvg_search_intra_chroma_rdo(state, x_px, y_px, depth, intra_mode, modes, num_modes, lcu, best_cclm);
  }

  return intra_mode_chroma;
}


/**
 * Update lcu to have best modes at this depth.
 * \return Cost of best mode.
 */
void uvg_search_cu_intra(encoder_state_t * const state,
                         const int x_px, const int y_px,
                         const int depth, lcu_t *lcu,
                         int8_t *mode_out, 
                         int8_t *trafo_out, 
                         double *cost_out,
                         uint8_t *multi_ref_idx_out,
                         bool *mip_flag_out,
                         bool * mip_transposed_out)
{
  const vector2d_t lcu_px = { SUB_SCU(x_px), SUB_SCU(y_px) };
  const int8_t cu_width = LCU_WIDTH >> depth;
  const int_fast8_t log2_width = LOG2_LCU_WIDTH - depth;

  cu_info_t *cur_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x, lcu_px.y);

  uvg_intra_references refs;

  int8_t candidate_modes[INTRA_MPM_COUNT];

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
  uvg_intra_get_dir_luma_predictor(x_px, y_px, candidate_modes, cur_cu, left_cu, above_cu);

  if (depth > 0) {
    const vector2d_t luma_px = { x_px, y_px };
    const vector2d_t pic_px = { state->tile->frame->width, state->tile->frame->height };

    // These references will only be used with rough search. No need for MRL stuff here.
    uvg_intra_build_reference(log2_width, COLOR_Y, &luma_px, &pic_px, lcu, &refs, state->encoder_control->cfg.wpp, NULL, 0);
  }

  int8_t modes[MAX_REF_LINE_IDX][67];
  int8_t trafo[MAX_REF_LINE_IDX][67] = { 0 };
  double costs[MAX_REF_LINE_IDX][67];

  bool enable_mip = state->encoder_control->cfg.mip;
  // The maximum number of mip modes is 32. Max modes can be less depending on block size.
  // Half of the possible modes are transposed, which is indicated by a separate transpose flag
  int8_t mip_modes[32]; 
  int8_t mip_trafo[32];
  double mip_costs[32];

  // The maximum number of possible MIP modes depend on block size & shape
  int width = LCU_WIDTH >> depth;
  int height = width; // TODO: proper height for non-square blocks.
  int num_mip_modes = 0;

  if (enable_mip) {
    for (int i = 0; i < 32; ++i) {
      mip_modes[i] = i;
      mip_costs[i] = MAX_INT;
    }
    // MIP is not allowed for 64 x 4 or 4 x 64 blocks
    if (!((width == 64 && height == 4) || (width == 4 && height == 64))) {
      num_mip_modes = NUM_MIP_MODES_FULL(width, height);
    }
  }

  // Find best intra mode for 2Nx2N.
  uvg_pixel *ref_pixels = &lcu->ref.y[lcu_px.x + lcu_px.y * LCU_WIDTH];

  int8_t number_of_modes[MAX_REF_LINE_IDX] = { 0 };
  bool skip_rough_search = (depth == 0 || state->encoder_control->cfg.rdo >= 4);
  if (!skip_rough_search) {
    number_of_modes[0] = search_intra_rough(state,
                                         ref_pixels, LCU_WIDTH,
                                         &refs,
                                         log2_width, candidate_modes,
                                         modes[0], costs[0]);
    // Copy rough results for other reference lines
    for (int line = 1; line < MAX_REF_LINE_IDX; ++line) {
      number_of_modes[line] = number_of_modes[0];
      for (int i = 0; i < number_of_modes[line]; ++i) {
        modes[line][i] = modes[0][i];
        costs[line][i] = costs[0][i];
      }
    }
  } else {
    for(int line = 0; line < MAX_REF_LINE_IDX; ++line) {
      number_of_modes[line] = 67;
      for (int i = 0; i < number_of_modes[line]; ++i) {
        modes[line][i] = i;
        costs[line][i] = MAX_INT;
      }
    }
  }

  uint8_t lines = 1;
  // Find modes with multiple reference lines if in use. Do not use if CU in first row.
  if (state->encoder_control->cfg.mrl && (y_px % LCU_WIDTH) != 0) {
    lines = MAX_REF_LINE_IDX;
  }

  // Set transform depth to current depth, meaning no transform splits.
  uvg_lcu_fill_trdepth(lcu, x_px, y_px, depth, depth);
  // Refine results with slower search or get some results if rough search was skipped.
  const int32_t rdo_level = state->encoder_control->cfg.rdo;
  if (rdo_level >= 2 || skip_rough_search) {
    int number_of_modes_to_search;
    if (rdo_level == 4) {
      number_of_modes_to_search = 67;
    } else if (rdo_level == 2 || rdo_level == 3) {
      number_of_modes_to_search = (cu_width == 4) ? 3 : 2;
    } else {
      // Check only the predicted modes.
      number_of_modes_to_search = 0;
    }
    
    for(int8_t line = 0; line < lines; ++line) {
      // For extra reference lines, only check predicted modes & no MIP search.
      if (line != 0) {
        number_of_modes_to_search = 0;
        num_mip_modes = 0;
      }
      int num_modes_to_check = MIN(number_of_modes[line], number_of_modes_to_search);
      uvg_sort_modes(modes[line], costs[line], number_of_modes[line]);
      // TODO: if rough search is implemented for MIP, sort mip_modes here.
      number_of_modes[line] = search_intra_rdo(state,
                            x_px, y_px, depth,
                            ref_pixels, LCU_WIDTH,
                            candidate_modes,
                            num_modes_to_check,
                            modes[line], trafo[line], costs[line],
                            num_mip_modes,
                            mip_modes, mip_trafo, mip_costs,
                            lcu, line);
    }
  }
  
  uint8_t best_line = 0;
  double best_line_mode_cost = costs[0][0];
  uint8_t best_mip_mode_idx = 0;
  uint8_t best_mode_indices[MAX_REF_LINE_IDX];

  int8_t tmp_best_mode;
  int8_t tmp_best_trafo;
  double tmp_best_cost;
  bool tmp_mip_flag = false;
  bool tmp_mip_transp = false;

  for (int line = 0; line < lines; ++line) {
    best_mode_indices[line] = select_best_mode_index(modes[line], costs[line], number_of_modes[line]);
    if (best_line_mode_cost > costs[line][best_mode_indices[line]]) {
      best_line_mode_cost = costs[line][best_mode_indices[line]];
      best_line = line;
    }
  }

  tmp_best_mode = modes[best_line][best_mode_indices[best_line]];
  tmp_best_trafo = trafo[best_line][best_mode_indices[best_line]];
  tmp_best_cost = costs[best_line][best_mode_indices[best_line]];

  if (num_mip_modes) {
    best_mip_mode_idx = select_best_mode_index(mip_modes, mip_costs, num_mip_modes);
    if (tmp_best_cost > mip_costs[best_mip_mode_idx]) {
      tmp_best_mode = mip_modes[best_mip_mode_idx];
      tmp_best_trafo = mip_trafo[best_mip_mode_idx];
      tmp_best_cost = mip_costs[best_mip_mode_idx];
      tmp_mip_flag = true;
      tmp_mip_transp = (tmp_best_mode >= (num_mip_modes >> 1)) ? 1 : 0;
    }
  }

  if (tmp_mip_flag) {
    // Transform best mode index to proper form.
    // Max mode index is half of max number of modes - 1 (i. e. for size id 2, max mode id is 5)
    tmp_best_mode = (tmp_mip_transp ? tmp_best_mode - (num_mip_modes >> 1) : tmp_best_mode);
  }

  *mode_out =  tmp_best_mode;
  *trafo_out = tmp_best_trafo;
  *cost_out =  tmp_best_cost;
  *mip_flag_out = tmp_mip_flag;
  *mip_transposed_out = tmp_mip_transp;
  *multi_ref_idx_out = tmp_mip_flag ? 0 : best_line;
}
