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
#include <math.h>


#include "cabac.h"
#include "encoder.h"
#include "encoderstate.h"
#include "encode_coding_tree.h"
#include "image.h"
#include "intra.h"
#include "uvg266.h"
#include "rdo.h"
#include "search.h"
#include "transform.h"
#include "strategies/strategies-picture.h"
#include "videoframe.h"
#include "strategies/strategies-quant.h"
#include "uvg_math.h"


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
// static double get_cost(encoder_state_t * const state, 
//                        uvg_pixel *pred, uvg_pixel *orig_block,
//                        cost_pixel_nxn_func *satd_func,
//                        cost_pixel_nxn_func *sad_func,
//                        int width)
// {
//   double satd_cost = satd_func(pred, orig_block);
//   if (TRSKIP_RATIO != 0 && width <= (1 << state->encoder_control->cfg.trskip_max_size) && state->encoder_control->cfg.trskip_enable) {
//     // If the mode looks better with SAD than SATD it might be a good
//     // candidate for transform skip. How much better SAD has to be is
//     // controlled by TRSKIP_RATIO.
//
//     // Add the offset bit costs of signaling 'luma and chroma use trskip',
//     // versus signaling 'luma and chroma don't use trskip' to the SAD cost.
//     const cabac_ctx_t *ctx = &state->search_cabac.ctx.transform_skip_model_luma;
//     double trskip_bits = CTX_ENTROPY_FBITS(ctx, 1) - CTX_ENTROPY_FBITS(ctx, 0);
//
//     
//     // ToDo: Check cost
//     if (state->encoder_control->chroma_format != KVZ_CSP_400) {
//       ctx = &state->search_cabac.ctx.transform_skip_model_chroma;
//       trskip_bits += 2.0 * (CTX_ENTROPY_FBITS(ctx, 1) - CTX_ENTROPY_FBITS(ctx, 0));
//     }
//     
//
//     double sad_cost = TRSKIP_RATIO * sad_func(pred, orig_block) + state->lambda_sqrt * trskip_bits;
//     if (sad_cost < satd_cost) {
//       return sad_cost;
//     }
//   }
//   return satd_cost;
// }


/**
 * \brief Calculate quality of the reconstruction.
 *
 * \param a  bc
 *
 * \return  
 */
static void get_cost_dual(
  encoder_state_t * const state,
  const pred_buffer preds,
  const uvg_pixel *orig_block,
  cost_pixel_nxn_multi_func *satd_twin_func,
  cost_pixel_nxn_multi_func *sad_twin_func,
  int width,
  int height,
  double *costs_out)
{
  #define PARALLEL_BLKS 2
  unsigned satd_costs[PARALLEL_BLKS] = { 0 };
  if (satd_twin_func != NULL) {
    satd_twin_func(preds, orig_block, PARALLEL_BLKS, satd_costs);
  } else {
    satd_costs[0] = uvg_satd_any_size_vtm(width, height, orig_block, width, preds[0], width);
    satd_costs[1] = uvg_satd_any_size_vtm(width, height, orig_block, width, preds[1], width);
  }
  unsigned unsigned_sad_costs[PARALLEL_BLKS] = { 0 };
  if (sad_twin_func != NULL) {
    sad_twin_func(preds, orig_block, PARALLEL_BLKS, unsigned_sad_costs);
  } else {
    unsigned_sad_costs[0] = uvg_reg_sad(preds[0], orig_block, width, height, width, width);
    unsigned_sad_costs[1] = uvg_reg_sad(preds[1], orig_block, width, height, width, width);
  }
  costs_out[0] = (double)MIN(satd_costs[0], unsigned_sad_costs[0] * 2);
  costs_out[1] = (double)MIN(satd_costs[1], unsigned_sad_costs[1] * 2);

  // TODO: width and height
  //if (TRSKIP_RATIO != 0 && width <= (1 << state->encoder_control->cfg.trskip_max_size) && state->encoder_control->cfg.trskip_enable) {
  //  // If the mode looks better with SAD than SATD it might be a good
  //  // candidate for transform skip. How much better SAD has to be is
  //  // controlled by TRSKIP_RATIO.

  //  // Add the offset bit costs of signaling 'luma and chroma use trskip',
  //  // versus signaling 'luma and chroma don't use trskip' to the SAD cost.
  //  const cabac_ctx_t *ctx = &state->cabac.ctx.transform_skip_model_luma;
  //  double trskip_bits = CTX_ENTROPY_FBITS(ctx, 1) - CTX_ENTROPY_FBITS(ctx, 0);

  //  
  //  // ToDo: check cost
  //  if (state->encoder_control->chroma_format != UVG_CSP_400) {
  //    ctx = &state->cabac.ctx.transform_skip_model_chroma;
  //    trskip_bits += 2.0 * (CTX_ENTROPY_FBITS(ctx, 1) - CTX_ENTROPY_FBITS(ctx, 0));
  //  }
  //  

  //  unsigned unsigned_sad_costs[PARALLEL_BLKS] = { 0 };
  //  double sad_costs[PARALLEL_BLKS] = { 0 };
  //  sad_twin_func(preds, orig_block, PARALLEL_BLKS, unsigned_sad_costs);
  //  for (int i = 0; i < PARALLEL_BLKS; ++i) {
  //    sad_costs[i] = TRSKIP_RATIO * (double)unsigned_sad_costs[i] + state->lambda_sqrt * trskip_bits;
  //    if (sad_costs[i] < (double)satd_costs[i]) {
  //      costs_out[i] = sad_costs[i];
  //    }
  //  }
  //}

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
                                   lcu_t *const lcu, const int width, const int height,
                                   const vector2d_t lcu_px)
{
  int8_t scan_idx = SCAN_DIAG;
  int32_t i;
  // ToDo: large block support in VVC?
  uint32_t sig_coeffgroup_flag[32 * 32] = { 0 };

  const uint32_t log2_block_width =  uvg_g_convert_to_log2[width];
  const uint32_t log2_block_height = uvg_g_convert_to_log2[height];
  const uint32_t log2_cg_size = uvg_g_log2_sbb_size[log2_block_width][log2_block_height][0]
    + uvg_g_log2_sbb_size[log2_block_width][log2_block_height][1];
  const uint32_t * const scan = uvg_get_scan_order_table(SCAN_GROUP_4X4, scan_idx, log2_block_width, log2_block_height);
  const uint32_t * const scan_cg = uvg_get_scan_order_table(SCAN_GROUP_UNGROUPED, scan_idx, log2_block_width, log2_block_height);

  coeff_t coeff_y[TR_MAX_WIDTH * TR_MAX_WIDTH];
  uvg_get_sub_coeff(coeff_y, lcu->coeff.y, lcu_px.x, lcu_px.y, width, height, LCU_WIDTH);

  signed scan_cg_last = -1;
  signed scan_pos_last = -1;

  for (int i = 0; i < width * height; i++) {
    if (coeff_y[scan[i]]) {
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
  const cu_loc_t* const cu_loc,
  double cost_treshold,
  intra_search_data_t *const search_data,
  lcu_t *const lcu,
  enum uvg_tree_type tree_type)
{
  const uint8_t width = cu_loc->width;
  const uint8_t height = cu_loc->height; // TODO: height for non-square blocks
  const uint8_t width_c = cu_loc->chroma_width;
  const uint8_t height_c = cu_loc->chroma_height;
  
  const vector2d_t lcu_px = { cu_loc->local_x, cu_loc->local_y };

  const bool reconstruct_chroma = false;// (depth != 4 || (depth == 4 && (x_px & 4 && y_px & 4))) && state->encoder_control->chroma_format != UVG_CSP_400;
  cu_info_t* pred_cu = &search_data->pred_cu;

  double split_cost = INT32_MAX;
  double nosplit_cost = INT32_MAX;

  cabac_data_t cabac_data;
  memcpy(&cabac_data, &state->search_cabac, sizeof(cabac_data_t));
  state->search_cabac.update = 1;

  if (width <= TR_MAX_WIDTH && height <= TR_MAX_WIDTH) {

    const bool mts_enabled = (state->encoder_control->cfg.mts == UVG_MTS_INTRA || state->encoder_control->cfg.mts == UVG_MTS_BOTH)
      && PU_IS_TU(pred_cu);

    nosplit_cost = 0.0;
    const bool has_been_split = 1 << pred_cu->log2_width != cu_loc->width ||
                                1 << pred_cu->log2_height != cu_loc->height;

    cbf_clear(&pred_cu->cbf, COLOR_Y);
    if (reconstruct_chroma) {
      cbf_clear(&pred_cu->cbf, COLOR_U);
      cbf_clear(&pred_cu->cbf, COLOR_V);
    }

    const int8_t chroma_mode = reconstruct_chroma ? (!pred_cu->intra.mip_flag ? pred_cu->intra.mode : 0) : -1;
    double best_rd_cost = MAX_INT;
    int best_tr_idx = 0;
    int best_lfnst_idx = 0;

    uint8_t trafo;
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
      // Do not do MTS search if ISP mode is used
      num_transforms = pred_cu->intra.isp_mode == ISP_MODE_NO_ISP ? num_transforms : 1;
    }
    const int mts_start = trafo;
    if (state->encoder_control->cfg.trskip_enable 
      && width <= (1 << state->encoder_control->cfg.trskip_max_size)
      && height <= (1 << state->encoder_control->cfg.trskip_max_size)
      && PU_IS_TU(pred_cu)
      && pred_cu->intra.isp_mode == ISP_MODE_NO_ISP) {
      num_transforms = MAX(num_transforms, 2);
    }
    pred_cu->intra.mode_chroma = -1;
    
    const int max_tb_size = TR_MAX_WIDTH;
    // LFNST search params
    int max_lfnst_idx = width > max_tb_size || height > max_tb_size ? 0 : 2;
    if(pred_cu->intra.mip_flag && (width < 16 || height < 16)) {
      max_lfnst_idx = 0;
    }
    
    int start_idx = 0;
    int end_lfnst_idx = state->encoder_control->cfg.lfnst && PU_IS_TU(pred_cu) &&
                  uvg_can_use_isp_with_lfnst(width, height, pred_cu->intra.isp_mode, tree_type) ? max_lfnst_idx : 0;
    for (int i = start_idx; i < end_lfnst_idx + 1; ++i) {
      search_data->lfnst_costs[i] = MAX_DOUBLE;
    }

    for (trafo = mts_start; trafo < num_transforms; trafo++) {
      for (int lfnst_idx = start_idx; lfnst_idx <= end_lfnst_idx; lfnst_idx++) {
        // Initialize lfnst variables
        search_data->best_isp_cbfs = 0;
        pred_cu->tr_idx = trafo;
        pred_cu->tr_skip = trafo == MTS_SKIP;
        pred_cu->lfnst_idx = lfnst_idx;
        pred_cu->violates_lfnst_constrained_luma = false;
        pred_cu->violates_lfnst_constrained_chroma = false;
        pred_cu->lfnst_last_scan_pos = false;

        bool constraints[2] = {false, false};
        if (mts_enabled) {
          pred_cu->mts_last_scan_pos = 0;
          pred_cu->violates_mts_coeff_constraint = 0;

          if (trafo == MTS_SKIP && ((width > (1 << state->encoder_control->cfg.trskip_max_size)
            || (height > (1 << state->encoder_control->cfg.trskip_max_size)))
            || !PU_IS_TU(pred_cu)
            || !state->encoder_control->cfg.trskip_enable)) {
            continue;
          }
        }
        // MTS and LFNST cannot be on at the same time
        if (pred_cu->lfnst_idx > 0 && pred_cu->tr_idx > 0) {
          continue;
        }

        if (!has_been_split && (lfnst_idx != 0 || trafo != 0)) {
          memcpy(&state->search_cabac, &cabac_data, sizeof(cabac_data));
          state->search_cabac.update = 1;
        }
        double rd_cost;
        if (pred_cu->intra.isp_mode != ISP_MODE_NO_ISP) {
          rd_cost = uvg_recon_and_estimate_cost_isp(
            state,
            cu_loc,
            cost_treshold,
            search_data,
            lcu,
            &constraints[0]
          );
          constraints[1] = search_data->best_isp_cbfs != 0;
        }
        else {
          uvg_intra_recon_cu(
            state,
            search_data,
            cu_loc,
            pred_cu,
            lcu,
            UVG_LUMA_T,
            true,
            false
          );
        }
        if (pred_cu->intra.isp_mode != ISP_MODE_NO_ISP && search_data->best_isp_cbfs == 0) continue;

        if ((trafo != 0 || lfnst_idx != 0) && !cbf_is_set(pred_cu->cbf, COLOR_Y)) continue;
        
        derive_mts_constraints(pred_cu, lcu, width, height, lcu_px);
        if (pred_cu->tr_idx > 1) {
          if (pred_cu->violates_mts_coeff_constraint || !pred_cu->
              mts_last_scan_pos) {
            continue;
          }
        }
        
        if (trafo != MTS_SKIP && end_lfnst_idx != 0 && pred_cu->intra.isp_mode == ISP_MODE_NO_ISP) {
          uvg_derive_lfnst_constraints(
            pred_cu,
            constraints,
            lcu->coeff.y,
            width,
            height,
            &lcu_px,
            COLOR_Y);
        }

        if (!constraints[1] && cbf_is_set(pred_cu->cbf, COLOR_Y)) {
          //end_idx = 0;
          if (pred_cu->lfnst_idx > 0) {
            continue;
          }
        }


        if (pred_cu->intra.isp_mode == ISP_MODE_NO_ISP) {
          rd_cost = uvg_cu_rd_cost_luma(
            state,
            cu_loc,
            pred_cu,
            lcu,
            search_data->best_isp_cbfs);
        }
        double transform_bits = 0;
        if (state->encoder_control->cfg.lfnst && PU_IS_TU(pred_cu) &&
          trafo != MTS_SKIP && end_lfnst_idx != 0 && (cbf_is_set(pred_cu->cbf, COLOR_Y) || search_data->best_isp_cbfs != 0)) {
          if ((!constraints[0] && (constraints[1] || pred_cu->intra.isp_mode != ISP_MODE_NO_ISP))) {
            transform_bits += CTX_ENTROPY_FBITS(
              &state->search_cabac.ctx.lfnst_idx_model[tree_type == UVG_LUMA_T],
              lfnst_idx != 0);
            if (lfnst_idx > 0) {
              transform_bits += CTX_ENTROPY_FBITS(
                &state->search_cabac.ctx.lfnst_idx_model[2],
                lfnst_idx == 2);
            }
          }
        }
        if (num_transforms > 2 && trafo != MTS_SKIP
            && (cbf_is_set(pred_cu->cbf, COLOR_Y) || search_data->best_isp_cbfs != 0)
            && pred_cu->intra.isp_mode == ISP_MODE_NO_ISP
            && lfnst_idx == 0
            && width <= 32
            && height <= 32
            && !pred_cu->violates_mts_coeff_constraint && pred_cu->
            mts_last_scan_pos) {

          bool symbol = trafo != 0;
          int ctx_idx = 0;
          transform_bits += CTX_ENTROPY_FBITS(
            &state->search_cabac.ctx.mts_idx_model[ctx_idx],
            symbol);

          ctx_idx++;
          for (int i = 0; i < 3 && symbol; i++, ctx_idx++) {
            symbol = trafo > i + MTS_DST7_DST7 ? 1 : 0;
            transform_bits += CTX_ENTROPY_FBITS(
              &state->search_cabac.ctx.mts_idx_model[ctx_idx],
              symbol);
          }

        }
        rd_cost += transform_bits * state->lambda;

        search_data->lfnst_costs[lfnst_idx] = MIN(
          search_data->lfnst_costs[lfnst_idx],
          rd_cost);
        if (rd_cost < best_rd_cost) {
          best_rd_cost = rd_cost;
          best_lfnst_idx = pred_cu->lfnst_idx;
          best_tr_idx = pred_cu->tr_idx;
          if (best_tr_idx == MTS_SKIP) break;
          // Very unlikely that further search is necessary if skip seems best option
        }
      } // end mts index loop (tr_idx)
      if (reconstruct_chroma) {
        int8_t luma_mode = pred_cu->intra.mode;
        pred_cu->intra.mode_chroma = chroma_mode;
        // TODO: Maybe check the jccr mode here also but holy shit is the interface of search_intra_rdo bad currently
        uvg_intra_recon_cu(
          state,
          search_data,
          cu_loc,
          pred_cu,
          lcu,
          UVG_BOTH_T,
          false,
          true
          );
        best_rd_cost += uvg_cu_rd_cost_chroma(
          state,
          pred_cu,
          lcu,
          cu_loc);
        pred_cu->intra.mode = luma_mode;

        // Check lfnst constraints for chroma
        if (pred_cu->lfnst_idx > 0) {
          // Temp constraints. Updating the actual pred_cu constraints here will break things later
          bool constraints[2] = {pred_cu->violates_lfnst_constrained_chroma,
                                 pred_cu->lfnst_last_scan_pos};
          uvg_derive_lfnst_constraints(
            pred_cu,
            constraints,
            lcu->coeff.u,
            width_c,
            height_c,
            &lcu_px,
            COLOR_U);
          if (constraints[0] || !constraints[1]) {
            best_lfnst_idx = 0;
            continue;
          }
          uvg_derive_lfnst_constraints(
            pred_cu,
            constraints,
            lcu->coeff.u,
            width_c,
            height_c,
            &lcu_px,
            COLOR_U);
          if (constraints[0] || !constraints[1]) {
            best_lfnst_idx = 0;
            continue;
          }
        }
      }
      if (best_tr_idx == MTS_SKIP) break;
    }
    if(reconstruct_chroma) {
      int8_t luma_mode = pred_cu->intra.mode;
      pred_cu->intra.mode_chroma = chroma_mode;
      uvg_intra_recon_cu(state,
                         search_data, cu_loc,
                         pred_cu, lcu,
                         UVG_BOTH_T,
                         false,
                         true);
      best_rd_cost += uvg_cu_rd_cost_chroma(state, pred_cu, lcu, cu_loc);
      pred_cu->intra.mode = luma_mode;
    }
    pred_cu->tr_skip = best_tr_idx == MTS_SKIP;
    pred_cu->tr_idx = best_tr_idx;
    pred_cu->lfnst_idx = best_lfnst_idx;
    pred_cu->lfnst_last_scan_pos = false;
    pred_cu->violates_lfnst_constrained_luma = false;
    nosplit_cost += best_rd_cost;

    // Early stop condition for the recursive search.
    // If the cost of any 1/4th of the transform is already larger than the
    // whole transform, assume that splitting further is a bad idea.
    if (nosplit_cost <= cost_treshold) {
      memcpy(&state->search_cabac, &cabac_data, sizeof(cabac_data));
      return nosplit_cost;
    }
  }
    
  
  // Recurse further if all of the following:
  // - Current depth is less than maximum depth of the search (max_depth).
  //   - Maximum transform hierarchy depth is constrained by clipping
  //     max_depth.
  // - Min transform size hasn't been reached (MAX_PU_DEPTH).
  else {
    split_cost = 0;


    enum split_type split;
    if (cu_loc->width > TR_MAX_WIDTH && cu_loc->height > TR_MAX_WIDTH) {
      split = QT_SPLIT;
    }
    else if (cu_loc->width > TR_MAX_WIDTH) {
      split = BT_VER_SPLIT;
    }
    else {
      split = BT_HOR_SPLIT;
    }

    cu_loc_t split_cu_loc[4];
    const int split_count = uvg_get_split_locs(cu_loc, split, split_cu_loc,NULL);
    for (int i = 0; i < split_count; ++i) {
      split_cost += search_intra_trdepth(state, &split_cu_loc[i], nosplit_cost, search_data, lcu, tree_type);
    }
  }
  memcpy(&state->search_cabac, &cabac_data, sizeof(cabac_data));

  if (!PU_IS_TU(pred_cu) || split_cost < nosplit_cost) {
    return split_cost;
  } else {
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

static int search_intra_chroma_rough(
  encoder_state_t * const state,
  intra_search_data_t* chroma_data,
  lcu_t* lcu,
  int8_t luma_mode,
  enum uvg_tree_type tree_type,
  const cu_loc_t* const cu_loc)
{
  const int_fast8_t log2_width_c = uvg_g_convert_to_log2[cu_loc->chroma_width];
  const vector2d_t pic_px = { state->tile->frame->width, state->tile->frame->height };
  const vector2d_t luma_px = { cu_loc->x, cu_loc->y};
  const int width = 1 << log2_width_c;
  const int height = width; // TODO: height for non-square blocks

  const cu_loc_t loc = { luma_px.x, luma_px.y, width, height, width, height };

  uvg_intra_references refs_u;
  uvg_intra_build_reference(state, &loc, &loc, COLOR_U, &luma_px, &pic_px, lcu, &refs_u, state->encoder_control->cfg.wpp, NULL, 0, 0);

  uvg_intra_references refs_v;
  uvg_intra_build_reference(state, &loc, &loc, COLOR_V, &luma_px, &pic_px, lcu, &refs_v, state->encoder_control->cfg.wpp, NULL, 0, 0);

  vector2d_t lcu_cpx = { (cu_loc->local_x) / 2, (cu_loc->local_y) / 2 };
  uvg_pixel* orig_u = &lcu->ref.u[lcu_cpx.x + lcu_cpx.y * LCU_WIDTH_C];
  uvg_pixel* orig_v = &lcu->ref.v[lcu_cpx.x + lcu_cpx.y * LCU_WIDTH_C];
  
  //cost_pixel_nxn_func *const sad_func = uvg_pixels_get_sad_func(width);
    
  uvg_pixel _pred[32 * 32 + SIMD_ALIGNMENT];
  uvg_pixel *pred = ALIGNED_POINTER(_pred, SIMD_ALIGNMENT);

  uvg_pixel _orig_block[32 * 32 + SIMD_ALIGNMENT];
  uvg_pixel *orig_block = ALIGNED_POINTER(_orig_block, SIMD_ALIGNMENT);

  uvg_pixels_blit(orig_u, orig_block, width, height, LCU_WIDTH_C, width);
  int modes_count = (state->encoder_control->cfg.cclm ? 8 : 5);
  for (int i = 0; i < modes_count; ++i) {
    const int8_t mode_chroma = chroma_data[i].pred_cu.intra.mode_chroma;
    if (mode_chroma == luma_mode || mode_chroma == 0 || mode_chroma >= 81) continue;
    uvg_intra_predict(state, &refs_u, cu_loc, &loc, COLOR_U, pred, &chroma_data[i], lcu);
    //costs[i] += get_cost(encoder_state, pred, orig_block, satd_func, sad_func, width);
    switch (width) {
      case 4: chroma_data[i].cost += uvg_satd_4x4(pred, orig_block);
        break;
      case 8: chroma_data[i].cost += uvg_satd_8x8(pred, orig_block);
        break;
      case 16: chroma_data[i].cost += uvg_satd_16x16(pred, orig_block);
        break;
      case 32: chroma_data[i].cost += uvg_satd_32x32(pred, orig_block);
        break;
      default: assert(0);
    }
  }

  uvg_pixels_blit(orig_v, orig_block, width, height, LCU_WIDTH_C, width);
  for (int i = 0; i < modes_count; ++i) {
    const int8_t mode_chroma = chroma_data[i].pred_cu.intra.mode_chroma;
    if (mode_chroma == luma_mode || mode_chroma == 0 || mode_chroma >= 81) continue;
    uvg_intra_predict(state, &refs_v, cu_loc, &loc, COLOR_V, pred, &chroma_data[i], lcu);
    //costs[i] += get_cost(encoder_state, pred, orig_block, satd_func, sad_func, width);
    switch (width) {
      case 4: chroma_data[i].cost += uvg_satd_4x4(pred, orig_block);
        break;
      case 8: chroma_data[i].cost += uvg_satd_8x8(pred, orig_block);
        break;
      case 16: chroma_data[i].cost += uvg_satd_16x16(pred, orig_block);
        break;
      case 32: chroma_data[i].cost += uvg_satd_32x32(pred, orig_block);
        break;
      default: assert(0);
    }
  }
  sort_modes(chroma_data, modes_count);
  if (modes_count > 5 && chroma_data[7].pred_cu.intra.mode_chroma > 81) modes_count--;
  if (modes_count > 5 && chroma_data[6].pred_cu.intra.mode_chroma > 81) modes_count--;
  return modes_count;
  // for (int i = 0; i < modes_count; ++i) {
  //   const double bits = kvz_chroma_mode_bits(state, chroma_data[i].pred_cu.intra.mode_chroma, chroma_data[i].pred_cu.intra.mode);
  //   chroma_data[i].bits = bits;
  //   chroma_data[i].cost = bits * state->lambda_sqrt;
  // }
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
/*
static int16_t search_intra_rough(
  encoder_state_t * const state,
  uvg_pixel *orig,
  int32_t origstride,
  uvg_intra_references *refs,
  int log2_width,
  int8_t *intra_preds,
  intra_search_data_t* modes_out,
  cu_info_t* const pred_cu,
  uint8_t mip_ctx)
{
  #define PARALLEL_BLKS 2 // TODO: use 4 for AVX-512 in the future?
  assert(log2_width >= 2 && log2_width <= 5);
  int_fast8_t width = 1 << log2_width;
  cost_pixel_nxn_func *satd_func = uvg_pixels_get_satd_func(width);
  cost_pixel_nxn_func *sad_func = uvg_pixels_get_sad_func(width);
  cost_pixel_nxn_multi_func *satd_dual_func = uvg_pixels_get_satd_dual_func(width);
  cost_pixel_nxn_multi_func *sad_dual_func = uvg_pixels_get_sad_dual_func(width);
  int8_t modes[UVG_NUM_INTRA_MODES];
  double costs[UVG_NUM_INTRA_MODES];

  // const uvg_config *cfg = &state->encoder_control->cfg;
  // const bool filter_boundary = !(cfg->lossless && cfg->implicit_rdpcm);

  // Temporary block arrays
  uvg_pixel _preds[PARALLEL_BLKS * 32 * 32 + SIMD_ALIGNMENT];
  pred_buffer preds = ALIGNED_POINTER(_preds, SIMD_ALIGNMENT);
  
  uvg_pixel _orig_block[32 * 32 + SIMD_ALIGNMENT];
  uvg_pixel *orig_block = ALIGNED_POINTER(_orig_block, SIMD_ALIGNMENT);

  // Store original block for SAD computation
  uvg_pixels_blit(orig, orig_block, width, height, origstride, width);

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
  cu_loc_t loc = { 0, 0, width, height, width, height };
  intra_search_data_t search_proxy;
  FILL(search_proxy, 0);
  search_proxy.pred_cu = *pred_cu;

  for (int mode = 2; mode <= 66; mode += PARALLEL_BLKS * offset) {
    
    double costs_out[PARALLEL_BLKS] = { 0 };
    for (int i = 0; i < PARALLEL_BLKS; ++i) {
      if (mode + i * offset <= 66) {
        search_proxy.pred_cu.intra.mode = mode + i*offset;
        uvg_intra_predict(state, refs, &loc, COLOR_Y, preds[i], &search_proxy, NULL);
      }
    }
    
    //TODO: add generic version of get cost  multi
    get_cost_dual(state, preds, orig_block, satd_dual_func, sad_dual_func, width, costs_out);

    for (int i = 0; i < PARALLEL_BLKS; ++i) {
      if (mode + i * offset <= 66) {
        costs[modes_selected] = costs_out[i];
        modes[modes_selected] = mode + i * offset;
        min_cost = (int32_t)MIN(min_cost, costs[modes_selected]);
        max_cost = (int32_t)MAX(max_cost, costs[modes_selected]);
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
            uvg_intra_predict(state, refs, &loc, COLOR_Y, preds[i], &search_proxy, NULL);
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
      uvg_intra_predict(state, refs, &loc, COLOR_Y, preds[0], &search_proxy, NULL);
      costs[modes_selected] = get_cost(state, preds[0], orig_block, satd_func, sad_func, width);
      modes[modes_selected] = mode;
      ++modes_selected;
    }
  }

  // Add prediction mode coding cost as the last thing. We don't want this
  // affecting the halving search.
  const double not_mrl = state->encoder_control->cfg.mrl ? CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.multi_ref_line[0]), 0) : 0;
  const double not_mip = state->encoder_control->cfg.mip ? CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.mip_flag[mip_ctx]), 0) : 0;
  const double mpm_mode_bit = CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.intra_luma_mpm_flag_model), 1);
  const double not_mpm_mode_bit = CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.intra_luma_mpm_flag_model), 0);
  const double planar_mode_flag = CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.luma_planar_model[1]), 0);
  const double not_planar_mode_flag = CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.luma_planar_model[1]), 1);
  for (int mode_i = 0; mode_i < modes_selected; ++mode_i) {
    int i = 0;
    int smaller_than_pred = 0;
    double bits;
    for (; i < INTRA_MPM_COUNT; i++) {
      if (intra_preds[i] == modes[mode_i]) {
        break;
      }
      if (modes[mode_i] < intra_preds[i]) {
        smaller_than_pred += 1;
      }
    }
    if (i == 0) {
      bits = planar_mode_flag + mpm_mode_bit;
    }
    else if (i < INTRA_MPM_COUNT) {
      bits = not_planar_mode_flag + mpm_mode_bit + MIN(i, 4);
    }
    else {
      bits = not_mpm_mode_bit + 5 + (modes[mode_i] - smaller_than_pred > 3);
    }
    bits += not_mrl + not_mip;
    costs[mode_i] += state->lambda_sqrt * bits;
    modes_out[mode_i].cost = costs[mode_i];
    modes_out[mode_i].pred_cu = *pred_cu;
    modes_out[mode_i].pred_cu.intra.mode = modes[mode_i];
    modes_out[mode_i].pred_cu.intra.mode_chroma = modes[mode_i];
  }

  #undef PARALLEL_BLKS
  return modes_selected;
}*/


static INLINE double count_bits(
  encoder_state_t* const state,
  int8_t* intra_preds,
  const double not_mrl,
  const double not_mip,
  const double mpm_mode_bit,
  const double not_mpm_mode_bit,
  const double planar_mode_flag,
  const double not_planar_mode_flag,
  const double not_isp_flag,
  int8_t mode
)
{
  int i = 0;
  int smaller_than_pred = 0;
  double bits;
  for (; i < INTRA_MPM_COUNT; i++) {
    if (intra_preds[i] == mode) {
      break;
    }
    if (mode > intra_preds[i]) {
      smaller_than_pred += 1;
    }
  }
  if (i == 0) {
    bits = planar_mode_flag + mpm_mode_bit;
  }
  else if (i < INTRA_MPM_COUNT) {
    bits = not_planar_mode_flag + mpm_mode_bit + MIN(i, 4);
  }
  else {
    bits = not_mpm_mode_bit + 5 + (mode - smaller_than_pred > 2);
  }
  bits += not_mrl + not_mip + not_isp_flag;
  return bits;
}

static uint8_t search_intra_rough(
  encoder_state_t * const state,
  const cu_loc_t* const cu_loc,
  uvg_pixel *orig,
  int32_t origstride,
  uvg_intra_references *refs,
  int width,
  int height,
  int8_t *intra_preds,
  intra_search_data_t* modes_out,
  cu_info_t* const pred_cu,
  uint8_t mip_ctx)
{
  #define PARALLEL_BLKS 2 // TODO: use 4 for AVX-512 in the future?
  assert(width >= 4 && width <= 32);
  // cost_pixel_nxn_func *satd_func = kvz_pixels_get_satd_func(width);
  // cost_pixel_nxn_func *sad_func = kvz_pixels_get_sad_func(width);
  cost_pixel_nxn_multi_func *satd_dual_func = uvg_pixels_get_satd_dual_func(width, height);
  cost_pixel_nxn_multi_func *sad_dual_func = uvg_pixels_get_sad_dual_func(width, height);
  bool mode_checked[UVG_NUM_INTRA_MODES] = {0};
  double costs[UVG_NUM_INTRA_MODES];

  // const kvz_config *cfg = &state->encoder_control->cfg;
  // const bool filter_boundary = !(cfg->lossless && cfg->implicit_rdpcm);

  // Temporary block arrays
  uvg_pixel _preds[PARALLEL_BLKS * 32 * 32 + SIMD_ALIGNMENT];
  pred_buffer preds = ALIGNED_POINTER(_preds, SIMD_ALIGNMENT);
  
  uvg_pixel _orig_block[32 * 32 + SIMD_ALIGNMENT];
  uvg_pixel *orig_block = ALIGNED_POINTER(_orig_block, SIMD_ALIGNMENT);

  // Store original block for SAD computation
  uvg_pixels_blit(orig, orig_block, width, height, origstride, width);

  int8_t modes_selected = 0;
  // Note: get_cost and get_cost_dual may return negative costs.
  double min_cost;
  double max_cost;

  struct mode_cost {
    int8_t mode;
    double cost;
  };
  
  const double not_mrl = state->encoder_control->cfg.mrl && (cu_loc->y % LCU_WIDTH) ? CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.multi_ref_line[0]), 0) : 0;
  const double not_mip = state->encoder_control->cfg.mip ? CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.mip_flag[mip_ctx]), 0) : 0;
  const double mpm_mode_bit = CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.intra_luma_mpm_flag_model), 1);
  const double not_mpm_mode_bit = CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.intra_luma_mpm_flag_model), 0);
  const double planar_mode_flag = CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.luma_planar_model[1]), 0);
  const double not_planar_mode_flag = CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.luma_planar_model[1]), 1);
  const double not_isp_flag = state->encoder_control->cfg.isp && uvg_can_use_isp(width, height) ? CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.intra_subpart_model[0]), 0) : 0;

  const uint8_t mode_list_size = state->encoder_control->cfg.mip ? 6 : 3;
  struct mode_cost best_six_modes[6];
  // Initial offset decides how many modes are tried before moving on to the
  // recursive search.

  // Calculate SAD for evenly spaced modes to select the starting point for 
  // the recursive search.
  intra_search_data_t search_proxy;
  FILL(search_proxy, 0);
  search_proxy.pred_cu = *pred_cu;

  int offset = 1 << state->encoder_control->cfg.intra_rough_search_levels;
  search_proxy.pred_cu.intra.mode = 0;
  uvg_intra_predict(state, refs, cu_loc, cu_loc, COLOR_Y, preds[0], &search_proxy, NULL);
  search_proxy.pred_cu.intra.mode = 1;
  uvg_intra_predict(state, refs, cu_loc, cu_loc, COLOR_Y, preds[1], &search_proxy, NULL);
  get_cost_dual(state, preds, orig_block, satd_dual_func, sad_dual_func, width, height, costs);
  mode_checked[0] = true;
  mode_checked[1] = true;
  costs[0] += count_bits(
    state,
    intra_preds,
    not_mrl,
    not_mip,
    mpm_mode_bit,
    not_mpm_mode_bit,
    planar_mode_flag,
    not_planar_mode_flag,
    not_isp_flag, 0) * state->lambda_sqrt;
  costs[1] += count_bits(
    state,
    intra_preds,
    not_mrl,
    not_mip,
    mpm_mode_bit,
    not_mpm_mode_bit,
    planar_mode_flag,
    not_planar_mode_flag,
    not_isp_flag, 1) * state->lambda_sqrt;
  if(costs[0] < costs[1]) {
    min_cost = costs[0];
    max_cost = costs[1];
    best_six_modes[0].mode = 0;
    best_six_modes[0].cost = costs[0];
    best_six_modes[1].mode = 1;
    best_six_modes[1].cost = costs[1];
  }
  else {
    min_cost = costs[1];
    max_cost = costs[0];
    best_six_modes[1].mode = 0;
    best_six_modes[1].cost = costs[0];
    best_six_modes[0].mode = 1;
    best_six_modes[0].cost = costs[1];    
  }
  best_six_modes[2].cost = MAX_DOUBLE;
  best_six_modes[3].cost = MAX_DOUBLE;
  best_six_modes[4].cost = MAX_DOUBLE;
  best_six_modes[5].cost = MAX_DOUBLE;
  for (int mode = 2 + offset / 2; mode <= 66; mode += PARALLEL_BLKS * offset) {
    
    double costs_out[PARALLEL_BLKS] = { 0 };
    for (int i = 0; i < PARALLEL_BLKS; ++i) {
      if (mode + i * offset <= 66) {
        search_proxy.pred_cu.intra.mode = mode + i*offset;
        uvg_intra_predict(state, refs, cu_loc, cu_loc, COLOR_Y, preds[i], &search_proxy, NULL);
      }
    }
    
    //TODO: add generic version of get cost  multi
    get_cost_dual(state, preds, orig_block, satd_dual_func, sad_dual_func, width, height, costs_out);
    for (int i = 0; i < PARALLEL_BLKS; ++i) {
      if (mode + i * offset <= 66) {
        costs_out[i] += count_bits(
          state,
          intra_preds,
          not_mrl,
          not_mip,
          mpm_mode_bit,
          not_mpm_mode_bit,
          planar_mode_flag,
          not_planar_mode_flag,
          not_isp_flag, mode + i * offset) * state->lambda_sqrt;
      }
    }

    for (int i = 0; i < PARALLEL_BLKS; ++i) {
      int8_t mode_i = mode + i* offset;
      if (mode_i <= 66) {
        costs[mode_i] = costs_out[i];
        mode_checked[mode_i] = true;
        min_cost = MIN(min_cost, costs[mode_i]);
        max_cost = MAX(max_cost, costs[mode_i]);
        ++modes_selected;
        for (int j = 0; j < mode_list_size; j++) {
          if (costs[mode_i] < best_six_modes[j].cost) {
            for(int k = mode_list_size - 1; k > j; k--) {
              best_six_modes[k] = best_six_modes[k - 1];
            }
            best_six_modes[j].cost = costs[mode_i];
            best_six_modes[j].mode = mode_i;
            break;
          }
        }
      }
    }
  }
  offset >>= 1;
  // Skip recursive search if all modes have the same cost.
  if (min_cost != max_cost) {
    // Do a recursive search to find the best mode, always centering on the
    // current best mode.
    for (; offset > 0; offset >>= 1) {

      struct mode_cost temp_best_six_modes[6];
      memcpy(temp_best_six_modes, best_six_modes, sizeof(temp_best_six_modes));
      int8_t modes_to_check[12];
      int num_modes_to_check = 0;
      for(int i = 0; i < mode_list_size; i++) {
        int8_t center_node = best_six_modes[i].mode;
        if(offset != 0 && (center_node < 3 || center_node > 65)) continue;
        int8_t test_modes[] = { center_node - offset, center_node + offset };
        for(int j = 0; j < 2; j++) {
          if((test_modes[j] >= 2 && test_modes[j] <= 66) && mode_checked[test_modes[j]] == false) {
            modes_to_check[num_modes_to_check++] = test_modes[j];
            mode_checked[test_modes[j]] = true;
          }
        }
      }
      while (num_modes_to_check & (PARALLEL_BLKS - 1)) {
        modes_to_check[num_modes_to_check++] = 1;
      } 
      for (int i = 0; i < num_modes_to_check; i += PARALLEL_BLKS) {
        double costs_out[PARALLEL_BLKS] = { 0 };        
      
        for (int block = 0; block < PARALLEL_BLKS; ++block) {
          search_proxy.pred_cu.intra.mode = modes_to_check[block + i];
          uvg_intra_predict(state, refs, cu_loc, cu_loc, COLOR_Y, preds[block], &search_proxy, NULL);
        
        }

        //TODO: add generic version of get cost multi
        get_cost_dual(state, preds, orig_block, satd_dual_func, sad_dual_func, width, height, costs_out);
        for (int block = 0; block < PARALLEL_BLKS; ++block) {
            costs_out[block] += count_bits(
              state,
              intra_preds,
              not_mrl,
              not_mip,
              mpm_mode_bit,
              not_mpm_mode_bit,
              planar_mode_flag,
              not_planar_mode_flag,
              not_isp_flag, modes_to_check[block + i]) * state->lambda_sqrt;
          
        }

        for (int block = 0; block < PARALLEL_BLKS; ++block) {
          int8_t mode = modes_to_check[i + block];
          if (mode == 1) continue;
          costs[mode] = costs_out[block];
          for (int j = 0; j < mode_list_size; j++) {
            if (costs[mode] < best_six_modes[j].cost) {
              for (int k = mode_list_size - 1; k > j; k--) {
                best_six_modes[k] = best_six_modes[k - 1];
              }
              best_six_modes[j].cost = costs[mode];
              best_six_modes[j].mode = mode;
              break;
            }
          }          
          
        }
      }
    }
  }

  // Add prediction mode coding cost as the last thing. We don't want this
  // affecting the halving search.
  for(int i=0; i < mode_list_size; i++) {
    const int8_t mode = best_six_modes[i].mode;
    modes_out[i].cost = costs[mode];
    modes_out[i].pred_cu = *pred_cu;
    modes_out[i].pred_cu.intra.mode = mode;
    modes_out[i].pred_cu.intra.mode_chroma = mode;

  }
  
  #undef PARALLEL_BLKS
  return mode_list_size;
}


static void get_rough_cost_for_2n_modes(
  encoder_state_t* const state,
  uvg_intra_references* refs,
  const cu_loc_t* const cu_loc,
  uvg_pixel *orig,
  int orig_stride,
  intra_search_data_t *search_data,
  int num_modes,
  uint8_t mip_ctx)
{
#define PARALLEL_BLKS 2
  assert(num_modes % 2 == 0 && "passing odd number of modes to get_rough_cost_for_2n_modes");
  const int width = cu_loc->width;
  const int height = cu_loc->height;
  cost_pixel_nxn_multi_func* satd_dual_func;
  cost_pixel_nxn_multi_func* sad_dual_func;
  satd_dual_func = uvg_pixels_get_satd_dual_func(width, height);
  sad_dual_func = uvg_pixels_get_sad_dual_func(width, height);


  uvg_pixel _preds[PARALLEL_BLKS * MIN(LCU_WIDTH, 64)* MIN(LCU_WIDTH, 64)+ SIMD_ALIGNMENT];
  pred_buffer preds = ALIGNED_POINTER(_preds, SIMD_ALIGNMENT);

  uvg_pixel _orig_block[MIN(LCU_WIDTH, 64) * MIN(LCU_WIDTH, 64) + SIMD_ALIGNMENT];
  uvg_pixel* orig_block = ALIGNED_POINTER(_orig_block, SIMD_ALIGNMENT);

  uvg_pixels_blit(orig, orig_block, width, height, orig_stride, width);
  
  const double mrl = state->encoder_control->cfg.mrl && (cu_loc->y % LCU_WIDTH) ? CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.multi_ref_line[0]), 1) : 0;
  const double not_mip = state->encoder_control->cfg.mip ? CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.mip_flag[mip_ctx]), 0) : 0;
  const double mip = state->encoder_control->cfg.mip ? CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.mip_flag[mip_ctx]), 1) : 0;
  double costs_out[PARALLEL_BLKS] = { 0 };
  double bits[PARALLEL_BLKS] = { 0 };
  for(int mode = 0; mode < num_modes; mode += PARALLEL_BLKS) {
    for (int i = 0; i < PARALLEL_BLKS; ++i) {
      uvg_intra_predict(state, &refs[search_data[mode + i].pred_cu.intra.multi_ref_idx], cu_loc, cu_loc, COLOR_Y, preds[i], &search_data[mode + i], NULL);
    }
    get_cost_dual(state, preds, orig_block, satd_dual_func, sad_dual_func, width, height, costs_out);

    for(int i = 0; i < PARALLEL_BLKS; ++i) {
      uint8_t multi_ref_idx = search_data[mode + i].pred_cu.intra.multi_ref_idx;
      if(multi_ref_idx) {
        bits[i] = mrl + not_mip;
        bits[i] += CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.multi_ref_line[1]), multi_ref_idx != 1);
        bits[i] += MIN(((mode + i) % 5) + 1, 4);
      }
      else if(search_data[mode + i].pred_cu.intra.mip_flag) {
        bits[i] = mip + 1;
        bits[i] += num_modes == 32 ? 4 : (num_modes == 16 ? 3 : (((mode + i) % 6) < 2 ? 2 : 3));
      }
      else {
        assert(0 && "get_rough_cost_for_2n_modes supports only mrl and mip mode cost calculation");
      }
    }
    search_data[mode].cost = costs_out[0];
    search_data[mode + 1].cost = costs_out[1];

    search_data[mode].cost += bits[0] * state->lambda_sqrt;
    search_data[mode + 1].cost += bits[1] * state->lambda_sqrt;
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
  int modes_to_check,
  intra_search_data_t *search_data,
  lcu_t *lcu,
  enum uvg_tree_type tree_type,
  const cu_loc_t* const cu_loc)
{
  const int width = cu_loc->width;
  const int height = cu_loc->height; // TODO: height for non-square blocks
  
  for (int mode = 0; mode < modes_to_check; mode++) {
    bool can_do_isp_search = search_data[mode].pred_cu.intra.mip_flag ? false : true; // Cannot use ISP with MIP
    // can_do_isp_search = search_data[mode].pred_cu.intra.multi_ref_idx == 0 ? can_do_isp_search : false; // Cannot use ISP with MRL
    const uint8_t mrl_idx = search_data[mode].pred_cu.intra.multi_ref_idx;
    double best_isp_cost = MAX_DOUBLE;
    double best_bits = MAX_DOUBLE;
    int8_t best_isp_mode = 0;
    int max_isp_modes = can_do_isp_search && uvg_can_use_isp(width, height) && state->encoder_control->cfg.isp ? NUM_ISP_MODES : 1;

    //
    uint8_t best_mts_mode_for_isp[NUM_ISP_MODES] = {0};
    uint8_t best_lfnst_mode_for_isp[NUM_ISP_MODES] = {0};
    for (int isp_mode = 0; isp_mode < max_isp_modes; ++isp_mode) {
       

      search_data[mode].pred_cu.intra.isp_mode = isp_mode;
      search_data[mode].pred_cu.intra.multi_ref_idx = isp_mode == ISP_MODE_NO_ISP ? mrl_idx : 0;
      double rdo_bitcost = uvg_luma_mode_bits(state, &search_data[mode].pred_cu, cu_loc, lcu);
      search_data[mode].pred_cu.tr_idx = MTS_TR_NUM;
      search_data[mode].bits = rdo_bitcost;
      search_data[mode].cost = rdo_bitcost * state->lambda;

      double mode_cost = search_intra_trdepth(state, cu_loc, MAX_INT, &search_data[mode], lcu, tree_type);
      best_mts_mode_for_isp[isp_mode] = search_data[mode].pred_cu.tr_idx;
      best_lfnst_mode_for_isp[isp_mode] = search_data[mode].pred_cu.lfnst_idx;
      search_data[mode].cost += mode_cost;
      if (search_data[mode].cost < best_isp_cost) {
        best_isp_cost = search_data[mode].cost;
        best_isp_mode = isp_mode;
        best_bits = search_data[mode].bits;
      }
      if (state->encoder_control->cfg.intra_rdo_et && !cbf_is_set_any(search_data[mode].pred_cu.cbf)) {
        modes_to_check = mode + 1;
        break;
      }
    }
    search_data[mode].cost = best_isp_cost;
    search_data[mode].bits = best_bits;
    search_data[mode].pred_cu.intra.isp_mode = best_isp_mode;
    search_data[mode].pred_cu.intra.multi_ref_idx = best_isp_mode == ISP_MODE_NO_ISP ? mrl_idx : 0;
    search_data[mode].pred_cu.tr_idx = best_mts_mode_for_isp[best_isp_mode];
    search_data[mode].pred_cu.tr_skip = best_mts_mode_for_isp[best_isp_mode] == MTS_SKIP;
    search_data[mode].pred_cu.lfnst_idx = best_lfnst_mode_for_isp[best_isp_mode];
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


double uvg_luma_mode_bits(const encoder_state_t *state, const cu_info_t* const cur_cu, const cu_loc_t*
                          const cu_loc,
                          const lcu_t* lcu)
{
  cabac_data_t* cabac = (cabac_data_t *)&state->search_cabac;
  double mode_bits = 0;
  cabac_data_t cabac_copy;
  memcpy(&cabac_copy, cabac, sizeof cabac_copy);
  uvg_encode_intra_luma_coding_unit(
    state,
    &cabac_copy, cur_cu,
    cu_loc, lcu, &mode_bits
    );

  return mode_bits;
}


double uvg_chroma_mode_bits(const encoder_state_t *state, int8_t chroma_mode, int8_t luma_mode)
{
  cabac_data_t* cabac = (cabac_data_t*)&state->search_cabac;
  const cabac_ctx_t* ctx;
  double mode_bits = 0;

  if (state->encoder_control->cfg.cclm) {
    ctx = &(cabac->ctx.cclm_flag);
    mode_bits += CTX_ENTROPY_FBITS(ctx, chroma_mode > 67);
  }
  ctx = &(cabac->ctx.chroma_pred_model);
  if (chroma_mode == luma_mode) {
    mode_bits += CTX_ENTROPY_FBITS(ctx, 0);
  } else {
    if(chroma_mode < 67) {
      mode_bits += 2.0 + CTX_ENTROPY_FBITS(ctx, 1);
    }
    else {
      ctx = &(cabac->ctx.cclm_model);
      mode_bits += CTX_ENTROPY_FBITS(ctx, chroma_mode != 81);
      if (chroma_mode != 81) mode_bits += 1;
    }
  }

  if(cabac->update) {
    if(chroma_mode != luma_mode) {
      // Again it does not matter what we actually write here
      CABAC_BINS_EP(cabac, 0, 2, "intra_chroma_pred_mode");      
    }
  }

  return mode_bits;
}

int8_t uvg_search_intra_chroma_rdo(
  encoder_state_t * const state,
  int8_t num_modes,
  lcu_t *const lcu,
  const cu_loc_t* const cu_loc,
  intra_search_data_t* chroma_data,
  int8_t luma_mode,
  enum uvg_tree_type tree_type,
  bool is_separate)
{
  const bool reconstruct_chroma = true;
  
  const int chroma_width  = cu_loc->chroma_width;
  const int chroma_height = cu_loc->chroma_height;
  uvg_intra_references refs[2];
  const vector2d_t luma_px = { cu_loc->x, cu_loc->y };
  const vector2d_t pic_px = {
    state->tile->frame->width,
    state->tile->frame->height,
  };


  if (reconstruct_chroma) {
    uvg_intra_build_reference(state, cu_loc, cu_loc, COLOR_U, &luma_px, &pic_px, lcu, &refs[0], state->encoder_control->cfg.wpp, NULL, 0, 0);
    uvg_intra_build_reference(state, cu_loc, cu_loc, COLOR_V, &luma_px, &pic_px, lcu, &refs[1], state->encoder_control->cfg.wpp, NULL, 0, 0);
    
    const vector2d_t lcu_px = { cu_loc->local_x, cu_loc->local_y };
    cabac_data_t temp_cabac;
    memcpy(&temp_cabac, &state->search_cabac, sizeof(cabac_data_t));
    
    const int offset = ((cu_loc->local_x) >> 1) + ((cu_loc->local_y) >> 1)* LCU_WIDTH_C;

    int lfnst_modes_to_check[3];
    if((is_separate || tree_type == UVG_CHROMA_T) && state->encoder_control->cfg.lfnst && PU_IS_TU(&chroma_data->pred_cu) && chroma_height >= 4 && chroma_width >= 4) {
      for (int i = 0; i < 3; ++i) {
        lfnst_modes_to_check[i] = i;
      }
    }
    else {
      lfnst_modes_to_check[0] = 0;
      lfnst_modes_to_check[1] = -1;
      lfnst_modes_to_check[2] = -1;
    }
    ALIGNED(64) uvg_pixel u_pred[LCU_WIDTH_C * LCU_WIDTH_C];
    ALIGNED(64) uvg_pixel v_pred[LCU_WIDTH_C * LCU_WIDTH_C];
    ALIGNED(64) int16_t u_resi[LCU_WIDTH_C * LCU_WIDTH_C];
    ALIGNED(64) int16_t v_resi[LCU_WIDTH_C * LCU_WIDTH_C];

    double original_c_lambda = state->c_lambda;
    state->quant_blocks[2].needs_init = true;
    state->rate_estimator[1].needs_init = true;

    for (int8_t mode_i = 0; mode_i < num_modes; ++mode_i) {
      const uint8_t mode = chroma_data[mode_i].pred_cu.intra.mode_chroma;
      double mode_bits = uvg_chroma_mode_bits(state, mode, luma_mode);
      chroma_data[mode_i].cost = mode_bits * state->c_lambda;
      chroma_data[mode_i].bits = mode_bits;
      cu_info_t* pred_cu = &chroma_data[mode_i].pred_cu;
      uint8_t best_lfnst_index = 0;
      for (int lfnst_i = 0; lfnst_i < 3; ++lfnst_i) {
        const int lfnst = lfnst_modes_to_check[lfnst_i];
        if (lfnst == -1) {
          continue;
        }
        state->c_lambda = original_c_lambda * (state->encoder_control->cfg.jccr && state->qp > 18 ? 1.3 : 1.0);
        pred_cu->cr_lfnst_idx = lfnst;
        chroma_data[mode_i].lfnst_costs[lfnst] += mode_bits * state->c_lambda;
        if (PU_IS_TU(pred_cu) && (tree_type != UVG_CHROMA_T || (pred_cu->log2_chroma_width < 5 && pred_cu->log2_chroma_height < 5))) {
          uvg_intra_predict(
            state,
            &refs[COLOR_U - 1],
            cu_loc,
            cu_loc,
            COLOR_U,
            u_pred,
            &chroma_data[mode_i],
            lcu);
          uvg_intra_predict(
            state,
            &refs[COLOR_V - 1],
            cu_loc,
            cu_loc,
            COLOR_V,
            v_pred,
            &chroma_data[mode_i],
            lcu);
          uvg_generate_residual(
            &lcu->ref.u[offset],
            u_pred,
            u_resi,
            chroma_width,
            chroma_height,
            LCU_WIDTH_C,
            chroma_width);
          uvg_generate_residual(
            &lcu->ref.v[offset],
            v_pred,
            v_resi,
            chroma_width,
            chroma_height,
            LCU_WIDTH_C,
            chroma_width);
          uvg_chorma_ts_out_t chorma_ts_out;
          uvg_chroma_transform_search(
            state,
            lcu,
            &temp_cabac,
            cu_loc,
            offset,
            pred_cu,
            u_pred,
            v_pred,
            u_resi,
            v_resi,
            &chorma_ts_out,
            is_separate ? UVG_CHROMA_T : tree_type);

          // LFNST constraint failed
          if(chorma_ts_out.best_u_index == -1 && chorma_ts_out.best_combined_index == -1) {
            chroma_data[mode_i].lfnst_costs[lfnst] = MAX_DOUBLE;
            continue;
          }

          double actual_cost = state->lambda * (chorma_ts_out.u_bits + chorma_ts_out.v_bits + mode_bits) + (chorma_ts_out.u_distortion + chorma_ts_out.v_distortion);
          if(chorma_ts_out.best_u_cost + chorma_ts_out.best_v_cost < chorma_ts_out.best_combined_cost) {
            chroma_data[mode_i].lfnst_costs[lfnst] = actual_cost;
            if( chroma_data[mode_i].lfnst_costs[lfnst] 
                < chroma_data[mode_i].lfnst_costs[best_lfnst_index] || lfnst_i == 0) {
              chroma_data[mode_i].pred_cu.joint_cb_cr = 0;
              chroma_data[mode_i].pred_cu.tr_skip &= 1;
              chroma_data[mode_i].pred_cu.tr_skip |= (chorma_ts_out.best_u_index == CHROMA_TS) << COLOR_U;
              chroma_data[mode_i].pred_cu.tr_skip |= (chorma_ts_out.best_v_index == CHROMA_TS) << COLOR_V;
              best_lfnst_index = lfnst;
              chroma_data[mode_i].cost = chroma_data[mode_i].lfnst_costs[lfnst];
            }
          }
          else {
            chroma_data[mode_i].lfnst_costs[lfnst] = actual_cost;
            if (chroma_data[mode_i].lfnst_costs[lfnst]
              < chroma_data[mode_i].lfnst_costs[best_lfnst_index] || lfnst_i == 0) {
              chroma_data[mode_i].pred_cu.joint_cb_cr = chorma_ts_out.best_combined_index;
              chroma_data[mode_i].pred_cu.tr_skip &= 1;
              best_lfnst_index = lfnst;
              chroma_data[mode_i].cost = chroma_data[mode_i].lfnst_costs[lfnst];
            }
          }

        }
        else {
          state->search_cabac.update = 1;
          chroma_data[mode_i].cost = mode_bits * state->c_lambda;
          uvg_intra_recon_cu(state,
                             &chroma_data[mode_i], cu_loc,
                             pred_cu, lcu,
                             tree_type,
                             false,
                             true);
          chroma_data[mode_i].cost += uvg_cu_rd_cost_chroma(state, pred_cu, lcu, cu_loc);
          memcpy(&state->search_cabac, &temp_cabac, sizeof(cabac_data_t));
        }
      }
      
      pred_cu->cr_lfnst_idx = best_lfnst_index;
    }
    sort_modes(chroma_data, num_modes);
    
    state->c_lambda = original_c_lambda;
    return chroma_data[0].pred_cu.intra.mode_chroma;
  }

  return 100;
}

#undef IS_JCCR_MODE

int8_t uvg_search_cu_intra_chroma(
  encoder_state_t * const state,
  const cu_loc_t* const cu_loc,
  lcu_t *lcu,
  intra_search_data_t *search_data,
  int8_t luma_mode,
  enum uvg_tree_type tree_type,
  bool is_separate)
{

  const cu_info_t *cur_pu = &search_data->pred_cu;
  
  int8_t modes[8] = { 0, 50, 18, 1, luma_mode, 81, 82, 83 };
  uint8_t total_modes = (state->encoder_control->cfg.cclm && uvg_cclm_is_allowed(state, cu_loc, cur_pu, tree_type) ? 8 : 5);
  for(int i = 0; i < 4; i++) {
    if (modes[i] == luma_mode) {
      modes[i] = 66;
      break;
    }
  }
  

  // The number of modes to select for slower chroma search. Luma mode
  // is always one of the modes, so 2 means the final decision is made
  // between luma mode and one other mode that looks the best
  // according to search_intra_chroma_rough.
  // const int8_t modes_in_depth[5] = { 1, 1, 1, 1, 2 };
  int num_modes = 1;

  if (state->encoder_control->cfg.rdo >= 2 || tree_type == UVG_CHROMA_T) {
    num_modes = total_modes;
  }

  intra_search_data_t chroma_data[8];
  FILL(chroma_data, 0);
  for (int i = 0; i < num_modes; i++) {
    chroma_data[i].pred_cu = *cur_pu;
    chroma_data[i].pred_cu.intra.mode_chroma = num_modes == 1 ? luma_mode : modes[i];
    chroma_data[i].cost = 0;
    if(!is_separate && tree_type == UVG_BOTH_T) {
      memcpy(chroma_data[i].lfnst_costs, search_data->lfnst_costs, sizeof(double) * 3);
    }
    chroma_data[i].best_isp_cbfs = search_data->best_isp_cbfs; //Copy isp cbfs so they don't get overwriten later
  }
  // Don't do rough mode search if all modes are selected.
  // FIXME: It might make more sense to only disable rough search if
  // num_modes is 0.is 0.

  if(state->encoder_control->cfg.cclm && 0){
    

    num_modes = search_intra_chroma_rough(state, chroma_data, lcu, luma_mode,
                                          tree_type,
                                          cu_loc);
  }
  
  if (num_modes > 1 || state->encoder_control->cfg.jccr) {
    uvg_search_intra_chroma_rdo(state, num_modes, lcu, cu_loc, chroma_data, luma_mode, tree_type, is_separate);
  }
  else if(cur_pu->lfnst_idx) {
    chroma_data[0].pred_cu.cr_lfnst_idx = cur_pu->lfnst_idx;
  }
  *search_data = chroma_data[0];
  return chroma_data[0].pred_cu.intra.mode_chroma;
}


static int select_candidates_for_further_search(const encoder_state_t * const state,
  intra_search_data_t *search_data,
  uint8_t regular_modes,
  uint8_t mip_modes,
  int width,
  int height
)
{
  const double threshold_cost = 1.0 + 1.4 / sqrt(width * height);
  const int max_cand_per_type = regular_modes >> 1;
  bool keepOneMip = search_data[regular_modes - 1].cost < search_data[regular_modes].cost;
  const int maxNumConv = 3;

  intra_search_data_t temp_mip_modes[3];
  const int transp_offset = mip_modes / 2;
  for(int i = 0; i <3; i++) {
    const bool     is_transp = search_data[regular_modes + i].cost > search_data[regular_modes + i + transp_offset].cost;
    temp_mip_modes[i] = search_data[regular_modes + i + (is_transp ? transp_offset : 0)];
  }
  sort_modes(search_data, regular_modes + mip_modes);
  const double minCost = search_data[0].cost;
  
  intra_search_data_t temp_list_out[9];
  int selected_modes = 0;
  int numConv = 0;
  int numMip = 0;
  for (int idx = 0; idx < regular_modes + keepOneMip; idx++)
  {
    bool addMode = false;

    if (!search_data[idx].pred_cu.intra.mip_flag)
    {
      addMode = (numConv < maxNumConv);
      numConv += addMode ? 1 : 0;
    }
    else
    {
      addMode = (numMip < max_cand_per_type || (search_data[idx].cost < threshold_cost * minCost) || keepOneMip);
      keepOneMip = false;
      numMip += addMode ? 1 : 0;
    }
    if (addMode)
    {
      temp_list_out[selected_modes++] = search_data[idx];
    }
  }

  if (width> 8 && height > 8)
  {
    // Sort MIP candidates by Hadamard cost
    // Append MIP mode to RD mode list
    for (int idx = 0; idx < 3; idx++)
    {
      bool alreadyIncluded = false;
      for (int list_idx = 0; list_idx < selected_modes; list_idx++)
      {
        if (temp_list_out[list_idx].pred_cu.intra.mip_flag &&
          temp_list_out[list_idx].pred_cu.intra.mip_is_transposed == temp_mip_modes[idx].pred_cu.intra.mip_is_transposed &&
          temp_list_out[list_idx].pred_cu.intra.mode == idx
          )
        {
          alreadyIncluded = true;
          break;
        }
      }

      if (!alreadyIncluded)
      {
        temp_list_out[selected_modes++] = temp_mip_modes[idx];
        // if (fastMip) break;
      }
    }
  }

  memcpy(search_data, temp_list_out, selected_modes * sizeof(intra_search_data_t));
  return selected_modes;
}



/**
 * Update lcu to have best modes at this depth.
 * \return Cost of best mode.
 */
void uvg_search_cu_intra(
  encoder_state_t * const state,
  intra_search_data_t* mode_out,
  lcu_t *lcu,
  enum uvg_tree_type tree_type,
  const cu_loc_t* const cu_loc)
{
  const vector2d_t lcu_px = { cu_loc->local_x, cu_loc->local_y };
  const int8_t log2_width = uvg_g_convert_to_log2[cu_loc->width];
  const int8_t log2_height = uvg_g_convert_to_log2[cu_loc->width];
  const vector2d_t luma_px = { cu_loc->x, cu_loc->y};
  const vector2d_t pic_px = { state->tile->frame->width, state->tile->frame->height };

  cu_info_t *cur_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x, lcu_px.y);

  uvg_intra_references refs[MAX_REF_LINE_IDX];

  int8_t candidate_modes[INTRA_MPM_COUNT];
  // Normal intra modes + mrl modes + mip modes
  intra_search_data_t search_data[UVG_NUM_INTRA_MODES +(MAX_REF_LINE_IDX - 1) * (INTRA_MPM_COUNT - 1) + 32];

  cu_info_t *left_cu = 0;
  cu_info_t *above_cu = 0;

  // Select left and top CUs if they are available.
  // Top CU is not available across LCU boundary.
  if (cu_loc->x >= SCU_WIDTH) {
    left_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x - 1, lcu_px.y+ cu_loc->height-1);
  }
  if (cu_loc->y >= SCU_WIDTH && lcu_px.y > 0) {
    above_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x+ cu_loc->width-1, lcu_px.y - 1);
  }
  int8_t num_cand = uvg_intra_get_dir_luma_predictor(cu_loc->x, cu_loc->y, candidate_modes, cur_cu, left_cu, above_cu);

  bool is_large = cu_loc->width > TR_MAX_WIDTH || cu_loc->height > TR_MAX_WIDTH;
  if (!is_large) {
    uvg_intra_build_reference(state, cu_loc, cu_loc, COLOR_Y, &luma_px, &pic_px, lcu, refs, state->encoder_control->cfg.wpp, NULL, 0, 0);
  }
  
  // This is needed for bit cost calculation and requires too many parameters to be
  // calculated inside the rough search functions
  uint8_t mip_ctx = uvg_get_mip_flag_context(cu_loc, lcu, NULL);

  // Find best intra mode for 2Nx2N.
  uvg_pixel *ref_pixels = &lcu->ref.y[lcu_px.x + lcu_px.y * LCU_WIDTH];

  // Need to set some data for all cus
  cu_info_t temp_pred_cu;
  temp_pred_cu = *cur_cu;
  temp_pred_cu.type = CU_INTRA;
  FILL(temp_pred_cu.intra, 0);
  // Find modes with multiple reference lines if in use. Do not use if CU in first row.
  uint8_t lines = state->encoder_control->cfg.mrl && lcu_px.y != 0 ? MAX_REF_LINE_IDX : 1;

  uint8_t number_of_modes;
  uint8_t num_regular_modes;
  bool skip_rough_search = (is_large || state->encoder_control->cfg.rdo >= 4);
  if (!skip_rough_search) {
    num_regular_modes = number_of_modes = search_intra_rough(
                          state,
                          cu_loc,
                          ref_pixels,
                          LCU_WIDTH,
                          refs,
                          cu_loc->width,
                          cu_loc->height,
                          candidate_modes,
                          search_data,
                          &temp_pred_cu,
                          mip_ctx);
    // if(lines == 1) sort_modes(search_data, number_of_modes);

  } else {
    for (int8_t i = 0; i < UVG_NUM_INTRA_MODES; i++) {
      search_data[i].pred_cu = temp_pred_cu;
      search_data[i].pred_cu.intra.mode = i;
      search_data[i].pred_cu.intra.mode_chroma = i;
      search_data[i].cost = MAX_INT;
    }
    number_of_modes = UVG_NUM_INTRA_MODES;
    num_regular_modes = UVG_NUM_INTRA_MODES;
  }

  uint8_t num_mrl_modes = 0;
  for(int line = 1; line < lines && !is_large; ++line) {
    uvg_pixel extra_refs[128 * MAX_REF_LINE_IDX] = { 0 };

    if (luma_px.x > 0 && lcu_px.x == 0 && lcu_px.y > 0) {
      videoframe_t* const frame = state->tile->frame;

      // Copy extra ref lines, including ref line 1 and top left corner.
      for (int i = 0; i < MAX_REF_LINE_IDX; ++i) {
        int height = (cu_loc->height) * 2 + MAX_REF_LINE_IDX;
        height = MIN(height, (LCU_WIDTH - lcu_px.y + MAX_REF_LINE_IDX)); // Cut short if on bottom LCU edge. Cannot take references from below since they don't exist.
        height = MIN(height, pic_px.y - luma_px.y + MAX_REF_LINE_IDX);
        uvg_pixels_blit(&frame->rec->y[(luma_px.y - MAX_REF_LINE_IDX) * frame->rec->stride + luma_px.x - (1 + i)],
          &extra_refs[i * 128],
          1, height,
          frame->rec->stride, 1);
      }
    }
    uvg_intra_build_reference(state, cu_loc, cu_loc, COLOR_Y, &luma_px, &pic_px, lcu, &refs[line], state->encoder_control->cfg.wpp, extra_refs, line, 0);
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
    get_rough_cost_for_2n_modes(state, refs, cu_loc,
                                ref_pixels,
                                LCU_WIDTH, search_data + number_of_modes, num_mrl_modes,
                                mip_ctx);
    sort_modes(search_data, number_of_modes + num_mrl_modes);
    number_of_modes = 6;
  }
  // number_of_modes += num_mrl_modes;
  // num_regular_modes += num_mrl_modes;

  int num_mip_modes = 0;
  if (state->encoder_control->cfg.mip) {
    // MIP is not allowed for 64 x 4 or 4 x 64 blocks
    if (!((cu_loc->height == 64 && cu_loc->width== 4) || (cu_loc->height== 4 && cu_loc->width == 64))) {
      num_mip_modes = NUM_MIP_MODES_FULL(cu_loc->width, cu_loc->height);

      for (int transpose = 0; transpose < 2; transpose++) {
        const int half_mip_modes = num_mip_modes / 2;
        for (int i = 0; i < half_mip_modes; ++i) {
          const int index = i + number_of_modes + transpose * half_mip_modes;
          search_data[index].pred_cu = temp_pred_cu;
          search_data[index].pred_cu.intra.mip_flag = 1;
          search_data[index].pred_cu.intra.mode = i;
          search_data[index].pred_cu.intra.mip_is_transposed = transpose;
          search_data[index].pred_cu.intra.mode_chroma = 0;
          search_data[index].cost = MAX_INT;
        }
      }
      if (!skip_rough_search) {
        get_rough_cost_for_2n_modes(state, refs, cu_loc,
          ref_pixels,
          LCU_WIDTH, search_data + number_of_modes, num_mip_modes,
          mip_ctx);
      }
    }
    number_of_modes += num_mip_modes;
  }

  // Refine results with slower search or get some results if rough search was skipped.
  const int32_t rdo_level = state->encoder_control->cfg.rdo;
  if (rdo_level >= 2 || skip_rough_search) {
    int number_of_modes_to_search;
    if (rdo_level == 4) {
      number_of_modes_to_search = number_of_modes;
    } else if (rdo_level == 2 || rdo_level == 3) {
      const uint8_t g_aucIntraModeNumFast_UseMPM_2D[7 - 2 + 1][7 - 2 + 1] =
      {
        {3, 3, 3, 3, 2, 2},  //   4x4,   4x8,   4x16,   4x32,   4x64,   4x128,
        {3, 3, 3, 3, 3, 2},  //   8x4,   8x8,   8x16,   8x32,   8x64,   8x128,
        {3, 3, 3, 3, 3, 2},  //  16x4,  16x8,  16x16,  16x32,  16x64,  16x128,
        {3, 3, 3, 3, 3, 2},  //  32x4,  32x8,  32x16,  32x32,  32x64,  32x128,
        {2, 3, 3, 3, 3, 2},  //  64x4,  64x8,  64x16,  64x32,  64x64,  64x128,
        {2, 2, 2, 2, 2, 3},  // 128x4, 128x8, 128x16, 128x32, 128x64, 128x128,
      };
      number_of_modes_to_search = g_aucIntraModeNumFast_UseMPM_2D[log2_width - 2][log2_height - 2];
    } else {
      // Check only the predicted modes.
      number_of_modes_to_search = 0;
    }
    if(!skip_rough_search) {
      if(state->encoder_control->cfg.mip) {
        number_of_modes_to_search = select_candidates_for_further_search(
          state,
          search_data,
          num_regular_modes,
          num_mip_modes,
          cu_loc->width,
          cu_loc->height
        );
      }
    }

    for(int pred_mode = 0; pred_mode < num_cand; ++pred_mode) {
      bool mode_found = false;
      for(int i = 0; i < number_of_modes_to_search; i++) {
        if(search_data[i].pred_cu.intra.mip_flag == 0 &&
          search_data[i].pred_cu.intra.multi_ref_idx == 0 &&
          search_data[i].pred_cu.intra.mode == candidate_modes[pred_mode]) {
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

    state->quant_blocks[0].needs_init = 1;
    state->rate_estimator[0].needs_init = 1;
    search_intra_rdo(
      state,
      number_of_modes_to_search,
      search_data,
      lcu,
      tree_type,
      cu_loc);
    search_data[0].pred_cu.mts_last_scan_pos = false;
    search_data[0].pred_cu.violates_mts_coeff_constraint = false;
  }

  *mode_out = search_data[0];
}
