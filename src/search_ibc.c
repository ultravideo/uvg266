/*****************************************************************************
 * This file is part of uvg266 VVC encoder.
 *
 * Copyright (c) 2022, Tampere University, ITU/ISO/IEC, project contributors
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

#include "search_ibc.h"
#include "search_inter.h"

#include <limits.h>
#include <stdlib.h>

#include "cabac.h"
#include "encoder.h"
#include "encode_coding_tree.h"
#include "image.h"
#include "imagelist.h"
#include "inter.h"
#include "uvg266.h"
#include "rdo.h"
#include "search.h"
#include "strategies/strategies-ipol.h"
#include "strategies/strategies-picture.h"
#include "transform.h"
#include "videoframe.h"

typedef struct {
  encoder_state_t *state;

  /**
   * \brief Current frame
   */
  const uvg_picture *pic;

  /**
   * \brief Top-left corner of the PU
   */
  vector2d_t origin;
  int32_t width;
  int32_t height;

  mv_t mv_cand[2][2];
  inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS];
  int32_t num_merge_cand;

  uvg_mvd_cost_func *mvd_cost_func;

  /**
   * \brief Possible optimized SAD implementation for the width, leave as
   *        NULL for arbitrary-width blocks
   */
  optimized_sad_func_ptr_t optimized_sad;

  lcu_t                   *lcu;

} ibc_search_info_t;





/**
 * \return  True if referred block is within current tile.
 */
static INLINE bool intmv_within_ibc_range(const ibc_search_info_t *info, int x, int y)
{  
  bool negative_values = x <= 0 && y <= 0;
  bool mv_range_valid = ((-y >= info->height) || (-x >= info->width)) && // Must be block height/width away from the block
                        SUB_SCU(info->origin.y) >= -y && // Y vector must be inside the current CTU
                        (-x <= IBC_BUFFER_WIDTH-LCU_WIDTH) && // X must be inside the buffer
                        info->origin.x + x >= 0; // Don't go outside of the frame


  return negative_values && mv_range_valid;
}

static INLINE bool fracmv_within_ibc_range(const ibc_search_info_t *info, int x, int y)
{
  return intmv_within_ibc_range(
    info,
    x >> INTERNAL_MV_PREC,
    y >> INTERNAL_MV_PREC);
}


static uint32_t calculate_ibc_cost_satd(const encoder_state_t *state, lcu_t* lcu, int32_t x, int32_t y, int32_t width, int32_t mv_x, int32_t mv_y)
{  
  const int x_scu    = SUB_SCU(x);
  const int y_scu    = SUB_SCU(y);

  cu_info_t *cur_cu    = LCU_GET_CU_AT_PX(lcu, x_scu, y_scu);

  cu_info_t cu_backup  = *cur_cu;
  uint32_t       cost  = MAX_INT;

  
  const uint32_t offset = x_scu + y_scu * LCU_WIDTH;
  const uint32_t offset_c = x_scu / 2 + y_scu / 2 * LCU_WIDTH_C;

  cur_cu->type    = CU_IBC;
  cur_cu->inter.mv_dir   = 1;
  cur_cu->skipped                         = false;
  cur_cu->merged                          = false;
  cur_cu->inter.mv_cand0                  = 0;
  cur_cu->joint_cb_cr                     = 0;
  cur_cu->inter.mv[0][0]                  = mv_x * (1 << INTERNAL_MV_PREC);;
  cur_cu->inter.mv[0][1]                  = mv_y * (1 << INTERNAL_MV_PREC);;

  uvg_inter_recon_cu(state, lcu, x, y, width, true, state->encoder_control->chroma_format != UVG_CSP_400);
  
  *cur_cu = cu_backup;

  cost = uvg_satd_any_size(width,
                           width,
                           lcu->rec.y + offset,
                           LCU_WIDTH,
                           &state->tile->frame->source->y[y * state->tile->frame->source->stride + x],
                           state->tile->frame->source->stride) >> (UVG_BIT_DEPTH - 8);

  if(state->encoder_control->chroma_format != UVG_CSP_400) {    
    cost += uvg_satd_any_size(width / 2,
                              width / 2,
                              lcu->rec.u + offset_c,
                              LCU_WIDTH_C,
                              &state->tile->frame->source->u[(y / 2) * (state->tile->frame->source->stride / 2) + (x / 2)],
                              state->tile->frame->source->stride / 2) >> (UVG_BIT_DEPTH - 8);
    cost += uvg_satd_any_size(width / 2,
                              width / 2,
                              lcu->rec.v + offset_c,
                              LCU_WIDTH_C,
                              &state->tile->frame->source->v[(y / 2) * (state->tile->frame->source->stride / 2) + (x / 2)],
                              state->tile->frame->source->stride / 2) >> (UVG_BIT_DEPTH - 8);
  }

  return cost;
}


static uint32_t calculate_ibc_cost_sad(const encoder_state_t *state, optimized_sad_func_ptr_t optimized_sad, lcu_t* lcu, int32_t x, int32_t y, int32_t width, int32_t mv_x, int32_t mv_y)
{  
  cu_info_t *cur_cu    = LCU_GET_CU_AT_PX(lcu, SUB_SCU(x), SUB_SCU(y));

  cu_info_t cu_backup  = *cur_cu;
  uint32_t       cost  = MAX_INT;

  const int x_scu    = SUB_SCU(x);
  const int y_scu    = SUB_SCU(y);
  const uint32_t offset = x_scu + y_scu * LCU_WIDTH;
  const uint32_t offset_c = x_scu / 2 + y_scu / 2 * LCU_WIDTH_C;

  cur_cu->type    = CU_IBC;
  cur_cu->inter.mv_dir   = 1;
  cur_cu->skipped                         = false;
  cur_cu->merged                          = false;
  cur_cu->inter.mv_cand0                  = 0;
  cur_cu->joint_cb_cr                     = 0;
  cur_cu->inter.mv[0][0]                  = mv_x * (1 << INTERNAL_MV_PREC);;
  cur_cu->inter.mv[0][1]                  = mv_y * (1 << INTERNAL_MV_PREC);;

  uvg_inter_recon_cu(state, lcu, x, y, width, true, state->encoder_control->chroma_format != UVG_CSP_400);
  
  *cur_cu = cu_backup;

  if (optimized_sad != NULL) {
    cost = optimized_sad(lcu->rec.y + offset, &state->tile->frame->source->y[y * state->tile->frame->source->stride + x], width, LCU_WIDTH, state->tile->frame->source->stride);
    if(state->encoder_control->chroma_format != UVG_CSP_400) {
      cost += optimized_sad(lcu->rec.u + offset_c, &state->tile->frame->source->u[(y / 2) * state->tile->frame->source->stride / 2 + x / 2], width / 2, LCU_WIDTH_C, state->tile->frame->source->stride / 2);
      cost += optimized_sad(lcu->rec.v + offset_c, &state->tile->frame->source->v[(y / 2) * state->tile->frame->source->stride / 2 + x / 2], width / 2, LCU_WIDTH_C, state->tile->frame->source->stride / 2);
    }
  } else {
    cost = uvg_reg_sad(lcu->rec.y + offset, &state->tile->frame->source->y[y * state->tile->frame->source->stride + x], width,width, LCU_WIDTH, state->tile->frame->source->stride);
    if(state->encoder_control->chroma_format != UVG_CSP_400) {
      cost += uvg_reg_sad(lcu->rec.u + offset_c, &state->tile->frame->source->u[(y / 2) * state->tile->frame->source->stride / 2 + x / 2], width / 2, width / 2, LCU_WIDTH_C, state->tile->frame->source->stride / 2);
      cost += uvg_reg_sad(lcu->rec.v + offset_c, &state->tile->frame->source->v[(y / 2) * state->tile->frame->source->stride / 2 + x / 2], width / 2, width / 2, LCU_WIDTH_C, state->tile->frame->source->stride / 2);
    }
  }

  return cost;
}

static bool check_mv_cost_satd(ibc_search_info_t *info,
                          int x,
                          int y,
                          double *best_cost,
                          double* best_bits,
                          vector2d_t *best_mv)
{

 }
/**
 * \brief Calculate cost for an integer motion vector.
 *
 * Updates best_mv, best_cost and best_bitcost to the new
 * motion vector if it yields a lower cost than the current one.
 *
 * If the motion vector violates the MV constraints for tiles or WPP, the
 * cost is not set.
 *
 * \return true if best_mv was changed, false otherwise
 */
static bool check_mv_cost(ibc_search_info_t *info,
                          int x,
                          int y,
                          double *best_cost,
                          double* best_bits,
                          vector2d_t *best_mv)
{
  if (!intmv_within_ibc_range(info, x, y)) return false;

  double bitcost = 0;
  double cost    = MAX_DOUBLE;

  cost = calculate_ibc_cost_sad(info->state, info->optimized_sad, info->lcu, info->origin.x, info->origin.y, info->width, x, y);

  if (cost >= *best_cost) return false;

  cost += info->mvd_cost_func(
      info->state,
      x, y, INTERNAL_MV_PREC,
      info->mv_cand,
      NULL,
      0,
      NULL,
      &bitcost
  );

  if (cost >= *best_cost) return false;

  // Set to motion vector in internal pixel precision.
  best_mv->x = x * (1 << INTERNAL_MV_PREC);
  best_mv->y = y * (1 << INTERNAL_MV_PREC);
  *best_cost = cost;
  *best_bits = bitcost;

  return true;
}


static unsigned get_ep_ex_golomb_bitcost(unsigned symbol)
{
  // Calculate 2 * log2(symbol )

  unsigned bins = 0;
  symbol += 0;
  if (symbol >= 1 << 8) { bins += 16; symbol >>= 8; }
  if (symbol >= 1 << 4) { bins += 8; symbol >>= 4; }
  if (symbol >= 1 << 2) { bins += 4; symbol >>= 2; }
  if (symbol >= 1 << 1) { bins += 2; }

  // TODO: It might be a good idea to put a small slope on this function to
  // make sure any search function that follows the gradient heads towards
  // a smaller MVD, but that would require fractinal costs and bits being
  // used everywhere in inter search.
  // return num_bins + 0.001 * symbol;

  return bins;
}


/**
 * \brief Checks if mv is one of the merge candidates.
 * \return true if found else return false
 */
static bool mv_in_merge(const ibc_search_info_t *info, vector2d_t mv)
{
  for (int i = 0; i < info->num_merge_cand; ++i) {
    if (info->merge_cand[i].dir == 3) continue;
    const vector2d_t merge_mv = {
      info->merge_cand[i].mv[info->merge_cand[i].dir - 1][0],
      info->merge_cand[i].mv[info->merge_cand[i].dir - 1][1]
    };
    if (merge_mv.x == mv.x * (1 << (INTERNAL_MV_PREC)) && merge_mv.y == mv.y * (1 << (INTERNAL_MV_PREC))) {
      return true;
    }
  }
  return false;
}


/**
 * \brief Select starting point for integer motion estimation search.
 *
 * Checks the zero vector, extra_mv and merge candidates and updates
 * best_mv to the best one.
 */
static void select_starting_point(ibc_search_info_t *info,
                                  vector2d_t extra_mv,
                                  double *best_cost,
                                  double* best_bits,
                                  vector2d_t *best_mv)
{
  // Check the 0-vector, so we can ignore all 0-vectors in the merge cand list.
  check_mv_cost(info, -info->width, 0, best_cost, best_bits, best_mv);

  // Change to integer precision.
  extra_mv.x >>= INTERNAL_MV_PREC;
  extra_mv.y >>= INTERNAL_MV_PREC;

  // Check mv_in if it's not one of the merge candidates.
  if ((extra_mv.x != 0 || extra_mv.y != 0) && !mv_in_merge(info, extra_mv)) {
    check_mv_cost(info, extra_mv.x, extra_mv.y, best_cost, best_bits, best_mv);
  }

  // Go through candidates
  for (int32_t i = 0; i < info->num_merge_cand; ++i) {
    int32_t x = (info->merge_cand[i].mv[info->merge_cand[i].dir - 1][0] + (1 << (INTERNAL_MV_PREC - 1)) ) >> INTERNAL_MV_PREC;
    int32_t y = (info->merge_cand[i].mv[info->merge_cand[i].dir - 1][1] + (1 << (INTERNAL_MV_PREC - 1)) ) >> INTERNAL_MV_PREC;

    check_mv_cost(info, x, y, best_cost, best_bits, best_mv);
  }
}

static double get_ibc_mvd_coding_cost(const encoder_state_t* state,
  const cabac_data_t* cabac,
  const int32_t mvd_hor,
  const int32_t mvd_ver)
{
  double bitcost = 4 << CTX_FRAC_BITS;
  const vector2d_t abs_mvd = { abs(mvd_hor), abs(mvd_ver) };
  bitcost += abs_mvd.x == 1 ? 1 << CTX_FRAC_BITS : (0 * (1 << CTX_FRAC_BITS));
  bitcost += abs_mvd.y == 1 ? 1 << CTX_FRAC_BITS : (0 * (1 << CTX_FRAC_BITS));

  bitcost += get_ep_ex_golomb_bitcost(abs_mvd.x) << CTX_FRAC_BITS;
  bitcost += get_ep_ex_golomb_bitcost(abs_mvd.y) << CTX_FRAC_BITS;

  // Round and shift back to integer bits.
  return bitcost / (1 << CTX_FRAC_BITS);
}


static int select_ibc_mv_cand(const encoder_state_t *state,
                          mv_t mv_cand[2][2],
                          int32_t mv_x,
                          int32_t mv_y,
                          double*cost_out)
{
  const bool same_cand =
    (mv_cand[0][0] == mv_cand[1][0] && mv_cand[0][1] == mv_cand[1][1]);

  if (same_cand && !cost_out) {
    // Pick the first one if both candidates are the same.
    return 0;
  }

  double (*mvd_coding_cost)(const encoder_state_t * const state,
                              const cabac_data_t*,
                              int32_t, int32_t);
  if (state->encoder_control->cfg.mv_rdo) {
    mvd_coding_cost = uvg_get_mvd_coding_cost_cabac;
  } else {
    mvd_coding_cost = get_ibc_mvd_coding_cost;
  }

  vector2d_t mvd = { mv_x - mv_cand[0][0], mv_y - mv_cand[0][1] };

  uvg_change_precision_vector2d(INTERNAL_MV_PREC, UVG_IMV_FPEL, &mvd);

  double cand1_cost = mvd_coding_cost(
      state, &state->cabac,
      mvd.x,
      mvd.y);

  double cand2_cost;
  if (same_cand) {
    cand2_cost = cand1_cost;
  } else {
    vector2d_t mvd2 = { mv_x - mv_cand[1][0], mv_y - mv_cand[1][1] };
    uvg_change_precision_vector2d(INTERNAL_MV_PREC, UVG_IMV_FPEL, &mvd2);
    cand2_cost = mvd_coding_cost(
      state, &state->cabac,
      mvd2.x,
      mvd2.y);
  }

  if (cost_out) {
    *cost_out = MIN(cand1_cost, cand2_cost);
  }

  // Pick the second candidate if it has lower cost.
  return cand2_cost < cand1_cost ? 1 : 0;
}

static double calc_ibc_mvd_cost(const encoder_state_t *state,
                            int x,
                            int y,
                            int mv_shift,
                            mv_t mv_cand[2][2],
                            inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS],
                            int16_t num_cand,
                            int32_t ref_idx,
                            double* bitcost)
{
  double temp_bitcost = 0;
  uint32_t merge_idx;
  int8_t merged      = 0;

  x *= 1 << mv_shift;
  y *= 1 << mv_shift;

  // Check every candidate to find a match
  for(merge_idx = 0; merge_idx < (uint32_t)num_cand; merge_idx++) {
    if (merge_cand[merge_idx].mv[merge_cand[merge_idx].dir - 1][0] == x &&
        merge_cand[merge_idx].mv[merge_cand[merge_idx].dir - 1][1] == y) {
      temp_bitcost += merge_idx;
      merged = 1;
      break;
    }
  }

  // Check mvd cost only if mv is not merged
  if (!merged) {
    double mvd_cost = 0;
    select_ibc_mv_cand(state, mv_cand, x, y, &mvd_cost);
    temp_bitcost += mvd_cost;
  }
  *bitcost = temp_bitcost;
  return temp_bitcost * state->lambda_sqrt;
}


static bool early_terminate(ibc_search_info_t *info,
                            double *best_cost,
                            double* best_bits,
                            vector2d_t *best_mv)
{
  static const vector2d_t small_hexbs[7] = {
      { 0, -1 }, { -1, 0 }, { 0, 1 }, { 1, 0 },
      { 0, -1 }, { -1, 0 }, { 0, 0 },
  };

  vector2d_t mv = { best_mv->x >> INTERNAL_MV_PREC, best_mv->y >> INTERNAL_MV_PREC };

  int first_index = 0;
  int last_index = 3;

  for (int k = 0; k < 2; ++k) {
    double threshold;
    if (info->state->encoder_control->cfg.me_early_termination ==
        UVG_ME_EARLY_TERMINATION_SENSITIVE)
    {
      threshold = *best_cost * 0.95;
    } else {
      threshold = *best_cost;
    }

    int best_index = 6;
    for (int i = first_index; i <= last_index; i++) {
      int x = mv.x + small_hexbs[i].x;
      int y = mv.y + small_hexbs[i].y;

      if (check_mv_cost(info, x, y, best_cost, best_bits, best_mv)) {
        best_index = i;
      }
    }

    // Adjust the movement vector
    mv.x += small_hexbs[best_index].x;
    mv.y += small_hexbs[best_index].y;

    // If best match is not better than threshold, we stop the search.
    if (*best_cost >= threshold) {
      return true;
    }

    first_index = (best_index + 3) % 4;
    last_index = first_index + 2;
  }
  return false;
}



/**
 * \brief Do motion search using the HEXBS algorithm.
 *
 * \param info      search info
 * \param extra_mv  extra motion vector to check
 * \param steps     how many steps are done at maximum before exiting, does not affect the final step
 *
 * Motion vector is searched by first searching iteratively with the large
 * hexagon pattern until the best match is at the center of the hexagon.
 * As a final step a smaller hexagon is used to check the adjacent pixels.
 *
 * If a non 0,0 predicted motion vector predictor is given as extra_mv,
 * the 0,0 vector is also tried. This is hoped to help in the case where
 * the predicted motion vector is way off. In the future even more additional
 * points like 0,0 might be used, such as vectors from top or left.
 */
static void hexagon_search(ibc_search_info_t *info,
                           vector2d_t extra_mv,
                           uint32_t steps,
                           double *best_cost,
                           double* best_bits,
                           vector2d_t *best_mv)
{
  // The start of the hexagonal pattern has been repeated at the end so that
  // the indices between 1-6 can be used as the start of a 3-point list of new
  // points to search.
  //   6--1,7
  //  /     \    =)
  // 5   0  2,8
  //  \     /
  //   4---3
  static const vector2d_t large_hexbs[9] = {
      { 0, 0 },
      { 1, -2 }, { 2, 0 }, { 1, 2 }, { -1, 2 }, { -2, 0 }, { -1, -2 },
      { 1, -2 }, { 2, 0 }
  };
  // This is used as the last step of the hexagon search.
  //   1
  // 2 0 3
  //   4
  static const vector2d_t small_hexbs[9] = {
      { 0, 0 },
      { 0, -1 }, { -1, 0 }, { 1, 0 }, { 0, 1 },
      { -1, -1 }, { 1, -1 }, { -1, 1 }, { 1, 1 }
  };

  vector2d_t mv = { best_mv->x >> INTERNAL_MV_PREC, best_mv->y >> INTERNAL_MV_PREC };

  // Current best index, either to merge_cands, large_hexbs or small_hexbs.
  int best_index = 0;

  // Search the initial 7 points of the hexagon.
  for (int i = 1; i < 7; ++i) {
    if (check_mv_cost(info, mv.x + large_hexbs[i].x, mv.y + large_hexbs[i].y, best_cost, best_bits, best_mv)) {
      best_index = i;
    }
  }

  // Iteratively search the 3 new points around the best match, until the best
  // match is in the center.
  while (best_index != 0 && steps != 0) {
    // decrement count if enabled
    if (steps > 0) steps -= 1;

    // Starting point of the 3 offsets to be searched.
    unsigned start;
    if (best_index == 1) {
      start = 6;
    } else if (best_index == 8) {
      start = 1;
    } else {
      start = best_index - 1;
    }

    // Move the center to the best match.
    mv.x += large_hexbs[best_index].x;
    mv.y += large_hexbs[best_index].y;
    best_index = 0;

    // Iterate through the next 3 points.
    for (int i = 0; i < 3; ++i) {
      vector2d_t offset = large_hexbs[start + i];
      if (check_mv_cost(info, mv.x + offset.x, mv.y + offset.y, best_cost, best_bits, best_mv)) {
        best_index = start + i;
      }
    }
  }

  // Move the center to the best match.
  //mv.x += large_hexbs[best_index].x;
  //mv.y += large_hexbs[best_index].y;

  // Do the final step of the search with a small pattern.
  for (int i = 1; i < 9; ++i) {
    check_mv_cost(info, mv.x + small_hexbs[i].x, mv.y + small_hexbs[i].y, best_cost, best_bits, best_mv);
  }
}

/**
* \brief Do motion search using the diamond algorithm.
*
* \param info      search info
* \param extra_mv  extra motion vector to check
* \param steps     how many steps are done at maximum before exiting
*
* Motion vector is searched by searching iteratively with a diamond-shaped
* pattern. We take care of not checking the direction we came from, but
* further checking for avoiding visits to already visited points is not done.
*
* If a non 0,0 predicted motion vector predictor is given as extra_mv,
* the 0,0 vector is also tried. This is hoped to help in the case where
* the predicted motion vector is way off. In the future even more additional
* points like 0,0 might be used, such as vectors from top or left.
**/
static void diamond_search(ibc_search_info_t *info,
                           vector2d_t extra_mv,
                           uint32_t steps,
                           double *best_cost,
                           double* best_bits,
                           vector2d_t *best_mv)
{
  enum diapos {
    DIA_UP = 0,
    DIA_RIGHT = 1,
    DIA_LEFT = 2,
    DIA_DOWN = 3,
    DIA_CENTER = 4,
  };

  // a diamond shape with the center included
  //   0
  // 2 4 1
  //   3
  static const vector2d_t diamond[5] = {
    {0, -1}, {1, 0}, {0, 1}, {-1, 0},
    {0, 0}
  };
  
  // current motion vector
  vector2d_t mv = { best_mv->x >> INTERNAL_MV_PREC, best_mv->y >> INTERNAL_MV_PREC };

  // current best index
  enum diapos best_index = DIA_CENTER;

  // initial search of the points of the diamond
  for (int i = 0; i < 5; ++i) {
    if (check_mv_cost(info, mv.x + diamond[i].x, mv.y + diamond[i].y, best_cost, best_bits, best_mv)) {
      best_index = i;
    }
  }

  if (best_index == DIA_CENTER) {
    // the center point was the best in initial check
    return;
  }

  // Move the center to the best match.
  mv.x += diamond[best_index].x;
  mv.y += diamond[best_index].y;

  // the arrival direction, the index of the diamond member that will be excluded
  enum diapos from_dir = DIA_CENTER;

  // whether we found a better candidate this iteration
  uint8_t better_found;

  do {
    better_found = 0;
    // decrement count if enabled
    if (steps > 0) steps -= 1;

    // search the points of the diamond
    for (int i = 0; i < 4; ++i) {
      // this is where we came from so it's checked already
      if (i == from_dir) continue;

      if (check_mv_cost(info, mv.x + diamond[i].x, mv.y + diamond[i].y, best_cost, best_bits, best_mv)) {
        best_index = i;
        better_found = 1;
      }
    }

    if (better_found) {
      // Move the center to the best match.
      mv.x += diamond[best_index].x;
      mv.y += diamond[best_index].y;

      // record where we came from to the next iteration
      // the xor operation flips the orientation
      from_dir = best_index ^ 0x3;
    }
  } while (better_found && steps != 0);
  // and we're done
}


/**
 * \brief Check if an identical merge candidate exists in a list
 *
 * \param all_cand        Full list of available merge candidates
 * \param cand_to_add     Merge candidate to be checked for duplicates
 * \param added_idx_list  List of indices of unique merge candidates
 * \param list_size       Size of the list
 *
 * \return                Does an identical candidate exist in list
 */
static bool merge_candidate_in_list(inter_merge_cand_t *all_cands,
                                    inter_merge_cand_t *cand_to_add,
                                    unit_stats_map_t *merge)
{
  bool found = false;
  for (int i = 0; i < merge->size && !found; ++i) {
    int key = merge->keys[i];
    inter_merge_cand_t * list_cand = &all_cands[merge->unit[key].merge_idx];

    found = 
        cand_to_add->mv[0][0] == list_cand->mv[0][0] &&
        cand_to_add->mv[0][1] == list_cand->mv[0][1];
  }

  return found;
}

/**
 * \brief Collect PU parameters and costs at this depth.
 *
 * \param state       encoder state
 * \param x_cu        x-coordinate of the containing CU
 * \param y_cu        y-coordinate of the containing CU
 * \param depth       depth of the CU in the quadtree
 * \param part_mode   partition mode of the CU
 * \param i_pu        index of the PU in the CU
 * \param lcu         containing LCU
 *
 * \param amvp        Return searched AMVP PUs sorted by costs
 * \param merge       Return searched Merge PUs sorted by costs
 */
static void search_pu_ibc(encoder_state_t * const state,
  int x_cu, int y_cu,
  int depth,
  part_mode_t part_mode,
  int i_pu,
  unit_stats_map_t *amvp,
  unit_stats_map_t *merge,
  ibc_search_info_t *info)
{
  const uvg_config *cfg = &state->encoder_control->cfg;
  const videoframe_t * const frame = state->tile->frame;
  const int width_cu = LCU_WIDTH >> depth;
  const int x = PU_GET_X(part_mode, width_cu, x_cu, i_pu);
  const int y = PU_GET_Y(part_mode, width_cu, y_cu, i_pu);
  const int width = PU_GET_W(part_mode, width_cu, i_pu);
  const int height = PU_GET_H(part_mode, width_cu, i_pu);

  // Merge candidate A1 may not be used for the second PU of Nx2N, nLx2N and
  // nRx2N partitions.
  const bool merge_a1 = i_pu == 0 || width >= height;
  // Merge candidate B1 may not be used for the second PU of 2NxN, 2NxnU and
  // 2NxnD partitions.
  const bool merge_b1 = i_pu == 0 || width <= height;


  lcu_t                     *lcu      = info->lcu;
  const int x_local = SUB_SCU(x);
  const int y_local = SUB_SCU(y);
  cu_info_t *cur_pu = LCU_GET_CU_AT_PX(lcu, x_local, y_local);
  cur_pu->type = CU_IBC;
  cur_pu->part_size = part_mode;
  cur_pu->depth = depth;
  cur_pu->tr_depth = depth;
  cur_pu->qp = state->qp;

  // Default to candidate 0
  CU_SET_MV_CAND(cur_pu, 0, 0);
  
  FILL(*info, 0);

  info->state          = state;
  info->pic            = frame->source;
  info->origin.x       = x;
  info->origin.y       = y;
  info->width          = width;
  info->height         = height;
  info->mvd_cost_func  = cfg->mv_rdo ? uvg_calc_ibc_mvd_cost_cabac : calc_ibc_mvd_cost;
  info->optimized_sad  = uvg_get_optimized_sad(width);
  info->lcu            = lcu;

  // Search for merge mode candidates
  info->num_merge_cand = uvg_inter_get_merge_cand(
                          state,
                          x, y,
                          width, height,
                          merge_a1, merge_b1,
                          info->merge_cand,
                          lcu);

  // Merge Analysis starts here
  merge->size = 0;
  for (int i = 0; i < MRG_MAX_NUM_CANDS; ++i) {
    merge->keys[i] = -1;
    merge->cost[i] = MAX_DOUBLE;
  }

  const double merge_flag_cost = CTX_ENTROPY_FBITS(&state->search_cabac.ctx.cu_merge_flag_ext_model, 1);
#ifdef COMPLETE_PRED_MODE_BITS
  // Technically counting these bits would be correct, however counting
  // them universally degrades quality so this block is disabled by default
  const double no_skip_flag = CTX_ENTROPY_FBITS(&state->search_cabac.ctx.cu_skip_flag_model[uvg_get_skip_context(x, y, lcu, NULL)], 0);
#else
  const double no_skip_flag = 0;
#endif
  // Check motion vector constraints and perform rough search
  for (int merge_idx = 0; merge_idx < info->num_merge_cand; ++merge_idx) {

    inter_merge_cand_t *cur_cand = &info->merge_cand[merge_idx];
    cur_pu->inter.mv_dir = cur_cand->dir;
    cur_pu->inter.mv[0][0] = cur_cand->mv[0][0];
    cur_pu->inter.mv[0][1] = cur_cand->mv[0][1];


    bool is_duplicate = merge_candidate_in_list(info->merge_cand, cur_cand, merge);

    // Don't try merge candidates that don't satisfy mv constraints.
    // Don't add duplicates to list
    if ((!fracmv_within_ibc_range(info, cur_pu->inter.mv[0][0], cur_pu->inter.mv[0][1])) ||
        is_duplicate)
    {
      continue;
    }
    uvg_inter_pred_pu(state, info->lcu, x_cu, y_cu, width_cu, true, false, i_pu);
    merge->unit[merge->size] = *cur_pu;
    merge->unit[merge->size].type = CU_IBC;
    merge->unit[merge->size].merge_idx = merge_idx;
    merge->unit[merge->size].merged = true;
    merge->unit[merge->size].skipped = false;

    double bits = merge_flag_cost + merge_idx + CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.cu_merge_idx_ext_model), merge_idx != 0);
    if(state->encoder_control->cfg.rdo >= 2 && cur_pu->part_size == SIZE_2Nx2N) {
      uvg_cu_cost_inter_rd2(state, x, y, depth, &merge->unit[merge->size], lcu, &merge->cost[merge->size], &bits);
    }
    else {
      merge->cost[merge->size] = uvg_satd_any_size(width, height,
        lcu->rec.y + y_local * LCU_WIDTH + x_local, LCU_WIDTH,
        lcu->ref.y + y_local * LCU_WIDTH + x_local, LCU_WIDTH);
      bits += no_skip_flag;
      merge->cost[merge->size] += bits * info->state->lambda_sqrt;
    }
    // Add cost of coding the merge index
    merge->bits[merge->size] = bits;
    merge->keys[merge->size] = merge->size;


    merge->size++;
  }

  assert(merge->size <= MAX_UNIT_STATS_MAP_SIZE);
  uvg_sort_keys_by_cost(merge);

  // Try early skip decision on just one merge candidate if available
  int num_rdo_cands = MIN(1, merge->size);
    
  // Early Skip Mode Decision
  bool has_chroma = state->encoder_control->chroma_format != UVG_CSP_400;
  if (cfg->early_skip && cur_pu->part_size == SIZE_2Nx2N) {
    for (int merge_key = 0; merge_key < num_rdo_cands; ++merge_key) {
      if(cfg->rdo >= 2 && merge->unit[merge->keys[merge_key]].skipped) {
        merge->size = 1;
        merge->bits[0] = merge->bits[merge->keys[merge_key]];
        merge->cost[0] = merge->cost[merge->keys[merge_key]];
        merge->unit[0] = merge->unit[merge->keys[merge_key]];
        merge->keys[0] = 0;
      }
      else if(cfg->rdo < 2) {
        // Reconstruct blocks with merge candidate.
        // Check luma CBF. Then, check chroma CBFs if luma CBF is not set
        // and chroma exists.
        // Early terminate if merge candidate with zero CBF is found.
        int merge_idx           = merge->unit[merge->keys[merge_key]].merge_idx;
        cur_pu->inter.mv_dir    = info->merge_cand[merge_idx].dir;
        cur_pu->inter.mv[0][0]  = info->merge_cand[merge_idx].mv[0][0];
        cur_pu->inter.mv[0][1]  = info->merge_cand[merge_idx].mv[0][1];
        uvg_lcu_fill_trdepth(lcu, x, y, depth, MAX(1, depth), UVG_BOTH_T);
        uvg_inter_recon_cu(state, lcu, x, y, width, true, false);
        uvg_quantize_lcu_residual(state, true, false, false, x, y, depth, cur_pu, lcu, true, UVG_BOTH_T);

        if (cbf_is_set(cur_pu->cbf, depth, COLOR_Y)) {
          continue;
        }
        else if (has_chroma) {
          uvg_inter_recon_cu(state, lcu, x, y, width, false, has_chroma);
          uvg_quantize_lcu_residual(state, false, has_chroma, 
            false, /*we are only checking for lack of coeffs so no need to check jccr*/
            x, y, depth, cur_pu, lcu, true, UVG_BOTH_T);
          if (!cbf_is_set_any(cur_pu->cbf, depth)) {
            cur_pu->type = CU_IBC;
            cur_pu->merge_idx = merge_idx;
            cur_pu->skipped = true;

            merge->size = 1;
            merge->cost[0] = 0.0; // TODO: Check this
            merge->bits[0] = merge_idx; // TODO: Check this
            merge->unit[0] = *cur_pu;
            return;
          }
        }
      }
    }
  }

  // AMVP search starts here
  amvp[0].size    = 0;
  amvp[1].size    = 0;
  amvp[2].size    = 0;
  amvp[0].cost[0] = MAX_DOUBLE;


  // Do the motion search

  uvg_inter_get_mv_cand(info->state,
    info->origin.x,
    info->origin.y,
    info->width,
    info->height,
    info->mv_cand,
    cur_pu,
    lcu,
    NULL);

  vector2d_t best_mv = { 0, 0 };

  double best_cost = MAX_DOUBLE;
  double best_bits = MAX_INT;

  // Select starting point from among merge candidates. These should
  // include both mv_cand vectors and (0, 0).
  select_starting_point(info, best_mv, &best_cost, &best_bits, &best_mv);
  bool skip_me = early_terminate(info, &best_cost, &best_bits, &best_mv);
      
  if (!(info->state->encoder_control->cfg.me_early_termination && skip_me)) {

    switch (cfg->ime_algorithm) {
      case UVG_IME_DIA:
        diamond_search(info, best_mv, info->state->encoder_control->cfg.me_max_steps,
                       &best_cost, &best_bits, &best_mv);
        break;
      default:
        hexagon_search(info, best_mv, info->state->encoder_control->cfg.me_max_steps,
                       &best_cost, &best_bits, &best_mv);
        break;
    }
  }

  if (best_cost < MAX_DOUBLE) {
    // Recalculate inter cost with SATD.
    best_cost = calculate_ibc_cost_satd(
      info->state,
      lcu,
      info->origin.x,
      info->origin.y,
      info->width,
      (best_mv.x >> INTERNAL_MV_PREC),
      (best_mv.y >> INTERNAL_MV_PREC));
    best_cost += best_bits * info->state->lambda_sqrt;
  }


  int cu_mv_cand = select_ibc_mv_cand(info->state, info->mv_cand, best_mv.x, best_mv.y, NULL);

  // Update best unipreds for biprediction
  bool valid_mv = fracmv_within_ibc_range(info, best_mv.x, best_mv.y);
  if (valid_mv && best_cost < MAX_DOUBLE) {

    // Map reference index to L0/L1 pictures
    unit_stats_map_t *cur_map = &amvp[0];
    int entry = cur_map->size;
    cu_info_t *unipred_pu = &cur_map->unit[entry];
    *unipred_pu = *cur_pu;
    unipred_pu->type = CU_IBC;
    unipred_pu->merged  = false;
    unipred_pu->skipped = false;
    unipred_pu->inter.mv_dir = 1;
    unipred_pu->inter.mv[0][0] = (mv_t)best_mv.x;
    unipred_pu->inter.mv[0][1] = (mv_t)best_mv.y;
    CU_SET_MV_CAND(unipred_pu, 0, cu_mv_cand);

    cur_map->cost[entry] = best_cost;
    cur_map->bits[entry] = best_bits;
    cur_map->keys[entry] = entry;
    cur_map->size++;
  }


  assert(amvp[0].size <= MAX_UNIT_STATS_MAP_SIZE);
  uvg_sort_keys_by_cost(&amvp[0]);

  int best_keys[2] = { 
    amvp[0].size > 0 ? amvp[0].keys[0] : 0, 
    amvp[1].size > 0 ? amvp[1].keys[0] : 0
  };

  cu_info_t *best_unipred[2] = {
    &amvp[0].unit[best_keys[0]],
    &amvp[1].unit[best_keys[1]]
  };


  if (state->encoder_control->cfg.rdo >= 2 && cur_pu->part_size == SIZE_2Nx2N) {
    if (amvp[0].size) uvg_cu_cost_inter_rd2(state, x, y, depth, &amvp[0].unit[best_keys[0]], lcu, &amvp[0].cost[best_keys[0]], &amvp[0].bits[best_keys[0]]);    
  }


  if(cfg->rdo < 2) {
    int predmode_ctx;
    const int skip_contest = uvg_get_skip_context(x, y, lcu, NULL, &predmode_ctx);
    const double no_skip_flag = CTX_ENTROPY_FBITS(&state->search_cabac.ctx.cu_skip_flag_model[skip_contest], 0);

    const double pred_mode_bits = CTX_ENTROPY_FBITS(&state->search_cabac.ctx.cu_pred_mode_model[predmode_ctx], 0);
    const double total_bits = no_skip_flag + pred_mode_bits;
    if(amvp[0].size > 0) {
      const uint8_t best_key = amvp[0].keys[0];
      amvp[0].bits[best_key] += total_bits;
      amvp[0].cost[best_key] += (total_bits)* state->lambda_sqrt;
    }
  }
}



/**
 * \brief Update CU to have best modes at this depth.
 *
 * Only searches the 2Nx2N partition mode.
 *
 * \param state       encoder state
 * \param x           x-coordinate of the CU
 * \param y           y-coordinate of the CU
 * \param depth       depth of the CU in the quadtree
 * \param lcu         containing LCU
 *
 * \param inter_cost    Return inter cost
 * \param inter_bitcost Return inter bitcost
 */
void uvg_search_cu_ibc(encoder_state_t * const state,
                         int x, int y, int depth,
                         lcu_t *lcu,
                         double   *inter_cost,
                         double* inter_bitcost)
{
  *inter_cost = MAX_DOUBLE;
  *inter_bitcost = MAX_INT;

  // Store information of L0, L1, and bipredictions.
  // Best cost will be left at MAX_DOUBLE if no valid CU is found.
  // These will be initialized by the following function.
  unit_stats_map_t amvp[3];
  unit_stats_map_t merge;
  ibc_search_info_t info;

  info.lcu = lcu;

  search_pu_ibc(state,
                  x, y, depth,
                  SIZE_2Nx2N, 0,
                  amvp,
                  &merge,
                  &info);

  // Early Skip CU decision
  if (merge.size == 1 && merge.unit[0].skipped) {
    *inter_cost    = merge.cost[0];
    *inter_bitcost = merge.bits[0];
    return;
  }

  cu_info_t *best_inter_pu = NULL;

  // Find best AMVP PU
  for (int mv_dir = 1; mv_dir < 4; ++mv_dir) {

    int best_key = amvp[mv_dir - 1].keys[0];

    if (amvp[mv_dir - 1].size > 0 &&
        amvp[mv_dir - 1].cost[best_key] < *inter_cost) {

      best_inter_pu  = &amvp[mv_dir - 1].unit[best_key];
      *inter_cost    =  amvp[mv_dir - 1].cost[best_key];
      *inter_bitcost =  amvp[mv_dir - 1].bits[best_key];
    }
  }

  // Compare best AMVP against best Merge mode
  int best_merge_key = merge.keys[0];

  if (merge.size > 0 && merge.cost[best_merge_key] < *inter_cost) {

    best_inter_pu  = &merge.unit[best_merge_key];
    *inter_cost    =  merge.cost[best_merge_key];
    *inter_bitcost =  0; // TODO: Check this
  }

  if (*inter_cost == MAX_DOUBLE) {
    // Could not find any motion vector.
    *inter_cost = MAX_DOUBLE;
    *inter_bitcost = MAX_INT;
    return;
  }

  const int x_local = SUB_SCU(x);
  const int y_local = SUB_SCU(y);
  cu_info_t *cur_pu = LCU_GET_CU_AT_PX(lcu, x_local, y_local);
  *cur_pu = *best_inter_pu;
  cur_pu->type       = CU_IBC;

  uvg_inter_recon_cu(state, lcu, x, y, CU_WIDTH_FROM_DEPTH(depth),
    true, state->encoder_control->chroma_format != UVG_CSP_400);   

  if (*inter_cost < MAX_DOUBLE) {
    assert(fracmv_within_ibc_range(&info, cur_pu->inter.mv[0][0], cur_pu->inter.mv[0][1]));
  }
}
