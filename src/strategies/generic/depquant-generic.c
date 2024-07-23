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

#include "strategies/generic/depquant-generic.h"

#include "dep_quant.h"

#include "cu.h"
#include "encoderstate.h"
#include "intra.h"
#include "rdo.h"
#include "strategyselector.h"
#include "transform.h"
#include "uvg_math.h"
#include "generic/quant-generic.h"
static const int32_t g_goRiceBits[4][RICEMAX] = {
  {32768,  65536,  98304,  131072, 163840, 196608, 262144, 262144,
   327680, 327680, 327680, 327680, 393216, 393216, 393216, 393216,
   393216, 393216, 393216, 393216, 458752, 458752, 458752, 458752,
   458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752},
  {65536,  65536,  98304,  98304,  131072, 131072, 163840, 163840,
   196608, 196608, 229376, 229376, 294912, 294912, 294912, 294912,
   360448, 360448, 360448, 360448, 360448, 360448, 360448, 360448,
   425984, 425984, 425984, 425984, 425984, 425984, 425984, 425984},
  {98304,  98304,  98304,  98304,  131072, 131072, 131072, 131072,
   163840, 163840, 163840, 163840, 196608, 196608, 196608, 196608,
   229376, 229376, 229376, 229376, 262144, 262144, 262144, 262144,
   327680, 327680, 327680, 327680, 327680, 327680, 327680, 327680},
  {131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072,
   163840, 163840, 163840, 163840, 163840, 163840, 163840, 163840,
   196608, 196608, 196608, 196608, 196608, 196608, 196608, 196608,
   229376, 229376, 229376, 229376, 229376, 229376, 229376, 229376},
};


static INLINE void checkRdCostSkipSbbZeroOut(
  Decision* decision, 
  const all_depquant_states* const state,
  int decision_id, 
  int skip_offset) {
  int64_t rdCost = state->m_rdCost[decision_id + skip_offset] + state->m_sbbFracBits[decision_id + skip_offset][0];
  decision->rdCost[decision_id] = rdCost;
  decision->absLevel[decision_id] = 0;
  decision->prevId[decision_id] = 4 + state->m_stateId[decision_id + skip_offset];
}

static INLINE void checkRdCostSkipSbb(const all_depquant_states* const state, Decision * decisions, int decision_id, int skip_offset)
{
  int64_t rdCost = state->m_rdCost[skip_offset + decision_id] + state->m_sbbFracBits[skip_offset + decision_id][0];
  if (rdCost < decisions->rdCost[decision_id])
  {
    decisions->rdCost[decision_id] = rdCost;
    decisions->absLevel[decision_id] = 0;
    decisions->prevId[decision_id] = 4 + state->m_stateId[skip_offset + decision_id];
  }
}

static INLINE void checkRdCostStart(const depquant_state* const state, int32_t lastOffset, const PQData *pqData, Decision *decisions, int
                                    decision_id)
{
  int64_t rdCost = pqData->deltaDist[decision_id] + lastOffset;
  if (pqData->absLevel[decision_id] < 4) {
    rdCost += state->m_coeffFracBits[pqData->absLevel[decision_id]];
  }
  else {
    const coeff_t value = (pqData->absLevel[decision_id] - 4) >> 1;
    rdCost += state->m_coeffFracBits[pqData->absLevel[decision_id] - (value << 1)]
              + g_goRiceBits[state->m_goRicePar][value < RICEMAX ? value : RICEMAX - 1];
  }
  if (rdCost < decisions->rdCost[decision_id]) {
    decisions->rdCost[decision_id] = rdCost;
    decisions->absLevel[decision_id] = pqData->absLevel[decision_id];
    decisions->prevId[decision_id] = -1;
  }
}



static const Decision startDec = { .rdCost = {INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2},
  .absLevel = {-1, -1, -1, -1, 0, 0, 0, 0}, .prevId = {-2, -2, -2, -2, 4, 5, 6, 7} };

static INLINE void preQuantCoeff(const quant_block * const qp, const coeff_t absCoeff, PQData* pqData, coeff_t quanCoeff)
{
  int64_t scaledOrg = (int64_t)(absCoeff) * quanCoeff;
  coeff_t  qIdx = MAX(1, (coeff_t)MIN(qp->m_maxQIdx, ((scaledOrg + qp->m_QAdd) >> qp->m_QShift)));
  int64_t scaledAdd = qIdx * qp->m_DistStepAdd - scaledOrg * qp->m_DistOrgFact;
  int index = qIdx & 3;
  pqData->deltaDist[index] = (scaledAdd * qIdx + qp->m_DistAdd) >> qp->m_DistShift;
  pqData->absLevel[index] = (++qIdx) >> 1;
  scaledAdd += qp->m_DistStepAdd;
  index = qIdx & 3;
  pqData->deltaDist[index] = (scaledAdd * qIdx + qp->m_DistAdd) >> qp->m_DistShift;
  pqData->absLevel[index] = (++qIdx) >> 1;
  scaledAdd += qp->m_DistStepAdd;
  index = qIdx & 3;
  pqData->deltaDist[index] = (scaledAdd * qIdx + qp->m_DistAdd) >> qp->m_DistShift;
  pqData->absLevel[index] = (++qIdx) >> 1;
  scaledAdd += qp->m_DistStepAdd;
  index = qIdx & 3;
  pqData->deltaDist[index] = (scaledAdd * qIdx + qp->m_DistAdd) >> qp->m_DistShift;
  pqData->absLevel[index] = (++qIdx) >> 1;
}

static void xDecide(
  all_depquant_states* const all_states,
  depquant_state* const      m_startState,
  quant_block*               qp,
  const enum ScanPosType     spt,
  const coeff_t              absCoeff,
  const int                  lastOffset,
  Decision*                  decisions,
  bool                       zeroOut,
  coeff_t                    quanCoeff,
  const int                  skip_offset,
  const int                  prev_offset)
{
  memcpy(decisions, &startDec, sizeof(Decision));

  if (zeroOut) {
    if (spt == SCAN_EOCSBB) {
      checkRdCostSkipSbbZeroOut(decisions, all_states, 0, skip_offset);
      checkRdCostSkipSbbZeroOut(decisions, all_states, 1, skip_offset);
      checkRdCostSkipSbbZeroOut(decisions, all_states, 2, skip_offset);
      checkRdCostSkipSbbZeroOut(decisions, all_states, 3, skip_offset);
    }
    return;
  }

  PQData pqData;
  preQuantCoeff(qp, absCoeff, &pqData, quanCoeff);
  uvg_dep_quant_check_rd_costs(all_states, spt, &pqData, decisions, 0, 2, prev_offset + 0);
  uvg_dep_quant_check_rd_costs(all_states, spt, &pqData, decisions, 2, 0, prev_offset + 1);
  uvg_dep_quant_check_rd_costs(all_states, spt, &pqData, decisions, 1, 3, prev_offset + 2);
  uvg_dep_quant_check_rd_costs(all_states, spt, &pqData, decisions, 3, 1, prev_offset + 3);
  if (spt == SCAN_EOCSBB) {
    checkRdCostSkipSbb(all_states, decisions, 0, skip_offset);
    checkRdCostSkipSbb(all_states, decisions, 1, skip_offset);
    checkRdCostSkipSbb(all_states, decisions, 2, skip_offset);
    checkRdCostSkipSbb(all_states, decisions, 3, skip_offset);
  }

  checkRdCostStart(m_startState, lastOffset, &pqData, decisions, 0);
  checkRdCostStart(m_startState, lastOffset, &pqData, decisions, 2);
}


static void uvg_dep_quant_decide_and_update_generic(
  rate_estimator_t*                         re,
  context_store*                          ctxs,
  struct dep_quant_scan_info const* const scan_info,
  const coeff_t                           absCoeff,
  const uint32_t                          scan_pos,
  const uint32_t                          width_in_sbb,
  const uint32_t                          height_in_sbb,
  const NbInfoSbb                         next_nb_info_ssb,
  bool                                    zeroOut,
  coeff_t                                 quantCoeff,
  const uint32_t                          effWidth,
  const uint32_t                          effHeight,
  bool                                    is_chroma)
{
  Decision* decisions = &ctxs->m_trellis[scan_pos];
  SWAP(ctxs->m_curr_state_offset, ctxs->m_prev_state_offset, int);

  enum ScanPosType spt = 0;
  if ((scan_pos & 15) == 15 && scan_pos > 16 && scan_pos < effHeight * effWidth - 1)
  {
    spt = SCAN_SOCSBB;
  }
  else if ((scan_pos & 15) == 0 && scan_pos > 0 && scan_pos < effHeight * effWidth - 16)
  {
    spt = SCAN_EOCSBB;
  }

  xDecide(&ctxs->m_allStates, &ctxs->m_startState, ctxs->m_quant, spt, absCoeff, re->m_lastBitsX[scan_info->pos_x] + re->m_lastBitsY[scan_info->pos_y], decisions, zeroOut, quantCoeff,ctxs->m_skip_state_offset, ctxs->m_prev_state_offset);
  decisions->zero_out = zeroOut;

  if (scan_pos) {
    if (!(scan_pos & 15)) {
      SWAP(ctxs->m_common_context.m_curr_sbb_ctx_offset, ctxs->m_common_context.m_prev_sbb_ctx_offset, int);
      uvg_dep_quant_update_state_eos(ctxs, scan_pos, scan_info->cg_pos, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], width_in_sbb, height_in_sbb, scan_info->next_sbb_right, scan_info->next_sbb_below, decisions, 0);
      uvg_dep_quant_update_state_eos(ctxs, scan_pos, scan_info->cg_pos, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], width_in_sbb, height_in_sbb, scan_info->next_sbb_right, scan_info->next_sbb_below, decisions, 1);
      uvg_dep_quant_update_state_eos(ctxs, scan_pos, scan_info->cg_pos, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], width_in_sbb, height_in_sbb, scan_info->next_sbb_right, scan_info->next_sbb_below, decisions, 2);
      uvg_dep_quant_update_state_eos(ctxs, scan_pos, scan_info->cg_pos, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], width_in_sbb, height_in_sbb, scan_info->next_sbb_right, scan_info->next_sbb_below, decisions, 3);
      memcpy(decisions->prevId + 4, decisions->prevId, 4 * sizeof(int32_t));
      memcpy(decisions->absLevel + 4, decisions->absLevel, 4 * sizeof(int32_t));
      memcpy(decisions->rdCost + 4, decisions->rdCost, 4 * sizeof(int64_t));
    } else if (!zeroOut) {
      uvg_dep_quant_update_state(ctxs, next_nb_info_ssb.num, scan_pos, decisions, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], next_nb_info_ssb, 4, false, 0);
      uvg_dep_quant_update_state(ctxs, next_nb_info_ssb.num, scan_pos, decisions, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], next_nb_info_ssb, 4, false, 1);
      uvg_dep_quant_update_state(ctxs, next_nb_info_ssb.num, scan_pos, decisions, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], next_nb_info_ssb, 4, false, 2);
      uvg_dep_quant_update_state(ctxs, next_nb_info_ssb.num, scan_pos, decisions, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], next_nb_info_ssb, 4, false, 3);
    }

    if (spt == SCAN_SOCSBB) {
      SWAP(ctxs->m_skip_state_offset, ctxs->m_prev_state_offset, int);
    }
  }
}


void uvg_find_first_non_zero_generic(const coeff_t* srcCoeff, const bool enableScalingLists, const context_store * const dep_quant_context, const uint32_t* const scan, const int32_t* q_coeff, int* firstTestPos, int width, int height)
{
  const int default_quant_coeff = dep_quant_context->m_quant->m_QScale;
  const int32_t thres  = dep_quant_context->m_quant->m_thresLast;
  int temp = *firstTestPos;
  for (; temp >= 0; (temp)--) {
    coeff_t thresTmp = (enableScalingLists) ? (thres / (4 * q_coeff[scan[(temp)]])) : (thres / (4 * default_quant_coeff));
    if (abs(srcCoeff[scan[(temp)]]) > thresTmp) {
      break;
    }
  }
  *firstTestPos = temp;
}

int uvg_strategy_register_depquant_generic(void* opaque, uint8_t bitdepth)
{
  bool success = true;
  
  success &= uvg_strategyselector_register(opaque, "dep_quant_decide_and_update", "generic", 0, &uvg_dep_quant_decide_and_update_generic);
  success &= uvg_strategyselector_register(opaque, "find_first_non_zero_coeff", "generic", 0, &uvg_find_first_non_zero_generic);

  return success;
}
