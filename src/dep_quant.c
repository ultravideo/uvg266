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

#include "dep_quant.h"

#include "cu.h"
#include "encoderstate.h"
#include "intra.h"
#include "rdo.h"
#include "transform.h"
#include "uvg_math.h"
#include "generic/quant-generic.h"


#define sm_numCtxSetsSig 3
#define sm_numCtxSetsGtx 2
#define sm_maxNumSigSbbCtx 2
#define sm_maxNumSigCtx    12
#define sm_maxNumGtxCtx    21
#define SCALE_BITS 15


typedef struct {
  int     m_QShift;
  int64_t m_QAdd;
  int64_t m_QScale;
  coeff_t  m_maxQIdx;
  coeff_t m_thresLast;
  coeff_t  m_thresSSbb;
  // distortion normalization
  int     m_DistShift;
  int64_t m_DistAdd;
  int64_t m_DistStepAdd;
  int64_t m_DistOrgFact;
} quant_block;

typedef struct {
  uint8_t num;
  uint8_t inPos[5];
} NbInfoSbb;

typedef struct {
  uint16_t maxDist;
  uint16_t num;
  uint16_t outPos[5];
} NbInfoOut;

typedef struct {
  uint8_t* sbbFlags;
  uint8_t* levels;
} SbbCtx;

typedef struct {
  const NbInfoOut* m_nbInfo;
  uint32_t      m_sbbFlagBits[2][2];
  SbbCtx           m_allSbbCtx[8];
  SbbCtx*          m_currSbbCtx;
  SbbCtx*          m_prevSbbCtx;
  uint8_t          m_memory[8 * (TR_MAX_WIDTH * TR_MAX_WIDTH + 1024)];
} common_context;


typedef struct 
{
  int32_t            m_lastBitsX[TR_MAX_WIDTH];
  int32_t               m_lastBitsY[TR_MAX_WIDTH];
  uint32_t        m_sigSbbFracBits[sm_maxNumSigSbbCtx][2];
  uint32_t        m_sigFracBits[sm_numCtxSetsSig][sm_maxNumSigCtx][2];
  int32_t      m_gtxFracBits[sm_maxNumGtxCtx][6];
  
} rate_estimator;


typedef struct {
  int64_t                    m_rdCost;
  uint16_t                   m_absLevelsAndCtxInit[24]; // 16x8bit for abs levels + 16x16bit for ctx init id
  int8_t                     m_numSigSbb;
  int                        m_remRegBins;
  int8_t                     m_refSbbCtxId;
  uint32_t                   m_sbbFracBits[2];
  uint32_t                   m_sigFracBits[2];
  int32_t                   m_coeffFracBits[6];
  int8_t                     m_goRicePar;
  int8_t                     m_goRiceZero;
  int8_t               m_stateId;
  const uint32_t*       m_sigFracBitsArray;
  const uint32_t*       m_gtxFracBitsArray;
  common_context*            m_commonCtx;
  
  unsigned effWidth;
  unsigned effHeight;
} depquant_state;


static void init_quant_block(
  const encoder_state_t* state,
  quant_block*           qp,
  const cu_info_t* const cur_tu,
  unsigned               log2_width,
  unsigned               log2_height,
  color_t                color,
  const bool             needsSqrt2ScaleAdjustment,
  const int              gValue)
{
  double     lambda = state->lambda;

  const int  qpDQ = state->qp + 1;
  const int  qpPer = qpDQ / 6;
  const int  qpRem = qpDQ - 6 * qpPer;
  const int  channelBitDepth = state->encoder_control->bitdepth;
  const int  maxLog2TrDynamicRange = MAX_TR_DYNAMIC_RANGE;
  const int  nomTransformShift = MAX_TR_DYNAMIC_RANGE - channelBitDepth - ((log2_width + log2_height) >> 1);
  const bool clipTransformShift = (cur_tu->tr_skip >> color) & 1 && false; // extended precision
  const int  transformShift =
    (clipTransformShift ? MAX(0, nomTransformShift) :
                          nomTransformShift) +
    (needsSqrt2ScaleAdjustment ? -1 : 0);
  // quant parameters
  qp->m_QShift = QUANT_SHIFT - 1 + qpPer + transformShift;
  qp->m_QAdd = -((3 << qp->m_QShift) >> 1);
  int invShift = IQUANT_SHIFT + 1 - qpPer - transformShift;
  qp->m_QScale = uvg_g_quant_scales[needsSqrt2ScaleAdjustment ? 1 : 0][qpRem];
  const unsigned qIdxBD = MIN(
    maxLog2TrDynamicRange + 1,
    8 * sizeof(int) + invShift - IQUANT_SHIFT - 1);
  qp->m_maxQIdx = (1 << (qIdxBD - 1)) - 4;
  qp->m_thresLast = (coeff_t)(((int64_t)(4) << qp->m_QShift));
  qp->m_thresSSbb = (coeff_t)(((int64_t)(3) << qp->m_QShift));
  // distortion calculation parameters
  const int64_t qScale = (gValue == -1) ? qp->m_QScale : gValue;
  const int     nomDShift =
    15 -
    2 * (nomTransformShift) +
    qp->m_QShift + (needsSqrt2ScaleAdjustment ? 1 : 0);
  const double qScale2 = (double)(qScale * qScale);
  const double nomDistFactor =
    (nomDShift < 0 ?
       1.0 / ((double)((int64_t)(1) << (-nomDShift)) * qScale2 * lambda) :
       (double)((int64_t)(1) << nomDShift) / (qScale2 * lambda));
  const int64_t pow2dfShift = (int64_t)(nomDistFactor * qScale2) + 1;
  assert(pow2dfShift > 0xfffffffll);
  const int dfShift = uvg_math_ceil_log2(pow2dfShift);
  qp->m_DistShift = 62 + qp->m_QShift - 2 * maxLog2TrDynamicRange - dfShift;
  qp->m_DistAdd = ((int64_t)(1) << qp->m_DistShift) >> 1;
  qp->m_DistStepAdd = (int64_t)(nomDistFactor * (double)((int64_t)(1) << (qp->m_DistShift + qp->m_QShift)) + .5);
  qp->m_DistOrgFact = (int64_t)(nomDistFactor * (double)((int64_t)(1) << (qp->m_DistShift + 1)) + .5);
}

static void reset_common_context(common_context* ctx, const rate_estimator * rate_estimator, int numSbb, int num_coeff)
{
  memset(&ctx->m_nbInfo, 0, sizeof(ctx->m_nbInfo));
  memcpy(&ctx->m_sbbFlagBits, &rate_estimator->m_sigSbbFracBits, sizeof(rate_estimator->m_sigSbbFracBits));
  const int chunkSize = numSbb + num_coeff;
  uint8_t*  nextMem   = ctx->m_memory;
  for (int k = 0; k < 8; k++, nextMem += chunkSize) {
    ctx->m_allSbbCtx[k].sbbFlags = nextMem;
    ctx->m_allSbbCtx[k].levels   = nextMem + numSbb;
  }
}

static void init_rate_esimator(rate_estimator * rate_estimator, const cabac_data_t * const ctx, color_t color)
{
  const cabac_ctx_t * base_ctx = color == COLOR_Y ? ctx->ctx.sig_coeff_group_model : (ctx->ctx.sig_coeff_group_model + 2);
  for (unsigned ctxId = 0; ctxId < sm_maxNumSigSbbCtx; ctxId++) {
    rate_estimator->m_sigSbbFracBits[ctxId][0] = CTX_ENTROPY_BITS(&base_ctx[ctxId], 0);
    rate_estimator->m_sigSbbFracBits[ctxId][1] = CTX_ENTROPY_BITS(&base_ctx[ctxId], 1);
  }
  unsigned numCtx = (color == COLOR_Y ? 12 : 8);
  for (unsigned ctxSetId = 0; ctxSetId < sm_numCtxSetsSig; ctxSetId++) {
    base_ctx = color == COLOR_Y ? ctx->ctx.cu_sig_model_luma[ctxSetId] : ctx->ctx.cu_sig_model_chroma[ctxSetId];
    for (unsigned ctxId = 0; ctxId < numCtx; ctxId++) {
      rate_estimator->m_sigFracBits[ctxSetId][ctxId][0] = CTX_ENTROPY_BITS(&base_ctx[ctxId], 0);
      rate_estimator->m_sigFracBits[ctxSetId][ctxId][1] = CTX_ENTROPY_BITS(&base_ctx[ctxId], 1);
    }
  }
  
  numCtx    = (color == COLOR_Y? 21 : 11);
  for (unsigned ctxId = 0; ctxId < numCtx; ctxId++) {
    const cabac_ctx_t * par_ctx = color == COLOR_Y ? &ctx->ctx.cu_parity_flag_model_luma[ctxId] : &ctx->ctx.cu_parity_flag_model_chroma[ctxId];
    const cabac_ctx_t * gt1_ctx = color == COLOR_Y ? &ctx->ctx.cu_gtx_flag_model_luma[0][ctxId] : &ctx->ctx.cu_gtx_flag_model_chroma[0][ctxId];
    const cabac_ctx_t * gt2_ctx = color == COLOR_Y ? &ctx->ctx.cu_gtx_flag_model_luma[1][ctxId] : &ctx->ctx.cu_gtx_flag_model_chroma[1][ctxId];

    int32_t* cb = &rate_estimator->m_gtxFracBits[ctxId];
    int32_t par0    = (1 << SCALE_BITS) + (int32_t)CTX_ENTROPY_BITS(par_ctx, 0);
    int32_t par1 = (1 << SCALE_BITS) + (int32_t)CTX_ENTROPY_BITS(par_ctx, 1);
    cb[0] = 0;
    cb[1] = CTX_ENTROPY_BITS(gt1_ctx, 0) + (1 << SCALE_BITS);
    cb[2] = CTX_ENTROPY_BITS(gt1_ctx, 1) + par0 + CTX_ENTROPY_BITS(gt2_ctx, 0);
    cb[3] = CTX_ENTROPY_BITS(gt1_ctx, 1) + par1 + CTX_ENTROPY_BITS(gt2_ctx, 0);
    cb[4] = CTX_ENTROPY_BITS(gt1_ctx, 1) + par0 + CTX_ENTROPY_BITS(gt2_ctx, 1);
    cb[5] = CTX_ENTROPY_BITS(gt1_ctx, 1) + par1 + CTX_ENTROPY_BITS(gt2_ctx, 1);
  }
}


  static void xSetLastCoeffOffset(
  const encoder_state_t* const state,
  const cu_info_t* const       cur_tu,
  const cu_loc_t* const        cu_loc,
      rate_estimator* rate_estimator,
      const bool cb_cbf,
  const color_t compID)
{
  int32_t cbfDeltaBits = 0;
  if (compID == COLOR_Y && cur_tu->type != CU_INTRA /*&& !tu.depth*/) {
    cbfDeltaBits = (int32_t)CTX_ENTROPY_BITS(&state->search_cabac.ctx.cu_qt_root_cbf_model, 1) - (int32_t)CTX_ENTROPY_BITS(&state->search_cabac.ctx.cu_qt_root_cbf_model, 0);
  } else {
    bool        prevLumaCbf           = false;
    bool        lastCbfIsInferred     = false;
    bool        useIntraSubPartitions = cur_tu->type == CU_INTRA && cur_tu->intra.isp_mode && compID == COLOR_Y;
    if (useIntraSubPartitions) {
      bool     rootCbfSoFar       = false;
      bool     isLastSubPartition = false; //TODO: isp check
      uint32_t nTus = uvg_get_isp_split_num(cu_loc->width, cu_loc->height, cur_tu->intra.isp_mode, true);
      if (isLastSubPartition) {
        //TransformUnit* tuPointer = tu.cu->firstTU;
        //for (int tuIdx = 0; tuIdx < nTus - 1; tuIdx++) {
        //  rootCbfSoFar |= TU::getCbfAtDepth(*tuPointer, COMPONENT_Y, tu.depth);
        //  tuPointer = tuPointer->next;
        //}
        if (!rootCbfSoFar) {
          lastCbfIsInferred = true;
        }
      }
      if (!lastCbfIsInferred) {
        prevLumaCbf = false;
      }
      const cabac_ctx_t * const cbf_ctx = &state->search_cabac.ctx.qt_cbf_model_luma[2 + prevLumaCbf];
      cbfDeltaBits = lastCbfIsInferred ? 0 : (int32_t)CTX_ENTROPY_BITS(cbf_ctx, 1) - (int32_t)CTX_ENTROPY_BITS(cbf_ctx, 0);
    } else {
      const cabac_ctx_t* cbf_ctx;
      switch (compID) {
        case COLOR_Y:
          cbf_ctx = &state->search_cabac.ctx.qt_cbf_model_luma[0];
          break;
        case COLOR_U:
          cbf_ctx = &state->search_cabac.ctx.qt_cbf_model_cb[0];
          break;
        case COLOR_V:
          cbf_ctx = &state->search_cabac.ctx.qt_cbf_model_cr[cb_cbf];
          break;
      }
      cbfDeltaBits = (int32_t)CTX_ENTROPY_BITS(cbf_ctx, 1) - (int32_t)CTX_ENTROPY_BITS(cbf_ctx, 0);
    }
     
  }

static const unsigned prefixCtx[] = {0, 0, 0, 3, 6, 10, 15, 21};
  uint32_t              ctxBits[14];
  for (unsigned xy = 0; xy < 2; xy++) {
    int32_t        bitOffset  = (xy ? cbfDeltaBits : 0);
    int32_t*       lastBits   = (xy ? rate_estimator->m_lastBitsY : rate_estimator->m_lastBitsX);
    const unsigned size = (xy ? (compID == COLOR_Y ? cu_loc->height : cu_loc->chroma_height) : (compID == COLOR_Y ? cu_loc->width : cu_loc->chroma_width));
    const unsigned log2Size   = uvg_math_ceil_log2(size);
    const bool     useYCtx    = (xy != 0);
    const cabac_ctx_t* const ctxSetLast = useYCtx ?
        (compID == COLOR_Y ? state->search_cabac.ctx.cu_ctx_last_y_luma : state->search_cabac.ctx.cu_ctx_last_y_chroma) :
        (compID == COLOR_Y ? state->search_cabac.ctx.cu_ctx_last_x_luma : state->search_cabac.ctx.cu_ctx_last_x_chroma);
    const unsigned lastShift = (compID == COLOR_Y ? (log2Size + 1) >> 2 : CLIP(0, 2, size >> 3));
    const unsigned lastOffset = (compID == COLOR_Y ? (prefixCtx[log2Size]) : 0);
    uint32_t sumFBits = 0;
    unsigned maxCtxId = g_group_idx[MIN(32, size) - 1];
    for (unsigned ctxId = 0; ctxId < maxCtxId; ctxId++) {
      ctxBits[ctxId] = sumFBits
        + CTX_ENTROPY_BITS(&ctxSetLast[lastOffset + (ctxId >> lastShift)], 0)
        + (ctxId > 3 ? ((ctxId - 2) >> 1) << SCALE_BITS : 0)
        + bitOffset;
      sumFBits += CTX_ENTROPY_BITS(&ctxSetLast[lastOffset + (ctxId >> lastShift)], 1);
    }
    ctxBits[maxCtxId] = sumFBits + (maxCtxId > 3 ? ((maxCtxId - 2) >> 1) << SCALE_BITS : 0) + bitOffset;
    for (unsigned pos = 0; pos < MIN(32, size); pos++) {
      lastBits[pos] = ctxBits[g_group_idx[pos]];
    }
  }
}


static void depquant_state_init(depquant_state* state, uint32_t sig_frac_bits[2], uint32_t gtx_frac_bits[6])
{
  state->m_rdCost = INT64_MAX;
  state->m_numSigSbb = 0;
  state->m_remRegBins = 4; // just large enough for last scan pos
  state->m_refSbbCtxId = -1;
  state->m_sigFracBits[0] = sig_frac_bits[0];
  state->m_sigFracBits[1] = sig_frac_bits[1];
  memcpy(state->m_coeffFracBits, gtx_frac_bits, sizeof(gtx_frac_bits));
  state->m_goRicePar = 0;
  state->m_goRiceZero = 0;
}

uint8_t uvg_dep_quant(
  const encoder_state_t* const state,
  const cu_info_t* const cur_tu,
  const cu_loc_t* const cu_loc,
  const coeff_t* srcCoeff,
  const coeff_t* coeff_out,
  const color_t compID,
  enum uvg_tree_type tree_type,
  const double lambda,
  coeff_t* absSum,
  const bool enableScalingLists)
{
  const encoder_control_t* const encoder = state->encoder_control;
  //===== reset / pre-init =====
  const int baseLevel = 4;
  
  const uint32_t  width           = compID == COLOR_Y ? cu_loc->width : cu_loc->chroma_width;
  const uint32_t  height          = compID == COLOR_Y ? cu_loc->height : cu_loc->chroma_height;
  const uint32_t  lfnstIdx = tree_type != UVG_CHROMA_T  || compID == COLOR_Y ?
                               cur_tu->lfnst_idx :
                               cur_tu->cr_lfnst_idx;

  const int       numCoeff = width * height;

  memset(coeff_out, 0x00, width * height * sizeof(coeff_t));
  *absSum                    = 0;

  const bool      is_mts   = compID == COLOR_Y && cur_tu->tr_idx > MTS_SKIP;
  const bool      is_ts    = cur_tu->tr_skip >> compID & 1;

    const uint32_t  log2_tr_width  = uvg_g_convert_to_log2[width];
  const uint32_t  log2_tr_height = uvg_g_convert_to_log2[height];
  const uint32_t* const scan     = uvg_get_scan_order_table(SCAN_GROUP_4X4,0,log2_tr_width,log2_tr_height);

  int32_t qp_scaled = uvg_get_scaled_qp(compID, state->qp, (encoder->bitdepth - 8) * 6, encoder->qp_map[0]);
  qp_scaled = is_ts ? MAX(qp_scaled, 4 + 6 * MIN_QP_PRIME_TS) : qp_scaled;
  bool needs_block_size_trafo_scale = is_ts && ((log2_tr_height + log2_tr_width) % 2 == 1);
  needs_block_size_trafo_scale |= 0; // Non log2 block size

    const int32_t scalinglist_type = (cur_tu->type == CU_INTRA ? 0 : 3) + (int8_t)compID;
  const int32_t *q_coeff = encoder->scaling_list.quant_coeff[log2_tr_width][log2_tr_height][scalinglist_type][qp_scaled % 6];
  const int32_t transform_shift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - ((log2_tr_height + log2_tr_width) >> 1) - needs_block_size_trafo_scale; //!< Represents scaling through forward transform
  const int64_t q_bits = QUANT_SHIFT + qp_scaled / 6 + (is_ts ? 0 : transform_shift );
  const int32_t add = ((state->frame->slicetype == UVG_SLICE_I) ? 171 : 85) << (q_bits - 9);

  quant_block   quant_block;
  init_quant_block(state, &quant_block, cur_tu, log2_tr_width, log2_tr_height, compID, needs_block_size_trafo_scale, -1);

  //===== scaling matrix ====
  //const int         qpDQ = cQP.Qp + 1;
  //const int         qpPer = qpDQ / 6;
  //const int         qpRem = qpDQ - 6 * qpPer;

  //TCoeff thresTmp = thres;
  bool            zeroOut         = false;
  bool            zeroOutforThres = false;
  int             effWidth = width, effHeight = height;
  if (
    (is_mts ||
     (state->encoder_control->cfg.mts && 1 /*sbt not used by block*/ &&
      height <= 32 && width <= 32)) &&
    compID == COLOR_Y) {
    effHeight = (height == 32) ? 16 : height;
    effWidth  = (width == 32) ? 16 : width;
    zeroOut   = (effHeight < height || effWidth < width);
  }
  zeroOutforThres  = zeroOut || (32 < height || 32 < width);
  //===== find first test position =====
  int firstTestPos = numCoeff - 1;
  if (
    lfnstIdx > 0 && !is_ts && width >= 4 &&
    height >= 4) {
    firstTestPos =
      ((width == 4 && height == 4) || (width == 8 && height == 8)) ? 7 : 15;
  }
  const int32_t default_quant_coeff = uvg_g_quant_scales[needs_block_size_trafo_scale][qp_scaled % 6];
  const coeff_t thres                          = 4 << q_bits;
  for (; firstTestPos >= 0; firstTestPos--) {
    coeff_t thresTmp = (enableScalingLists) ? (thres / (4 * q_coeff[firstTestPos])) :(thres / (4 * default_quant_coeff));
    if (abs(srcCoeff[firstTestPos]) > thresTmp) {
      break;
    }
  }
  if (firstTestPos < 0) {
    return 0;
  }

  //===== real init =====
  rate_estimator rate_estimator;
  init_rate_esimator(&rate_estimator, &state->search_cabac, compID);
  xSetLastCoeffOffset(state, cur_tu, cu_loc, &rate_estimator, cbf_is_set(cur_tu->cbf, COLOR_U), compID);
  common_context common_context;
  reset_common_context(&common_context, &rate_estimator, (width * height) >> 4, numCoeff);
  depquant_state all_state[12];
  depquant_state start_state;


  for (int k = 0; k < 12; k++) {
    depquant_state_init(&all_state[k], rate_estimator.m_sigFracBits[0][0], rate_estimator.m_gtxFracBits[0]);
    all_state[k].effHeight = MIN(32, effHeight);
    all_state[k].effWidth = MIN(32, effWidth);
  }
  depquant_state_init(&start_state, rate_estimator.m_sigFracBits[0][0], rate_estimator.m_gtxFracBits[0]);
  start_state.effHeight = MIN(32, effHeight);
  start_state.effWidth = MIN(32, effWidth);
  
  //===== populate trellis =====
  for (int scanIdx = firstTestPos; scanIdx >= 0; scanIdx--) {
    const ScanInfo& scanInfo = tuPars.m_scanInfo[scanIdx];
    if (enableScalingLists) {
      m_quant.initQuantBlock(
        tu,
        compID,
        cQP,
        lambda,
        quantCoeff[scanInfo.rasterPos]);
      xDecideAndUpdate(
        abs(tCoeff[scanInfo.rasterPos]),
        scanInfo,
        (zeroOut && (scanInfo.posX >= effWidth || scanInfo.posY >= effHeight)),
        quantCoeff[scanInfo.rasterPos],
        effectWidth,
        effectHeight,
        tu.cu->slice->getReverseLastSigCoeffFlag());
    } else {
      xDecideAndUpdate(
        abs(tCoeff[scanInfo.rasterPos]),
        scanInfo,
        (zeroOut && (scanInfo.posX >= effWidth || scanInfo.posY >= effHeight)),
        default_quant_coeff,
        effectWidth,
        effectHeight,
        tu.cu->slice->getReverseLastSigCoeffFlag());
    }
  }

  //===== find best path =====
  Decision decision    = {std::numeric_limits<int64_t>::max(), -1, -2};
  int64_t  minPathCost = 0;
  for (int8_t stateId = 0; stateId < 4; stateId++) {
    int64_t pathCost = m_trellis[0][stateId].rdCost;
    if (pathCost < minPathCost) {
      decision.prevId = stateId;
      minPathCost     = pathCost;
    }
  }

  //===== backward scanning =====
  int scanIdx = 0;
  for (; decision.prevId >= 0; scanIdx++) {
    decision       = m_trellis[scanIdx][decision.prevId];
    int32_t blkpos = tuPars.m_scanId2BlkPos[scanIdx].idx;
    q_coeff[blkpos] =
      (tCoeff[blkpos] < 0 ? -decision.absLevel : decision.absLevel);
    absSum += decision.absLevel;
  }
}
