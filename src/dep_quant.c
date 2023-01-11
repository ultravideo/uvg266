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
#define RICEMAX 32

static const int32_t g_goRiceBits[4][RICEMAX] = {
    { 32768,  65536,  98304, 131072, 163840, 196608, 262144, 262144, 327680, 327680, 327680, 327680, 393216, 393216, 393216, 393216, 393216, 393216, 393216, 393216, 458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752},
    { 65536,  65536,  98304,  98304, 131072, 131072, 163840, 163840, 196608, 196608, 229376, 229376, 294912, 294912, 294912, 294912, 360448, 360448, 360448, 360448, 360448, 360448, 360448, 360448, 425984, 425984, 425984, 425984, 425984, 425984, 425984, 425984},
    { 98304,  98304,  98304,  98304, 131072, 131072, 131072, 131072, 163840, 163840, 163840, 163840, 196608, 196608, 196608, 196608, 229376, 229376, 229376, 229376, 262144, 262144, 262144, 262144, 327680, 327680, 327680, 327680, 327680, 327680, 327680, 327680},
    {131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072, 163840, 163840, 163840, 163840, 163840, 163840, 163840, 163840, 196608, 196608, 196608, 196608, 196608, 196608, 196608, 196608, 229376, 229376, 229376, 229376, 229376, 229376, 229376, 229376},
};

static const int g_riceT[4] = { 32,128, 512, 2048 };
static const int g_riceShift[5] = { 0, 2, 4, 6, 8 };

static const uint32_t g_goRiceParsCoeff[32] = { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2,
                                         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3 };

enum ScanPosType { SCAN_ISCSBB = 0, SCAN_SOCSBB = 1, SCAN_EOCSBB = 2 };

typedef struct
{
  int m_QShift;
  int64_t m_QAdd;
  int64_t m_QScale;
  coeff_t m_maxQIdx;
  coeff_t m_thresLast;
  coeff_t m_thresSSbb;
  // distortion normalization
  int m_DistShift;
  int64_t m_DistAdd;
  int64_t m_DistStepAdd;
  int64_t m_DistOrgFact;
} quant_block;


typedef struct
{
  uint8_t* sbbFlags;
  uint8_t* levels;
} SbbCtx;


typedef struct
{
  coeff_t absLevel;
  int64_t deltaDist;
} PQData;

typedef struct
{
  int64_t rdCost;
  coeff_t absLevel;
  int prevId;
} Decision;


typedef struct
{
  const NbInfoOut* m_nbInfo;
  uint32_t m_sbbFlagBits[2][2];
  SbbCtx m_allSbbCtx[8];
  SbbCtx* m_currSbbCtx;
  SbbCtx* m_prevSbbCtx;
  uint8_t m_memory[8 * (TR_MAX_WIDTH * TR_MAX_WIDTH + 1024)];
} common_context;


typedef struct
{
  int32_t m_lastBitsX[TR_MAX_WIDTH];
  int32_t m_lastBitsY[TR_MAX_WIDTH];
  uint32_t m_sigSbbFracBits[sm_maxNumSigSbbCtx][2];
  uint32_t m_sigFracBits[sm_numCtxSetsSig][sm_maxNumSigCtx][2];
  int32_t m_gtxFracBits[sm_maxNumGtxCtx][6];
} rate_estimator;


typedef struct
{
  int64_t m_rdCost;
  uint16_t m_absLevelsAndCtxInit[24]; // 16x8bit for abs levels + 16x16bit for ctx init id
  int8_t m_numSigSbb;
  int m_remRegBins;
  int8_t m_refSbbCtxId;
  uint32_t m_sbbFracBits[2];
  uint32_t m_sigFracBits[2];
  int32_t m_coeffFracBits[6];
  int8_t m_goRicePar;
  int8_t m_goRiceZero;
  int8_t m_stateId;
  const uint32_t* m_sigFracBitsArray[2];
  const uint32_t* m_gtxFracBitsArray[6];
  struct common_context* m_commonCtx;

  unsigned effWidth;
  unsigned effHeight;
} depquant_state;

typedef struct
{
    common_context   m_common_context;
    depquant_state       m_allStates[12];
    depquant_state* m_currStates;
    depquant_state* m_prevStates;
    depquant_state* m_skipStates;
    depquant_state       m_startState;
    quant_block   m_quant;
    Decision    m_trellis[TR_MAX_WIDTH * TR_MAX_WIDTH][8];
} context_store;


int uvg_init_nb_info(encoder_control_t * encoder) {
  memset(encoder->m_scanId2NbInfoSbbArray, 0, sizeof(encoder->m_scanId2NbInfoSbbArray));
  memset(encoder->m_scanId2NbInfoOutArray, 0, sizeof(encoder->m_scanId2NbInfoOutArray));
  for (int hd = 0; hd <= 7; hd++)
  {

    uint32_t raster2id[64 * 64] = {0};

    for (int vd = 0; vd <= 7; vd++)
    {
      if ((hd == 0 && vd <= 1) || (hd <= 1 && vd == 0))
      {
        continue;
      }
      const uint32_t      blockWidth = (1 << hd);
      const uint32_t      blockHeight = (1 << vd);
      const uint32_t      log2CGWidth = g_log2_sbb_size[hd][vd][0];
      const uint32_t      log2CGHeight = g_log2_sbb_size[hd][vd][1];
      const uint32_t      groupWidth = 1 << log2CGWidth;
      const uint32_t      groupHeight = 1 << log2CGHeight;
      const uint32_t      groupSize = groupWidth * groupHeight;
      const int           scanType = SCAN_DIAG;
      const uint32_t      blkWidthIdx = hd;
      const uint32_t      blkHeightIdx = vd;
      const uint32_t* scanId2RP = uvg_get_scan_order_table(SCAN_GROUP_4X4, scanType, blkWidthIdx, blkHeightIdx);
      NbInfoSbb** sId2NbSbb = &encoder->m_scanId2NbInfoSbbArray[hd][vd];
      NbInfoOut** sId2NbOut = &encoder->m_scanId2NbInfoOutArray[hd][vd];
      // consider only non-zero-out region
      const uint32_t      blkWidthNZOut = MIN(32, blockWidth);
      const uint32_t      blkHeightNZOut = MIN(32, blockHeight);
      const uint32_t      totalValues = blkWidthNZOut * blkHeightNZOut;

      *sId2NbSbb = MALLOC(NbInfoSbb, totalValues);
      if (*sId2NbSbb == NULL) {
        return 0;
      }
      *sId2NbOut = MALLOC(NbInfoOut, totalValues);
      if (*sId2NbOut == NULL) {
        return 0;
      }

      for (uint32_t scanId = 0; scanId < totalValues; scanId++)
      {
        raster2id[scanId2RP[scanId]] = scanId;
      }

      for (unsigned scanId = 0; scanId < totalValues; scanId++)
      {
        const int rpos = scanId2RP[scanId];
        uint32_t  pos_y = rpos >> hd;
        uint32_t  pos_x = rpos - (pos_y << hd); // TODO: height
        {
          //===== inside subband neighbours =====
          NbInfoSbb *nbSbb = &(*sId2NbSbb)[scanId];
          const int      begSbb = scanId - (scanId & (groupSize - 1)); // first pos in current subblock
          int            cpos[5];

          cpos[0] = (pos_x + 1 < blkWidthNZOut ? (raster2id[rpos + 1] < groupSize + begSbb ? raster2id[rpos + 1] - begSbb : 0) : 0);
          cpos[1] = (pos_x + 2 < blkWidthNZOut ? (raster2id[rpos + 2] < groupSize + begSbb ? raster2id[rpos + 2] - begSbb : 0) : 0);
          cpos[2] = (pos_x + 1 < blkWidthNZOut && pos_y + 1 < blkHeightNZOut ? (raster2id[rpos + 1 + blockWidth] < groupSize + begSbb ? raster2id[rpos + 1 + blockWidth] - begSbb : 0) : 0);
          cpos[3] = (pos_y + 1 < blkHeightNZOut ? (raster2id[rpos + blockWidth] < groupSize + begSbb ? raster2id[rpos + blockWidth] - begSbb : 0) : 0);
          cpos[4] = (pos_y + 2 < blkHeightNZOut ? (raster2id[rpos + 2 * blockWidth] < groupSize + begSbb ? raster2id[rpos + 2 * blockWidth] - begSbb : 0) : 0);

          for (nbSbb->num = 0; true; )
          {
            int nk = -1;
            for (int k = 0; k < 5; k++)
            {
              if (cpos[k] != 0 && (nk < 0 || cpos[k] < cpos[nk]))
              {
                nk = k;
              }
            }
            if (nk < 0)
            {
              break;
            }
            nbSbb->inPos[nbSbb->num++] = (uint8_t)(cpos[nk]);
            cpos[nk] = 0;
          }
          for (int k = nbSbb->num; k < 5; k++)
          {
            nbSbb->inPos[k] = 0;
          }
        }
        {
          //===== outside subband neighbours =====
          NbInfoOut *nbOut = &(*sId2NbOut)[scanId];
          const int      begSbb = scanId - (scanId & (groupSize - 1)); // first pos in current subblock
          int            cpos[5];

          cpos[0] = (pos_x + 1 < blkWidthNZOut ? (raster2id[rpos + 1] >= groupSize + begSbb ? raster2id[rpos + 1] : 0) : 0);
          cpos[1] = (pos_x + 2 < blkWidthNZOut ? (raster2id[rpos + 2] >= groupSize + begSbb ? raster2id[rpos + 2] : 0) : 0);
          cpos[2] = (pos_x + 1 < blkWidthNZOut && pos_y + 1 < blkHeightNZOut ? (raster2id[rpos + 1 + blockWidth] >= groupSize + begSbb ? raster2id[rpos + 1 + blockWidth] : 0) : 0);
          cpos[3] = (pos_y + 1 < blkHeightNZOut ? (raster2id[rpos + blockWidth] >= groupSize + begSbb ? raster2id[rpos + blockWidth] : 0) : 0);
          cpos[4] = (pos_y + 2 < blkHeightNZOut ? (raster2id[rpos + 2 * blockWidth] >= groupSize + begSbb ? raster2id[rpos + 2 * blockWidth] : 0) : 0);

          for (nbOut->num = 0; true; )
          {
            int nk = -1;
            for (int k = 0; k < 5; k++)
            {
              if (cpos[k] != 0 && (nk < 0 || cpos[k] < cpos[nk]))
              {
                nk = k;
              }
            }
            if (nk < 0)
            {
              break;
            }
            nbOut->outPos[nbOut->num++] = (uint16_t)(cpos[nk]);
            cpos[nk] = 0;
          }
          for (int k = nbOut->num; k < 5; k++)
          {
            nbOut->outPos[k] = 0;
          }
          nbOut->maxDist = (scanId == 0 ? 0 : (*sId2NbOut)[scanId - 1].maxDist);
          for (int k = 0; k < nbOut->num; k++)
          {
            if (nbOut->outPos[k] > nbOut->maxDist)
            {
              nbOut->maxDist = nbOut->outPos[k];
            }
          }
        }
      }

      // make it relative
      for (unsigned scanId = 0; scanId < totalValues; scanId++)
      {
        NbInfoOut *nbOut = &(*sId2NbOut)[scanId];
        const int  begSbb = scanId - (scanId & (groupSize - 1)); // first pos in current subblock
        for (int k = 0; k < nbOut->num; k++)
        {
          nbOut->outPos[k] -= begSbb;
        }
        nbOut->maxDist -= scanId;
      }
    }
  }
  return 1;
}

void uvg_dealloc_nb_info(encoder_control_t* encoder) {

  for (int hd = 0; hd <= 7; hd++) {
    for (int vd = 0; vd <= 7; vd++)
    {
      if ((hd == 0 && vd <= 1) || (hd <= 1 && vd == 0))
      {
        continue;
      }
      if(encoder->m_scanId2NbInfoOutArray[hd][vd]) FREE_POINTER(encoder->m_scanId2NbInfoOutArray[hd][vd]);
      if(encoder->m_scanId2NbInfoOutArray[hd][vd]) FREE_POINTER(encoder->m_scanId2NbInfoSbbArray[hd][vd]);
    }
  }
}


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
  //memset(&ctx->m_nbInfo, 0, sizeof(ctx->m_nbInfo));
  memcpy(&ctx->m_sbbFlagBits, &rate_estimator->m_sigSbbFracBits, sizeof(rate_estimator->m_sigSbbFracBits));
  const int chunkSize = numSbb + num_coeff;
  uint8_t*  nextMem   = ctx->m_memory;
  for (int k = 0; k < 8; k++, nextMem += chunkSize) {
    ctx->m_allSbbCtx[k].sbbFlags = nextMem;
    ctx->m_allSbbCtx[k].levels   = nextMem + numSbb;
  }
  ctx->m_currSbbCtx = &ctx->m_allSbbCtx[0];
  ctx->m_prevSbbCtx = &ctx->m_allSbbCtx[4];
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

static INLINE void checkRdCostSkipSbbZeroOut(Decision *decision, const depquant_state * const state) 
{
    int64_t rdCost = state->m_rdCost + state->m_sbbFracBits[0];
    decision->rdCost = rdCost;
    decision->absLevel = 0;
    decision->prevId = 4 + state->m_stateId;
}

static void checkRdCosts(const depquant_state * const state, const enum ScanPosType spt, const PQData *pqDataA, const PQData *pqDataB, Decision *decisionA, Decision *decisionB)
{
    const int32_t* goRiceTab = g_goRiceBits[state->m_goRicePar];
    int64_t         rdCostA = state->m_rdCost + pqDataA->deltaDist;
    int64_t         rdCostB = state->m_rdCost + pqDataB->deltaDist;
    int64_t         rdCostZ = state->m_rdCost;
    if (state->m_remRegBins >= 4)
    {
        if (pqDataA->absLevel < 4)
        {
            rdCostA += state->m_coeffFracBits[pqDataA->absLevel];
        }
        else
        {
            const coeff_t value = (pqDataA->absLevel - 4) >> 1;
            rdCostA +=
                state->m_coeffFracBits[pqDataA->absLevel - (value << 1)] + goRiceTab[value < RICEMAX ? value : RICEMAX - 1];
        }
        if (pqDataB->absLevel < 4)
        {
            rdCostB += state->m_coeffFracBits[pqDataB->absLevel];
        }
        else
        {
            const coeff_t value = (pqDataB->absLevel - 4) >> 1;
            rdCostB +=
                state->m_coeffFracBits[pqDataB->absLevel - (value << 1)] + goRiceTab[value < RICEMAX ? value : RICEMAX - 1];
        }
        if (spt == SCAN_ISCSBB)
        {
            rdCostA += state->m_sigFracBits[1];
            rdCostB += state->m_sigFracBits[1];
            rdCostZ += state->m_sigFracBits[0];
        }
        else if (spt == SCAN_SOCSBB)
        {
            rdCostA += state->m_sbbFracBits[1] + state->m_sigFracBits[1];
            rdCostB += state->m_sbbFracBits[1] + state->m_sigFracBits[1];
            rdCostZ += state->m_sbbFracBits[1] + state->m_sigFracBits[0];
        }
        else if (state->m_numSigSbb)
        {
            rdCostA += state->m_sigFracBits[1];
            rdCostB += state->m_sigFracBits[1];
            rdCostZ += state->m_sigFracBits[0];
        }
        else
        {
            rdCostZ = decisionA->rdCost;
        }
    }
    else
    {
        rdCostA +=
            (1 << SCALE_BITS)
            + goRiceTab[pqDataA->absLevel <= state->m_goRiceZero ? pqDataA->absLevel - 1
            : (pqDataA->absLevel < RICEMAX ? pqDataA->absLevel : RICEMAX - 1)];
        rdCostB +=
            (1 << SCALE_BITS)
            + goRiceTab[pqDataB->absLevel <= state->m_goRiceZero ? pqDataB->absLevel - 1
            : (pqDataB->absLevel < RICEMAX ? pqDataB->absLevel : RICEMAX - 1)];
        rdCostZ += goRiceTab[state->m_goRiceZero];
    }
    if (rdCostA < decisionA->rdCost)
    {
        decisionA->rdCost = rdCostA;
        decisionA->absLevel = pqDataA->absLevel;
        decisionA->prevId = state->m_stateId;
    }
    if (rdCostZ < decisionA->rdCost)
    {
        decisionA->rdCost = rdCostZ;
        decisionA->absLevel = 0;
        decisionA->prevId = state->m_stateId;
    }
    if (rdCostB < decisionB->rdCost)
    {
        decisionB->rdCost = rdCostB;
        decisionB->absLevel = pqDataB->absLevel;
        decisionB->prevId = state->m_stateId;
    }
}

static INLINE void checkRdCostSkipSbb(const depquant_state* const state, Decision *decision)
{
    int64_t rdCost = state->m_rdCost + state->m_sbbFracBits[0];
    if (rdCost < decision->rdCost)
    {
        decision->rdCost = rdCost;
        decision->absLevel = 0;
        decision->prevId = 4 + state->m_stateId;
    }
}

static INLINE void checkRdCostStart(const depquant_state* const state, int32_t lastOffset, const PQData *pqData, Decision *decision)
{
    int64_t rdCost = pqData->deltaDist + lastOffset;
    if (pqData->absLevel < 4)
    {
        rdCost += state->m_coeffFracBits[pqData->absLevel];
    }
    else
    {
        const coeff_t value = (pqData->absLevel - 4) >> 1;
        rdCost += state->m_coeffFracBits[pqData->absLevel - (value << 1)] + g_goRiceBits[state->m_goRicePar][value < RICEMAX ? value : RICEMAX - 1];
    }
    if (rdCost < decision->rdCost)
    {
        decision->rdCost = rdCost;
        decision->absLevel = pqData->absLevel;
        decision->prevId = -1;
    }
}


static INLINE void preQuantCoeff(const quant_block * const qp, const coeff_t absCoeff, PQData* pqData, coeff_t quanCoeff)
{
    int64_t scaledOrg = (int64_t)(absCoeff) * quanCoeff;
    coeff_t  qIdx = MAX(1, MIN(qp->m_maxQIdx, (coeff_t)((scaledOrg + qp->m_QAdd) >> qp->m_QShift)));
    int64_t scaledAdd = qIdx * qp->m_DistStepAdd - scaledOrg * qp->m_DistOrgFact;
    PQData *pq_a = &pqData[qIdx & 3];
    pq_a->deltaDist = (scaledAdd * qIdx + qp->m_DistAdd) >> qp->m_DistShift;
    pq_a->absLevel = (++qIdx) >> 1;
    scaledAdd += qp->m_DistStepAdd;
    PQData *pq_b = &pqData[qIdx & 3];
    pq_b->deltaDist = (scaledAdd * qIdx + qp->m_DistAdd) >> qp->m_DistShift;
    pq_b->absLevel = (++qIdx) >> 1;
    scaledAdd += qp->m_DistStepAdd;
    PQData *pq_c = &pqData[qIdx & 3];
    pq_c->deltaDist = (scaledAdd * qIdx + qp->m_DistAdd) >> qp->m_DistShift;
    pq_c->absLevel = (++qIdx) >> 1;
    scaledAdd += qp->m_DistStepAdd;
    PQData *pq_d = &pqData[qIdx & 3];
    pq_d->deltaDist = (scaledAdd * qIdx + qp->m_DistAdd) >> qp->m_DistShift;
    pq_d->absLevel = (++qIdx) >> 1;
}


#define DINIT(l,p) {INT64_MAX>>2,(l),(p)}
static const Decision startDec[8] = { DINIT(-1,-2),DINIT(-1,-2),DINIT(-1,-2),DINIT(-1,-2),DINIT(0,4),DINIT(0,5),DINIT(0,6),DINIT(0,7) };
#undef  DINIT


static void xDecide(
    depquant_state* const m_skipStates,
    depquant_state* const m_prevStates,
    depquant_state* const m_startState,
    quant_block *qp,
    const enum ScanPosType spt,
    const coeff_t absCoeff,
    const int lastOffset,
    Decision* decisions,
    bool zeroOut,
    coeff_t quanCoeff)
{
    memcpy(decisions, startDec, 8 * sizeof(Decision));

    if (zeroOut)
    {
        if (spt == SCAN_EOCSBB)
        {
            checkRdCostSkipSbbZeroOut(&decisions[0], &m_skipStates[0]);
            checkRdCostSkipSbbZeroOut(&decisions[1], &m_skipStates[1]);
            checkRdCostSkipSbbZeroOut(&decisions[2], &m_skipStates[2]);
            checkRdCostSkipSbbZeroOut(&decisions[3], &m_skipStates[3]);
        }
        return;
    }

    PQData  pqData[4];
    preQuantCoeff(qp, absCoeff, pqData, quanCoeff);
    checkRdCosts(&m_prevStates[0], spt, &pqData[0], &pqData[2], &decisions[0], &decisions[2]);
    checkRdCosts(&m_prevStates[1], spt, &pqData[0], &pqData[2], &decisions[2], &decisions[0]);
    checkRdCosts(&m_prevStates[2], spt, &pqData[3], &pqData[1], &decisions[1], &decisions[3]);
    checkRdCosts(&m_prevStates[3], spt, &pqData[3], &pqData[1], &decisions[3], &decisions[1]);
    if (spt == SCAN_EOCSBB)
    {
        checkRdCostSkipSbb(&m_skipStates[0], &decisions[0]);
        checkRdCostSkipSbb(&m_skipStates[1], &decisions[1]);
        checkRdCostSkipSbb(&m_skipStates[2], &decisions[2]);
        checkRdCostSkipSbb(&m_skipStates[3], &decisions[3]);
    }

    checkRdCostStart(m_startState, lastOffset, &pqData[0], &decisions[0]);
    checkRdCostStart(m_startState, lastOffset, &pqData[2], &decisions[2]);
}


unsigned templateAbsCompare(coeff_t sum)
{
    int rangeIdx = 0;
    if (sum < g_riceT[0])
    {
        rangeIdx = 0;
    }
    else if (sum < g_riceT[1])
    {
        rangeIdx = 1;
    }
    else if (sum < g_riceT[2])
    {
        rangeIdx = 2;
    }
    else if (sum < g_riceT[3])
    {
        rangeIdx = 3;
    }
    else
    {
        rangeIdx = 4;
    }
    return g_riceShift[rangeIdx];
}

static INLINE void update_common_context(
  common_context * cc,
  const uint32_t scan_pos,
  const uint32_t width_in_sbb,
  const uint32_t height_in_sbb,
  const int sigNSbb,
  const depquant_state* prevState,
  depquant_state *currState)
{
  const uint32_t numSbb = width_in_sbb * height_in_sbb;
  uint8_t* sbbFlags = cc->m_currSbbCtx[currState->m_stateId].sbbFlags;
  uint8_t* levels = cc->m_currSbbCtx[currState->m_stateId].levels;
  size_t setCpSize = cc->m_nbInfo[scan_pos - 1].maxDist * sizeof(uint8_t);
  if (prevState && prevState->m_refSbbCtxId >= 0) {
    memcpy(sbbFlags, cc->m_prevSbbCtx[prevState->m_refSbbCtxId].sbbFlags, numSbb * sizeof(uint8_t));
    memcpy(levels + scan_pos, cc->m_prevSbbCtx[prevState->m_refSbbCtxId].levels + scan_pos, setCpSize);
  }
  else {
    memset(sbbFlags, 0, numSbb * sizeof(uint8_t));
    memset(levels + scan_pos, 0, setCpSize);
  }
  sbbFlags[scan_pos >> 4] = !!currState->m_numSigSbb;
  memcpy(levels + scan_pos, currState->m_absLevelsAndCtxInit, 16 * sizeof(uint8_t));

  currState->m_numSigSbb = 0;
  if (prevState) {
    currState->m_remRegBins = prevState->m_remRegBins;
  }
  else {
    int ctxBinSampleRatio = 28;
    // (scanInfo.chType == COLOR_Y) ? MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_LUMA : MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_CHROMA;
    currState->m_remRegBins = (currState->effWidth * currState->effHeight * ctxBinSampleRatio) / 16;
  }
  currState->m_goRicePar = 0;
  currState->m_refSbbCtxId = currState->m_stateId;
  currState->m_sbbFracBits[0] = cc->m_sbbFlagBits[sigNSbb][0];
  currState->m_sbbFracBits[1] = cc->m_sbbFlagBits[sigNSbb][1];

  uint16_t templateCtxInit[16];
  const int scanBeg = scan_pos - 16;
  const NbInfoOut* nbOut = cc->m_nbInfo + scanBeg;
  const uint8_t* absLevels = levels + scanBeg;
  for (int id = 0; id < 16; id++, nbOut++) {
    if (nbOut->num) {
      coeff_t sumAbs = 0, sumAbs1 = 0, sumNum = 0;
#define UPDATE(k) {coeff_t t=absLevels[nbOut->outPos[k]]; sumAbs+=t; sumAbs1+=MIN(4+(t&1),t); sumNum+=!!t; }
      UPDATE(0);
      if (nbOut->num > 1) {
        UPDATE(1);
        if (nbOut->num > 2) {
          UPDATE(2);
          if (nbOut->num > 3) {
            UPDATE(3);
            if (nbOut->num > 4) {
              UPDATE(4);
            }
          }
        }
      }
#undef UPDATE
      templateCtxInit[id] = (uint16_t)(sumNum) + ((uint16_t)(sumAbs1) << 3) + ((uint16_t)MIN(127, sumAbs) << 8);
    }
    else {
      templateCtxInit[id] = 0;
    }
  }
  memset(currState->m_absLevelsAndCtxInit, 0, 16 * sizeof(uint8_t));
  memcpy(currState->m_absLevelsAndCtxInit + 8, templateCtxInit, 16 * sizeof(uint16_t));
}


static INLINE void updateStateEOS(
  depquant_state * state, 
  const uint32_t  scan_pos,
  const uint32_t sigCtxOffsetNext,
  const uint32_t gtxCtxOffsetNext,
  const uint32_t width_in_sbb,
  const uint32_t height_in_sbb,
  const uint32_t sigNSbb,
  const depquant_state* prevStates,
  const depquant_state* skipStates,
  const Decision *decision)
{
    state->m_rdCost = decision->rdCost;
    if (decision->prevId > -2)
    {
      const depquant_state* prvState = 0;
      if (decision->prevId >= 4)
      {
          prvState = skipStates + (decision->prevId - 4);
          state->m_numSigSbb = 0;
          memset(state->m_absLevelsAndCtxInit, 0, 16 * sizeof(uint8_t));
      }
      else if (decision->prevId >= 0)
      {
          prvState = prevStates + decision->prevId;
          state->m_numSigSbb = prvState->m_numSigSbb + !!decision->absLevel;
          memcpy(state->m_absLevelsAndCtxInit, prvState->m_absLevelsAndCtxInit, 16 * sizeof(uint8_t));
      }
      else
      {
          state->m_numSigSbb = 1;
          memset(state->m_absLevelsAndCtxInit, 0, 16 * sizeof(uint8_t));
      }
      uint8_t* temp = (uint8_t*)(state->m_absLevelsAndCtxInit[scan_pos & 15]);
      *temp = (uint8_t)MIN(255, decision->absLevel);

      update_common_context(state->m_commonCtx, scan_pos, width_in_sbb, height_in_sbb, sigNSbb, prvState, state);
      
      coeff_t  tinit = state->m_absLevelsAndCtxInit[8 + ((scan_pos - 1) & 15)];
      coeff_t  sumNum = tinit & 7;
      coeff_t  sumAbs1 = (tinit >> 3) & 31;
      coeff_t  sumGt1 = sumAbs1 - sumNum;
      state->m_sigFracBits[0] = state->m_sigFracBitsArray[sigCtxOffsetNext + MIN((sumAbs1 + 1) >> 1, 3)][0];
      state->m_sigFracBits[1] = state->m_sigFracBitsArray[sigCtxOffsetNext + MIN((sumAbs1 + 1) >> 1, 3)][1];
      
      memcpy(state->m_coeffFracBits, state->m_gtxFracBitsArray[gtxCtxOffsetNext + (sumGt1 < 4 ? sumGt1 : 4)], sizeof(state->m_coeffFracBits));
    }
}

static INLINE void updateState(
  depquant_state* state, 
  int numIPos, const uint32_t scan_pos,
  const depquant_state* prevStates,
  const Decision* decision,
  const uint32_t sigCtxOffsetNext,
  const uint32_t gtxCtxOffsetNext,
  const NbInfoSbb next_nb_info_ssb,
  const int baseLevel,
  const bool extRiceFlag) {
    state->m_rdCost = decision->rdCost;
    if (decision->prevId > -2)
    {
        if (decision->prevId >= 0)
        {
            const depquant_state* prvState = prevStates + decision->prevId;
            state->m_numSigSbb = prvState->m_numSigSbb + !!decision->absLevel;
            state->m_refSbbCtxId = prvState->m_refSbbCtxId;
            state->m_sbbFracBits[0] = prvState->m_sbbFracBits[0];
            state->m_sbbFracBits[1] = prvState->m_sbbFracBits[1];
            state->m_remRegBins = prvState->m_remRegBins - 1;
            state->m_goRicePar = prvState->m_goRicePar;
            if (state->m_remRegBins >= 4)
            {
                state->m_remRegBins -= (decision->absLevel < 2 ? (unsigned)decision->absLevel : 3);
            }
            memcpy(state->m_absLevelsAndCtxInit, prvState->m_absLevelsAndCtxInit, 48 * sizeof(uint8_t));
        }
        else
        {
            state->m_numSigSbb = 1;
            state->m_refSbbCtxId = -1;
            int ctxBinSampleRatio = 28; //(scanInfo.chType == CHANNEL_TYPE_LUMA) ? MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_LUMA : MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_CHROMA;
            state->m_remRegBins = (state->effWidth * state->effHeight * ctxBinSampleRatio) / 16 - (decision->absLevel < 2 ? (unsigned)decision->absLevel : 3);
            memset(state->m_absLevelsAndCtxInit, 0, 48 * sizeof(uint8_t));
        }

        uint8_t* levels = (uint8_t*)(state->m_absLevelsAndCtxInit);
        levels[scan_pos & 15] = (uint8_t)MIN(255, decision->absLevel);

        if (state->m_remRegBins >= 4)
        {
            coeff_t  tinit = state->m_absLevelsAndCtxInit[8 + ((scan_pos - 1) & 15)];
            coeff_t  sumAbs1 = (tinit >> 3) & 31;
            coeff_t sumNum = tinit & 7;
#define UPDATE(k) {coeff_t t=levels[next_nb_info_ssb.inPos[k]]; sumAbs1+=MIN(4+(t&1),t); sumNum+=!!t; }
            if (numIPos == 1)
            {
                UPDATE(0);
            }
            else if (numIPos == 2)
            {
                UPDATE(0);
                UPDATE(1);
            }
            else if (numIPos == 3)
            {
                UPDATE(0);
                UPDATE(1);
                UPDATE(2);
            }
            else if (numIPos == 4)
            {
                UPDATE(0);
                UPDATE(1);
                UPDATE(2);
                UPDATE(3);
            }
            else if (numIPos == 5)
            {
                UPDATE(0);
                UPDATE(1);
                UPDATE(2);
                UPDATE(3);
                UPDATE(4);
            }
#undef UPDATE
            coeff_t sumGt1 = sumAbs1 - sumNum;
            state->m_sigFracBits[0] = state->m_sigFracBitsArray[sigCtxOffsetNext + MIN((sumAbs1 + 1) >> 1, 3)][0];
            state->m_sigFracBits[1] = state->m_sigFracBitsArray[sigCtxOffsetNext + MIN((sumAbs1 + 1) >> 1, 3)][1];
            memcpy(state->m_coeffFracBits, state->m_gtxFracBitsArray[gtxCtxOffsetNext + (sumGt1 < 4 ? sumGt1 : 4)], sizeof(state->m_coeffFracBits));
            

            coeff_t  sumAbs = state->m_absLevelsAndCtxInit[8 + ((scan_pos - 1) & 15)] >> 8;
#define UPDATE(k) {coeff_t t=levels[next_nb_info_ssb.inPos[k]]; sumAbs+=t; }
            if (numIPos == 1)
            {
                UPDATE(0);
            }
            else if (numIPos == 2)
            {
                UPDATE(0);
                UPDATE(1);
            }
            else if (numIPos == 3)
            {
                UPDATE(0);
                UPDATE(1);
                UPDATE(2);
            }
            else if (numIPos == 4)
            {
                UPDATE(0);
                UPDATE(1);
                UPDATE(2);
                UPDATE(3);
            }
            else if (numIPos == 5)
            {
                UPDATE(0);
                UPDATE(1);
                UPDATE(2);
                UPDATE(3);
                UPDATE(4);
            }
#undef UPDATE
            if (extRiceFlag)
            {
                unsigned currentShift = templateAbsCompare(sumAbs);
                sumAbs = sumAbs >> currentShift;
                int sumAll = MAX(MIN(31, (int)sumAbs - (int)baseLevel), 0);
                state->m_goRicePar = g_goRiceParsCoeff[sumAll];
                state->m_goRicePar += currentShift;
            }
            else
            {
                int sumAll = MAX(MIN(31, (int)sumAbs - 4 * 5), 0);
                state->m_goRicePar = g_goRiceParsCoeff[sumAll];
            }
        }
        else
        {
            coeff_t  sumAbs = state->m_absLevelsAndCtxInit[8 + ((scan_pos - 1) & 15)] >> 8;
#define UPDATE(k) {coeff_t t=levels[next_nb_info_ssb.inPos[k]]; sumAbs+=t; }
            if (numIPos == 1)
            {
                UPDATE(0);
            }
            else if (numIPos == 2)
            {
                UPDATE(0);
                UPDATE(1);
            }
            else if (numIPos == 3)
            {
                UPDATE(0);
                UPDATE(1);
                UPDATE(2);
            }
            else if (numIPos == 4)
            {
                UPDATE(0);
                UPDATE(1);
                UPDATE(2);
                UPDATE(3);
            }
            else if (numIPos == 5)
            {
                UPDATE(0);
                UPDATE(1);
                UPDATE(2);
                UPDATE(3);
                UPDATE(4);
            }
#undef UPDATE
            if (extRiceFlag)
            {
                unsigned currentShift = templateAbsCompare(sumAbs);
                sumAbs = sumAbs >> currentShift;
                sumAbs = MIN(31, sumAbs);
                state->m_goRicePar = g_goRiceParsCoeff[sumAbs];
                state->m_goRicePar += currentShift;
            }
            else
            {
                sumAbs = MIN(31, sumAbs);
                state->m_goRicePar = g_goRiceParsCoeff[sumAbs];
            }
            state->m_goRiceZero = (state->m_stateId < 2 ? 1 : 2) << state->m_goRicePar; 
        }
    }
}

static void xDecideAndUpdate(
  rate_estimator* re,
  context_store* ctxs,
  const coeff_t absCoeff,
  const uint32_t scan_pos,
  const uint32_t pos_x,
  const uint32_t pos_y,
  const uint32_t sigCtxOffsetNext,
  const uint32_t gtxCtxOffsetNext,
  const uint32_t width_in_sbb,
  const uint32_t height_in_sbb,
  const uint32_t sigNSbb,
  const NbInfoSbb next_nb_info_ssb,
  bool zeroOut,
  coeff_t quantCoeff,
  int effWidth,
  int effHeight)
{
  Decision* decisions = ctxs->m_trellis[scan_pos];
  SWAP(ctxs->m_currStates, ctxs->m_prevStates, depquant_state*);

  enum ScanPosType spt = 0;
  if ((scan_pos & 15) == 15 && scan_pos > 16 && scan_pos < effHeight * effWidth - 1)
  {
    spt = SCAN_SOCSBB;
  }
  else if ((scan_pos & 15) == 0 && scan_pos > 0 && scan_pos < effHeight * effWidth - 16)
  {
    spt = SCAN_EOCSBB;
  }

  xDecide(ctxs->m_skipStates, ctxs->m_prevStates, &ctxs->m_startState, &ctxs->m_quant, spt, absCoeff, re->m_lastBitsX[pos_x] + re->m_lastBitsY[pos_y], decisions, zeroOut, quantCoeff);

  if (scan_pos) {
    if (!(scan_pos & 15)) {
      SWAP(ctxs->m_common_context.m_currSbbCtx, ctxs->m_common_context.m_prevSbbCtx, SbbCtx*);
      updateStateEOS(&ctxs->m_currStates[0], scan_pos, sigCtxOffsetNext, gtxCtxOffsetNext, width_in_sbb, height_in_sbb, sigNSbb, ctxs->m_prevStates, ctxs->m_skipStates, &decisions[0]);
      updateStateEOS(&ctxs->m_currStates[1], scan_pos, sigCtxOffsetNext, gtxCtxOffsetNext, width_in_sbb, height_in_sbb, sigNSbb, ctxs->m_prevStates, ctxs->m_skipStates, &decisions[1]);
      updateStateEOS(&ctxs->m_currStates[2], scan_pos, sigCtxOffsetNext, gtxCtxOffsetNext, width_in_sbb, height_in_sbb, sigNSbb, ctxs->m_prevStates, ctxs->m_skipStates, &decisions[2]);
      updateStateEOS(&ctxs->m_currStates[3], scan_pos, sigCtxOffsetNext, gtxCtxOffsetNext, width_in_sbb, height_in_sbb, sigNSbb, ctxs->m_prevStates, ctxs->m_skipStates, &decisions[3]);
      memcpy(decisions + 4, decisions, 4 * sizeof(Decision));
    } else if (!zeroOut) {

      updateState(&ctxs->m_currStates[0], next_nb_info_ssb.num, scan_pos, ctxs->m_prevStates, &decisions[0], sigCtxOffsetNext, gtxCtxOffsetNext, next_nb_info_ssb, 4, false);
      updateState(&ctxs->m_currStates[1], next_nb_info_ssb.num, scan_pos, ctxs->m_prevStates, &decisions[1], sigCtxOffsetNext, gtxCtxOffsetNext, next_nb_info_ssb, 4, false);
      updateState(&ctxs->m_currStates[2], next_nb_info_ssb.num, scan_pos, ctxs->m_prevStates, &decisions[2], sigCtxOffsetNext, gtxCtxOffsetNext, next_nb_info_ssb, 4, false);
      updateState(&ctxs->m_currStates[3], next_nb_info_ssb.num, scan_pos, ctxs->m_prevStates, &decisions[3], sigCtxOffsetNext, gtxCtxOffsetNext, next_nb_info_ssb, 4, false);
    }

    if (spt == SCAN_SOCSBB) {
      SWAP(ctxs->m_skipStates, ctxs->m_prevStates, depquant_state*);
    }
  }
}


int uvg_dep_quant(
  const encoder_state_t* const state,
  const cu_info_t* const cur_tu,
  const cu_loc_t* const cu_loc,
  const coeff_t* srcCoeff,
  coeff_t* coeff_out,
  const color_t compID,
  enum uvg_tree_type tree_type,
  int* absSum,
  const bool enableScalingLists)
{
  const encoder_control_t* const encoder = state->encoder_control;
  //===== reset / pre-init =====
  const int baseLevel = 4;
  context_store dep_quant_context;
  dep_quant_context.m_currStates = &dep_quant_context.m_allStates[0];
  dep_quant_context.m_prevStates = &dep_quant_context.m_allStates[4];
  dep_quant_context.m_skipStates = &dep_quant_context.m_allStates[8];
  
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
  const uint32_t* const cg_scan     = uvg_get_scan_order_table(SCAN_GROUP_UNGROUPED,0,log2_tr_width,log2_tr_height);

  int32_t qp_scaled = uvg_get_scaled_qp(compID, state->qp, (encoder->bitdepth - 8) * 6, encoder->qp_map[0]);
  qp_scaled = is_ts ? MAX(qp_scaled, 4 + 6 * MIN_QP_PRIME_TS) : qp_scaled;
  bool needs_block_size_trafo_scale = is_ts && ((log2_tr_height + log2_tr_width) % 2 == 1);
  needs_block_size_trafo_scale |= 0; // Non log2 block size

    const int32_t scalinglist_type = (cur_tu->type == CU_INTRA ? 0 : 3) + (int8_t)compID;
  const int32_t *q_coeff = encoder->scaling_list.quant_coeff[log2_tr_width][log2_tr_height][scalinglist_type][qp_scaled % 6];
  const int32_t transform_shift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - ((log2_tr_height + log2_tr_width) >> 1) - needs_block_size_trafo_scale; //!< Represents scaling through forward transform
  const int64_t q_bits = QUANT_SHIFT + qp_scaled / 6 + (is_ts ? 0 : transform_shift );
  const int32_t add = ((state->frame->slicetype == UVG_SLICE_I) ? 171 : 85) << (q_bits - 9);
  
  init_quant_block(state, &dep_quant_context.m_quant, cur_tu, log2_tr_width, log2_tr_height, compID, needs_block_size_trafo_scale, -1);
  
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
    firstTestPos =((width == 4 && height == 4) || (width == 8 && height == 8)) ? 7 : 15;
  }
  const int32_t default_quant_coeff = uvg_g_quant_scales[needs_block_size_trafo_scale][qp_scaled % 6];
  const coeff_t thres                          = 4 << q_bits;
  for (; firstTestPos >= 0; firstTestPos--) {
    coeff_t thresTmp = (enableScalingLists) ? (thres / (4 * q_coeff[scan[firstTestPos]])) : (thres / (4 * default_quant_coeff));
    if (abs(srcCoeff[scan[firstTestPos]]) > thresTmp) {
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

  reset_common_context(&dep_quant_context.m_common_context, &rate_estimator, (width * height) >> 4, numCoeff);
  dep_quant_context.m_common_context.m_nbInfo = encoder->m_scanId2NbInfoOutArray[log2_tr_width][log2_tr_height];
  

  int effectHeight = MIN(32, effHeight);
  int effectWidth = MIN(32, effWidth);
  for (int k = 0; k < 12; k++) {
    depquant_state_init(&dep_quant_context.m_allStates[k], rate_estimator.m_sigFracBits[0][0], rate_estimator.m_gtxFracBits[0]);
    dep_quant_context.m_allStates[k].effHeight = effectHeight;
    dep_quant_context.m_allStates[k].effWidth = effectWidth;
  }
  depquant_state_init(&dep_quant_context.m_startState, rate_estimator.m_sigFracBits[0][0], rate_estimator.m_gtxFracBits[0]);
  dep_quant_context.m_startState.effHeight = effectHeight;
  dep_quant_context.m_startState.effWidth = effectWidth;


  const uint32_t height_in_sbb = MAX(height >> 2, 1);
  const uint32_t width_in_sbb = MAX(width >> 2, 1);
  //===== populate trellis =====
  for (int scanIdx = firstTestPos; scanIdx >= 0; scanIdx--) {
    uint32_t blkpos = scan[scanIdx];
    uint32_t  pos_y = blkpos >> log2_tr_width;
    uint32_t  pos_x = blkpos - (pos_y << log2_tr_width);

    uint32_t cg_blockpos = scanIdx ? cg_scan[(scanIdx -1) >> 4] : 0;
    uint32_t cg_pos_y = cg_blockpos / height_in_sbb;
    uint32_t cg_pos_x = cg_blockpos - (cg_pos_y * height_in_sbb);
    uint32_t diag = cg_pos_y + cg_pos_x;

    uint32_t sig_ctx_offset = compID == COLOR_Y ? (diag < 2 ? 8 : diag < 5 ? 4 : 0) : (diag < 2 ? 4 : 0);
    uint32_t gtx_ctx_offset = compID == COLOR_Y ? (diag < 1 ? 16 : diag < 3 ? 11 : diag < 10 ? 6 : 1) : (diag < 1 ? 6 : 1);

    uint32_t nextSbbRight = (cg_pos_x < width_in_sbb - 1 ? cg_blockpos + 1 : 0);
    uint32_t nextSbbBelow = (cg_pos_y < height_in_sbb - 1 ? cg_blockpos + width_in_sbb : 0);

    if (enableScalingLists) {
      init_quant_block(state, &dep_quant_context.m_quant, cur_tu, log2_tr_width, log2_tr_height, compID, needs_block_size_trafo_scale, q_coeff[blkpos]);

      xDecideAndUpdate(
        &rate_estimator,
        &dep_quant_context,
        abs(srcCoeff[blkpos]),
        scanIdx,
        pos_x,
        pos_y,
        sig_ctx_offset,
        gtx_ctx_offset,
        width_in_sbb,
        height_in_sbb,
        nextSbbRight || nextSbbBelow,
        encoder->m_scanId2NbInfoSbbArray[log2_tr_width][log2_tr_height][scanIdx ? scanIdx - 1 : 0],
        (zeroOut && (pos_x >= effWidth || pos_y >= effHeight)),
        q_coeff[blkpos],
        effectWidth,
        effectHeight
        ); //tu.cu->slice->getReverseLastSigCoeffFlag());
    } else {
      xDecideAndUpdate(
        &rate_estimator,
        &dep_quant_context,
        abs(srcCoeff[blkpos]),
        scanIdx,
        pos_x,
        pos_y,
        sig_ctx_offset,
        gtx_ctx_offset,
        width_in_sbb,
        height_in_sbb,
        nextSbbRight || nextSbbBelow,
        encoder->m_scanId2NbInfoSbbArray[log2_tr_width][log2_tr_height][scanIdx ? scanIdx - 1 : 0],
        (zeroOut && (pos_x >= effWidth || pos_y >= effHeight)),
        default_quant_coeff,
        effectWidth,
        effectHeight); //tu.cu->slice->getReverseLastSigCoeffFlag());
      }
  }

  //===== find best path =====
  Decision decision    = {INT64_MAX, -1, -2};
  int64_t  minPathCost = 0;
  for (int8_t stateId = 0; stateId < 4; stateId++) {
    int64_t pathCost = dep_quant_context.m_trellis[0][stateId].rdCost;
    if (pathCost < minPathCost) {
      decision.prevId = stateId;
      minPathCost     = pathCost;
    }
  }

  //===== backward scanning =====
  int scanIdx = 0;
  for (; decision.prevId >= 0; scanIdx++) {
    decision       = dep_quant_context.m_trellis[scanIdx][decision.prevId];
    int32_t blkpos = scan[scanIdx];
    coeff_out[blkpos] = (srcCoeff[blkpos] < 0 ? -decision.absLevel : decision.absLevel);
    *absSum += decision.absLevel;
  }
  return *absSum;
}


void uvg_dep_quant_dequant(
  const encoder_state_t* const state,
  const cu_info_t* const cur_tu,
  const cu_loc_t* const cu_loc,
  const color_t compID,
  coeff_t* quant_coeff,
  coeff_t * coeff, 
  bool enableScalingLists)
{
  const encoder_control_t* const encoder = state->encoder_control;
  const uint32_t  width = compID == COLOR_Y ? cu_loc->width : cu_loc->chroma_width;
  const uint32_t  height = compID == COLOR_Y ? cu_loc->height : cu_loc->chroma_height;

  const int       numCoeff = width * height;
  
  const uint32_t  log2_tr_width = uvg_g_convert_to_log2[width];
  const uint32_t  log2_tr_height = uvg_g_convert_to_log2[height];
  const uint32_t* const scan = uvg_get_scan_order_table(SCAN_GROUP_4X4, 0, log2_tr_width, log2_tr_height);
  bool needs_block_size_trafo_scale =((log2_tr_height + log2_tr_width) % 2 == 1);
  needs_block_size_trafo_scale |= 0; // Non log2 block size

  //----- reset coefficients and get last scan index -----
  memset(coeff, 0, numCoeff * sizeof(coeff_t));
  int lastScanIdx = -1;
  for (int scanIdx = numCoeff - 1; scanIdx >= 0; scanIdx--)
  {
    if (quant_coeff[scan[scanIdx]])
    {
      lastScanIdx = scanIdx;
      break;
    }
  }
  if (lastScanIdx < 0)
  {
    return;
  }

  //----- set dequant parameters -----
  const int         qpDQ = uvg_get_scaled_qp(compID, state->qp, (encoder->bitdepth - 8) * 6, encoder->qp_map[0]); + 1;
  const int         qpPer = qpDQ / 6;
  const int         qpRem = qpDQ - 6 * qpPer;
  const int         channelBitDepth = encoder->bitdepth;
  const int         maxLog2TrDynamicRange = MAX_TR_DYNAMIC_RANGE;
  const coeff_t      minTCoeff = -(1 << maxLog2TrDynamicRange);
  const coeff_t      maxTCoeff = (1 << maxLog2TrDynamicRange) - 1;
  const int         transformShift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - ((log2_tr_height + log2_tr_width) >> 1) - needs_block_size_trafo_scale;
  int  shift = IQUANT_SHIFT + 1 - qpPer - transformShift + (enableScalingLists ? 4 : 0);
  int  invQScale = uvg_g_inv_quant_scales[needs_block_size_trafo_scale ? 1 : 0][qpRem];
  int  add = (shift < 0) ? 0 : ((1 << shift) >> 1);
  int32_t scalinglist_type = (cur_tu->type == CU_INTRA ? 0 : 3) + (int8_t)(compID);

  const int32_t* dequant_coef = encoder->scaling_list.de_quant_coeff[log2_tr_width][log2_tr_height][scalinglist_type][qpDQ % 6];
  //----- dequant coefficients -----
  for (int state = 0, scanIdx = lastScanIdx; scanIdx >= 0; scanIdx--)
  {
    const unsigned  rasterPos = scan[scanIdx];
    const coeff_t level = quant_coeff[rasterPos];
    if (level)
    {
      if (enableScalingLists)
      {
        invQScale = dequant_coef[rasterPos];//scalingfactor*levelScale
      }
      if (shift < 0 && (enableScalingLists || scanIdx == lastScanIdx))
      {
        invQScale <<= -shift;
      }
      int  qIdx = (level << 1) + (level > 0 ? -(state >> 1) : (state >> 1));
      int64_t  nomTCoeff = ((int64_t)qIdx * (int64_t)invQScale + add) >> ((shift < 0) ? 0 : shift);
      coeff[rasterPos] = (coeff_t)CLIP(minTCoeff, maxTCoeff, nomTCoeff);
    }
    state = (32040 >> ((state << 2) + ((level & 1) << 1))) & 3;   // the 16-bit value "32040" represent the state transition table
  }
}
