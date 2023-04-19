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
#include <immintrin.h>



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
  uint8_t* sbbFlags;
  uint8_t* levels;
} SbbCtx;


typedef struct
{
  int32_t absLevel[4];
  int64_t deltaDist[4];
} PQData;

typedef struct
{
  int64_t ALIGNED(32) rdCost[8];
  int32_t ALIGNED(32) absLevel[8];
  int32_t ALIGNED(32) prevId[8];
} Decision;


typedef struct
{
  const NbInfoOut* m_nbInfo;
  uint32_t m_sbbFlagBits[2][2];
  SbbCtx m_allSbbCtx[8];
  int m_curr_sbb_ctx_offset;
  int m_prev_sbb_ctx_offset;
  uint8_t sbb_memory[8 * 1024];
  uint8_t level_memory[8* TR_MAX_WIDTH * TR_MAX_WIDTH];
  int num_coeff;
} common_context;


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
  uint32_t *m_sigFracBitsArray[12];
  int32_t *m_gtxFracBitsArray[21];
  common_context* m_commonCtx;

  unsigned effWidth;
  unsigned effHeight;
} depquant_state;

typedef struct
{
  int64_t         ALIGNED(32) m_rdCost[12];
  uint16_t        ALIGNED(32) m_absLevelsAndCtxInit[12][24]; // 16x8bit for abs levels + 16x16bit for ctx init id
  int8_t          ALIGNED(16) m_numSigSbb[12];
  int             ALIGNED(32) m_remRegBins[12];
  int8_t          ALIGNED(16) m_refSbbCtxId[12];
  uint32_t        ALIGNED(32) m_sbbFracBits[12][2];
  uint32_t        ALIGNED(32) m_sigFracBits[12][2];
  int32_t         ALIGNED(32) m_coeffFracBits[12][6];
  int8_t          ALIGNED(16) m_goRicePar[12];
  int8_t          ALIGNED(16) m_goRiceZero[12];
  int8_t          ALIGNED(16) m_stateId[12];
  uint32_t        ALIGNED(32) m_sigFracBitsArray[12][12][2];
  int32_t         ALIGNED(32) m_gtxFracBitsArray[21][6];
  common_context* m_commonCtx;

  unsigned effWidth;
  unsigned effHeight;

  bool all_gte_four;
  bool all_lt_four;
} all_depquant_states;

typedef struct
{
    common_context  m_common_context;
    all_depquant_states m_allStates;
    int m_curr_state_offset;
    int m_prev_state_offset;
    int m_skip_state_offset;
    depquant_state       m_startState;
    quant_block*   m_quant;
    Decision    m_trellis[TR_MAX_WIDTH * TR_MAX_WIDTH];
} context_store;


int uvg_init_nb_info(encoder_control_t * encoder) {
  memset(encoder->m_scanId2NbInfoSbbArray, 0, sizeof(encoder->m_scanId2NbInfoSbbArray));
  memset(encoder->m_scanId2NbInfoOutArray, 0, sizeof(encoder->m_scanId2NbInfoOutArray));
  memset(encoder->scan_info, 0, sizeof(encoder->scan_info));
  for (int hd = 0; hd <= 6; hd++)
  {

    uint32_t raster2id[64 * 64] = {0};

    for (int vd = 0; vd <= 6; vd++)
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
      const uint32_t* const cg_scan = uvg_get_scan_order_table(SCAN_GROUP_UNGROUPED, 0, hd, vd);
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
      encoder->scan_info[hd][vd] = MALLOC(struct dep_quant_scan_info, totalValues);
      if (encoder->scan_info[hd][vd] == NULL) {
        return 0;
      }


      for (uint32_t scanId = 0; scanId < totalValues; scanId++)
      {
        raster2id[scanId2RP[scanId]] = scanId;
      }
      const uint32_t height_in_sbb = MAX(blockHeight >> 2, 1);
      const uint32_t width_in_sbb = MAX(blockWidth >> 2, 1);

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
        uint32_t cg_pos = cg_scan[scanId >> 4];

        uint32_t blkpos_next = scanId2RP[scanId ? scanId - 1 : 0];
        uint32_t  pos_y_next = blkpos_next >> hd;
        uint32_t  pos_x_next = blkpos_next - (pos_y_next << hd);
        uint32_t cg_blockpos_next = scanId ? cg_scan[(scanId - 1) >> 4] : 0;
        uint32_t cg_pos_y_next = cg_blockpos_next / width_in_sbb;
        uint32_t cg_pos_x_next = cg_blockpos_next - (cg_pos_y_next * width_in_sbb);
        uint32_t diag = pos_y_next + pos_x_next;
        

        uint32_t nextSbbRight = (cg_pos_x_next < width_in_sbb - 1 ? cg_blockpos_next + 1 : 0);
        uint32_t nextSbbBelow = (cg_pos_y_next < height_in_sbb - 1 ? cg_blockpos_next + width_in_sbb : 0);
        encoder->scan_info[hd][vd][scanId].pos_x = pos_x;
        encoder->scan_info[hd][vd][scanId].pos_y = pos_y;
        encoder->scan_info[hd][vd][scanId].sig_ctx_offset[0] = (diag < 2 ? 8 : diag < 5 ? 4 : 0);
        encoder->scan_info[hd][vd][scanId].sig_ctx_offset[1] = (diag < 2 ? 4 : 0);
        encoder->scan_info[hd][vd][scanId].gtx_ctx_offset[0] = (diag < 1 ? 16 : diag < 3 ? 11 : diag < 10 ? 6 : 1);
        encoder->scan_info[hd][vd][scanId].gtx_ctx_offset[1] = (diag < 1 ? 6 : 1);
        encoder->scan_info[hd][vd][scanId].cg_pos = cg_pos;
        encoder->scan_info[hd][vd][scanId].next_sbb_right = nextSbbRight;
        encoder->scan_info[hd][vd][scanId].next_sbb_below = nextSbbBelow;
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
      if(encoder->scan_info[hd][vd]) FREE_POINTER(encoder->scan_info[hd][vd]);
    }
  }
}


static INLINE int ceil_log2(uint64_t x)
{
  static const uint64_t t[6] = { 0xFFFFFFFF00000000ull, 0x00000000FFFF0000ull, 0x000000000000FF00ull, 0x00000000000000F0ull, 0x000000000000000Cull, 0x0000000000000002ull };
  int y = (((x & (x - 1)) == 0) ? 0 : 1);
  int j = 32;
  for (int i = 0; i < 6; i++)
  {
    int k = (((x & t[i]) == 0) ? 0 : j);
    y += k;
    x >>= k;
    j >>= 1;
  }
  return y;
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
  double     lambda = color == COLOR_Y ? state->lambda : state->c_lambda;

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
  qp->m_thresLast = (((int64_t)(4) << (int64_t)qp->m_QShift));
  qp->m_thresSSbb = (((int64_t)(3) << (int64_t)qp->m_QShift));
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
  const int dfShift = ceil_log2(pow2dfShift);
  qp->m_DistShift = 62 + qp->m_QShift - 2 * maxLog2TrDynamicRange - dfShift;
  qp->m_DistAdd = ((int64_t)(1) << qp->m_DistShift) >> 1;
  qp->m_DistStepAdd = (int64_t)(nomDistFactor * (double)((int64_t)(1) << (qp->m_DistShift + qp->m_QShift)) + .5);
  qp->m_DistOrgFact = (int64_t)(nomDistFactor * (double)((int64_t)(1) << (qp->m_DistShift + 1)) + .5);
  qp->needs_init = false;
}

static void reset_common_context(common_context* ctx, const rate_estimator_t * rate_estimator, int numSbb, int num_coeff)
{
  //memset(&ctx->m_nbInfo, 0, sizeof(ctx->m_nbInfo));
  memcpy(&ctx->m_sbbFlagBits, &rate_estimator->m_sigSbbFracBits, sizeof(rate_estimator->m_sigSbbFracBits));
  uint8_t*  next_sbb_memory   = ctx->sbb_memory;
  uint8_t*  next_level_memory   = ctx->level_memory;
  for (int k = 0; k < 8; k++, next_sbb_memory += numSbb, next_level_memory += num_coeff) {
    ctx->m_allSbbCtx[k].sbbFlags = next_sbb_memory;
    ctx->m_allSbbCtx[k].levels = next_level_memory;
  }
  ctx->m_curr_sbb_ctx_offset = 0;
  ctx->m_prev_sbb_ctx_offset = 4;
  ctx->num_coeff = num_coeff;
}

static void init_rate_esimator(rate_estimator_t * rate_estimator, const cabac_data_t * const ctx, color_t color)
{
  const cabac_ctx_t * base_ctx = color == COLOR_Y ? ctx->ctx.sig_coeff_group_model : (ctx->ctx.sig_coeff_group_model + 2);
  for (unsigned ctxId = 0; ctxId < SM_MAX_NUM_SIG_SBB_CTX; ctxId++) {
    rate_estimator->m_sigSbbFracBits[ctxId][0] = CTX_ENTROPY_BITS(&base_ctx[ctxId], 0);
    rate_estimator->m_sigSbbFracBits[ctxId][1] = CTX_ENTROPY_BITS(&base_ctx[ctxId], 1);
  }
  unsigned numCtx = (color == COLOR_Y ? 12 : 8);
  for (unsigned ctxSetId = 0; ctxSetId < SM_NUM_CTX_SETS_SIG; ctxSetId++) {
    base_ctx = color == COLOR_Y ? ctx->ctx.cu_sig_model_luma[ctxSetId] : ctx->ctx.cu_sig_model_chroma[ctxSetId];
    for (unsigned ctxId = 0; ctxId < numCtx; ctxId++) {
      rate_estimator->m_sigFracBits[ctxSetId][ctxId][0] = CTX_ENTROPY_BITS(&base_ctx[ctxId], 0);
      rate_estimator->m_sigFracBits[ctxSetId][ctxId][1] = CTX_ENTROPY_BITS(&base_ctx[ctxId], 1);
    }
  }
  
  numCtx    = (color == COLOR_Y? 21 : 11);
  for (unsigned ctxId = 0; ctxId < numCtx; ctxId++) {
    const cabac_ctx_t * par_ctx = color == COLOR_Y ? &ctx->ctx.cu_parity_flag_model_luma[ctxId] : &ctx->ctx.cu_parity_flag_model_chroma[ctxId];
    const cabac_ctx_t * gt2_ctx = color == COLOR_Y ? &ctx->ctx.cu_gtx_flag_model_luma[0][ctxId] : &ctx->ctx.cu_gtx_flag_model_chroma[0][ctxId];
    const cabac_ctx_t * gt1_ctx = color == COLOR_Y ? &ctx->ctx.cu_gtx_flag_model_luma[1][ctxId] : &ctx->ctx.cu_gtx_flag_model_chroma[1][ctxId];

    int32_t* cb = rate_estimator->m_gtxFracBits[ctxId];
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
  const cu_info_t* const cur_tu,
  const int width,
  const int height,
  rate_estimator_t* rate_estimator,
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
      uint32_t nTus = uvg_get_isp_split_num(1 << cur_tu->log2_width, 1 << cur_tu->log2_height, cur_tu->intra.isp_mode, true);
      bool     isLastSubPartition = cur_tu->intra.isp_index +1 == nTus; //TODO: isp check
      if (isLastSubPartition) {
        lastCbfIsInferred = cur_tu->intra.isp_cbfs == 0;
      }
      if (!lastCbfIsInferred) {
        prevLumaCbf = cur_tu->intra.isp_index != 0 && (cur_tu->intra.isp_cbfs & (1 << (cur_tu->intra.isp_index - 1)));
      }
      const cabac_ctx_t * const cbf_ctx = &state->search_cabac.ctx.qt_cbf_model_luma[2 + prevLumaCbf];
      cbfDeltaBits = lastCbfIsInferred ? 0 : (int32_t)CTX_ENTROPY_BITS(cbf_ctx, 1) - (int32_t)CTX_ENTROPY_BITS(cbf_ctx, 0);
    }
    else {
      const cabac_ctx_t* cbf_ctx;
      switch (compID) {
        case COLOR_Y:
          cbf_ctx = &state->search_cabac.ctx.qt_cbf_model_luma[0];
          break;
        case COLOR_U:
          cbf_ctx = &state->search_cabac.ctx.qt_cbf_model_cb[0];
          break;
        case COLOR_V:
          cbf_ctx = &state->search_cabac.ctx.qt_cbf_model_cr[cbf_is_set(cur_tu->cbf, COLOR_U)];
          break;
      }
      cbfDeltaBits = compID != COLOR_Y && cur_tu->joint_cb_cr ? 0 : (int32_t)CTX_ENTROPY_BITS(cbf_ctx, 1) - (int32_t)CTX_ENTROPY_BITS(cbf_ctx, 0);
    }
     
  }

  static const unsigned prefixCtx[] = {0, 0, 0, 3, 6, 10, 15, 21};
  uint32_t              ctxBits[14];
  for (unsigned xy = 0; xy < 2; xy++) {
    int32_t        bitOffset  = (xy ? cbfDeltaBits : 0);
    int32_t*       lastBits   = (xy ? rate_estimator->m_lastBitsY : rate_estimator->m_lastBitsX);
    const unsigned size = (xy ? (height) : (width));
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
  state->m_rdCost = INT64_MAX >> 1;
  state->m_numSigSbb = 0;
  state->m_remRegBins = 4; // just large enough for last scan pos
  state->m_refSbbCtxId = -1;
  state->m_sigFracBits[0] = sig_frac_bits[0];
  state->m_sigFracBits[1] = sig_frac_bits[1];
  memcpy(state->m_coeffFracBits, gtx_frac_bits, sizeof(state->m_coeffFracBits));
  state->m_goRicePar = 0;
  state->m_goRiceZero = 0;

  state->m_sbbFracBits[0] = 0;
  state->m_sbbFracBits[1] = 0;
}

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



static void check_rd_costs_avx2(const all_depquant_states* const state, const enum ScanPosType spt, const PQData* pqDataA, Decision* decisions, int start)
{
  int64_t temp_rd_cost_a[4] = {0, 0, 0, 0};
  int64_t temp_rd_cost_b[4] = {0, 0, 0, 0};
  int64_t temp_rd_cost_z[4] = {0, 0, 0, 0};

  __m256i pq_a_delta_dist = _mm256_setr_epi64x(pqDataA->deltaDist[0], pqDataA->deltaDist[0], pqDataA->deltaDist[3], pqDataA->deltaDist[3]);
  __m256i pq_b_delta_dist = _mm256_setr_epi64x(pqDataA->deltaDist[2], pqDataA->deltaDist[2], pqDataA->deltaDist[1], pqDataA->deltaDist[1]);

  __m256i rd_cost_a = _mm256_load_si256((__m256i const*)&state->m_rdCost[start]);
  __m256i rd_cost_b = rd_cost_a;
  __m256i rd_cost_z = rd_cost_a;

  rd_cost_a = _mm256_add_epi64(rd_cost_a, pq_a_delta_dist);
  rd_cost_b = _mm256_add_epi64(rd_cost_b, pq_b_delta_dist);


  if (state->all_gte_four) {
    if (pqDataA->absLevel[0] < 4 && pqDataA->absLevel[3] < 4) {
      __m128i offsets = _mm_set_epi32(18 + pqDataA->absLevel[3], 12 + pqDataA->absLevel[3], 6 + pqDataA->absLevel[0], 0 + pqDataA->absLevel[0]);
      __m128i coeff_frac_bits = _mm_i32gather_epi32(&state->m_coeffFracBits[start][0], offsets, 4);
      __m256i ext_frac_bits = _mm256_cvtepi32_epi64(coeff_frac_bits);
      rd_cost_a = _mm256_add_epi64(rd_cost_a, ext_frac_bits);
    } else if (pqDataA->absLevel[0] >= 4 && pqDataA->absLevel[3] >= 4) {
      __m128i value = _mm_set_epi32((pqDataA->absLevel[3] - 4) >> 1, (pqDataA->absLevel[3] - 4) >> 1, (pqDataA->absLevel[0] - 4) >> 1, (pqDataA->absLevel[0] - 4) >> 1);

      __m128i offsets = _mm_set_epi32(18 + pqDataA->absLevel[3], 12 + pqDataA->absLevel[3], 6 + pqDataA->absLevel[0], 0 + pqDataA->absLevel[0]);
      __m128i t = _mm_slli_epi32(value, 1);
      offsets = _mm_sub_epi32(offsets, t);
      __m128i coeff_frac_bits = _mm_i32gather_epi32(state->m_coeffFracBits[start], offsets, 4);

      __m128i max_rice = _mm_set1_epi32(31);
      value = _mm_min_epi32(value, max_rice);
      __m128i go_rice_tab = _mm_cvtepi8_epi32(_mm_loadu_si32(&state->m_goRicePar[start]));
      go_rice_tab = _mm_slli_epi32(go_rice_tab, 5);
      value = _mm_add_epi32(value, go_rice_tab);

      __m128i temp = _mm_add_epi32(coeff_frac_bits, _mm_i32gather_epi32(&g_goRiceBits[0][0], value, 4));
      rd_cost_a = _mm256_add_epi64(rd_cost_a, _mm256_cvtepi32_epi64(temp));
    } else {
      const int pqAs[4] = {0, 0, 3, 3};
      ALIGNED(32) int64_t rd_costs[4] = {0, 0, 0, 0}; 
      for (int i = 0; i < 4; i++) {
        const int      state_offset = start + i;
        const int      pqA = pqAs[i];
        const int32_t* goRiceTab = g_goRiceBits[state->m_goRicePar[state_offset]];
        if (pqDataA->absLevel[pqA] < 4) {
          rd_costs[i] = state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqA]];
        } else {
          const coeff_t value = (pqDataA->absLevel[pqA] - 4) >> 1;
          rd_costs[i] += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqA] - (value << 1)] + goRiceTab[value < RICEMAX ? value : RICEMAX - 1];
        }
      }
      rd_cost_a = _mm256_add_epi64(rd_cost_a, _mm256_loadu_si256((__m256i const *)&rd_costs[0]));
    }

    if (pqDataA->absLevel[1] < 4 && pqDataA->absLevel[2] < 4) {
      __m128i offsets = _mm_set_epi32(18 + pqDataA->absLevel[1], 12 + pqDataA->absLevel[1], 6 + pqDataA->absLevel[2], 0 + pqDataA->absLevel[2]);
      __m128i coeff_frac_bits = _mm_i32gather_epi32(state->m_coeffFracBits[start], offsets, 4);
      __m256i ext_frac_bits = _mm256_cvtepi32_epi64(coeff_frac_bits);
      rd_cost_b = _mm256_add_epi64(rd_cost_b, ext_frac_bits);
    } else if (pqDataA->absLevel[1] >= 4 && pqDataA->absLevel[2] >= 4) {
      __m128i value = _mm_set_epi32((pqDataA->absLevel[1] - 4) >> 1, (pqDataA->absLevel[1] - 4) >> 1, (pqDataA->absLevel[2] - 4) >> 1, (pqDataA->absLevel[2] - 4) >> 1);

      __m128i offsets = _mm_set_epi32(18 + pqDataA->absLevel[1], 12 + pqDataA->absLevel[1], 6 + pqDataA->absLevel[2], 0 + pqDataA->absLevel[2]);
      __m128i t = _mm_slli_epi32(value, 1);
      offsets = _mm_sub_epi32(offsets, t);
      __m128i coeff_frac_bits = _mm_i32gather_epi32(state->m_coeffFracBits[start], offsets, 4);

      __m128i max_rice = _mm_set1_epi32(31);
      value = _mm_min_epi32(value, max_rice);
      __m128i go_rice_tab = _mm_cvtepi8_epi32(_mm_loadu_si32(&state->m_goRicePar[start]));
      go_rice_tab = _mm_slli_epi32(go_rice_tab, 5);
      value = _mm_add_epi32(value, go_rice_tab);

      __m128i temp = _mm_add_epi32(coeff_frac_bits, _mm_i32gather_epi32(&g_goRiceBits[0][0], value, 4));
      rd_cost_b = _mm256_add_epi64(rd_cost_b, _mm256_cvtepi32_epi64(temp));
    } else {
      const int pqBs[4] = {2, 2, 1, 1};
      int64_t rd_costs[4] = {0, 0, 0, 0}; 
      for (int i = 0; i < 4; i++) {
        const int      state_offset = start + i;
        const int      pqB = pqBs[i];
        const int32_t* goRiceTab = g_goRiceBits[state->m_goRicePar[state_offset]];
        if (pqDataA->absLevel[pqB] < 4) {
          rd_costs[i] = state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqB]];
        } else {
          const coeff_t value = (pqDataA->absLevel[pqB] - 4) >> 1;
          rd_costs[i] += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqB] - (value << 1)] + goRiceTab[value < RICEMAX ? value : RICEMAX - 1];
        }
      }
      rd_cost_b =
        _mm256_add_epi64(rd_cost_b, _mm256_loadu_si256((__m256i const *) & rd_costs[0]));
    }

    if (spt == SCAN_ISCSBB) {
      __m256i original = _mm256_loadu_si256((__m256i const*)state->m_sigFracBits[start]);
      __m256i even_mask = _mm256_setr_epi32(0, 2, 4, 6, -1, -1, -1, -1);
      __m256i odd_mask = _mm256_setr_epi32(1, 3, 5, 7, -1, -1, -1, -1);
      __m256i even = _mm256_permutevar8x32_epi32(original, even_mask);
      __m256i odd = _mm256_permutevar8x32_epi32(original, odd_mask);
      __m256i even_64 = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(even, 0));
      __m256i odd_64 = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(odd, 0));
      rd_cost_a = _mm256_add_epi64(rd_cost_a, odd_64);
      rd_cost_b = _mm256_add_epi64(rd_cost_b, odd_64);
      rd_cost_z = _mm256_add_epi64(rd_cost_z, even_64);
    } else if (spt == SCAN_SOCSBB) {
      __m256i original = _mm256_loadu_si256((__m256i const*)state->m_sigFracBits[start]);
      __m256i even_mask = _mm256_setr_epi32(0, 2, 4, 6, -1, -1, -1, -1);
      __m256i odd_mask = _mm256_setr_epi32(1, 3, 5, 7, -1, -1, -1, -1);
      __m256i even = _mm256_permutevar8x32_epi32(original, even_mask);
      __m256i odd = _mm256_permutevar8x32_epi32(original, odd_mask);
      __m256i m_sigFracBits_0 = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(even, 0));
      __m256i m_sigFracBits_1 = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(odd, 0));

      original = _mm256_loadu_si256((__m256i const*)state->m_sbbFracBits[start]);
      odd = _mm256_permutevar8x32_epi32(original, odd_mask);
      __m256i m_sbbFracBits_1 = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(odd, 0));

      
      rd_cost_a = _mm256_add_epi64(rd_cost_a, m_sbbFracBits_1);
      rd_cost_b = _mm256_add_epi64(rd_cost_b, m_sbbFracBits_1);
      rd_cost_z = _mm256_add_epi64(rd_cost_z, m_sbbFracBits_1);

      rd_cost_a = _mm256_add_epi64(rd_cost_a, m_sigFracBits_1);
      rd_cost_b = _mm256_add_epi64(rd_cost_b, m_sigFracBits_1);
      rd_cost_z = _mm256_add_epi64(rd_cost_z, m_sigFracBits_0);
    }
    else {
      if (state->m_numSigSbb[start] && state->m_numSigSbb[start + 1] && state->m_numSigSbb[start + 2] && state->m_numSigSbb[start + 3]) {
        __m256i original = _mm256_loadu_si256((__m256i const*)state->m_sigFracBits[start]);
        __m256i even_mask = _mm256_setr_epi32(0, 2, 4, 6, -1, -1, -1, -1);
        __m256i odd_mask = _mm256_setr_epi32(1, 3, 5, 7, -1, -1, -1, -1);
        __m256i even = _mm256_permutevar8x32_epi32(original, even_mask);
        __m256i odd = _mm256_permutevar8x32_epi32(original, odd_mask);
        __m256i even_64 = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(even, 0));
        __m256i odd_64 = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(odd, 0));
        rd_cost_a = _mm256_add_epi64(rd_cost_a, odd_64);
        rd_cost_b = _mm256_add_epi64(rd_cost_b, odd_64);
        rd_cost_z = _mm256_add_epi64(rd_cost_z, even_64);     
      }
      else if (!state->m_numSigSbb[start] && !state->m_numSigSbb[start + 1] && !state->m_numSigSbb[start + 2] && !state->m_numSigSbb[start + 3]) {
        rd_cost_z = _mm256_setr_epi64x(decisions->rdCost[0], decisions->rdCost[0], decisions->rdCost[3], decisions->rdCost[3]);
      }

      else {
        const int ALIGNED(32) pqAs[4] = {0, 0, 3, 3};
        _mm256_store_si256((__m256i*)temp_rd_cost_a, rd_cost_a);
        _mm256_store_si256((__m256i*)temp_rd_cost_b, rd_cost_b);
        _mm256_store_si256((__m256i*)temp_rd_cost_z, rd_cost_z);
        for (int i = 0; i < 4; i++) {
          const int state_offset = start + i;
          if (state->m_numSigSbb[state_offset]) {
            temp_rd_cost_a[i] += state->m_sigFracBits[state_offset][1];
            temp_rd_cost_b[i] += state->m_sigFracBits[state_offset][1];
            temp_rd_cost_z[i] += state->m_sigFracBits[state_offset][0];
          } else {
            temp_rd_cost_z[i] = decisions->rdCost[pqAs[i]];
          }
        }
        rd_cost_a = _mm256_loadu_si256((__m256i*)temp_rd_cost_a);
        rd_cost_b = _mm256_loadu_si256((__m256i*)temp_rd_cost_b);
        rd_cost_z = _mm256_loadu_si256((__m256i*)temp_rd_cost_z);
      }
    }
  } else if (state->all_lt_four) {
    __m128i scale_bits = _mm_set1_epi32(1 << SCALE_BITS);
    __m128i max_rice = _mm_set1_epi32(31);
    __m128i go_rice_zero = _mm_cvtepi8_epi32(_mm_loadu_si128((const __m128i*)&state->m_goRiceZero[start]));
    // RD cost A
    {
      __m128i pq_abs_a = _mm_set_epi32(pqDataA->absLevel[3], pqDataA->absLevel[3], pqDataA->absLevel[0], pqDataA->absLevel[0]);
      __m128i cmp = _mm_cmpgt_epi32(pq_abs_a, go_rice_zero);
      
      __m128i go_rice_smaller = _mm_min_epi32(pq_abs_a, max_rice);

      __m128i other = _mm_sub_epi32(pq_abs_a, _mm_set1_epi32(1));

      __m128i selected = _mm_blendv_epi8(other, go_rice_smaller, cmp);


      __m128i go_rice_offset = _mm_cvtepi8_epi32(_mm_loadu_si128((__m128i*)&state->m_goRicePar[start]));
      go_rice_offset = _mm_slli_epi32(go_rice_offset, 5);

      __m128i offsets = _mm_add_epi32(selected, go_rice_offset);
      __m128i go_rice_tab = _mm_i32gather_epi32(&g_goRiceBits[0][0], offsets, 4);
      __m128i temp = _mm_add_epi32(go_rice_tab, scale_bits);

      rd_cost_a = _mm256_add_epi64(rd_cost_a, _mm256_cvtepi32_epi64(temp));
    }
    // RD cost b
    {
      __m128i pq_abs_b = _mm_set_epi32(pqDataA->absLevel[1], pqDataA->absLevel[1], pqDataA->absLevel[2], pqDataA->absLevel[2]);
      __m128i cmp = _mm_cmpgt_epi32(pq_abs_b, go_rice_zero);

      __m128i go_rice_smaller = _mm_min_epi32(pq_abs_b, max_rice);

      __m128i other = _mm_sub_epi32(pq_abs_b, _mm_set1_epi32(1));

      __m128i selected = _mm_blendv_epi8(other, go_rice_smaller, cmp);


      __m128i go_rice_offset = _mm_cvtepi8_epi32(_mm_loadu_si128((__m128i*)&state->m_goRicePar[start]));
      go_rice_offset = _mm_slli_epi32(go_rice_offset, 5);

      __m128i offsets = _mm_add_epi32(selected, go_rice_offset);
      __m128i go_rice_tab = _mm_i32gather_epi32(&g_goRiceBits[0][0], offsets, 4);
      __m128i temp = _mm_add_epi32(go_rice_tab, scale_bits);

      rd_cost_b = _mm256_add_epi64(rd_cost_b, _mm256_cvtepi32_epi64(temp));
    }
    // RD cost Z
    {
      __m128i go_rice_offset = _mm_cvtepi8_epi32(_mm_loadu_si128((__m128i*)&state->m_goRicePar[start]));
      go_rice_offset = _mm_slli_epi32(go_rice_offset, 5);

      go_rice_offset = _mm_add_epi32(go_rice_offset, go_rice_zero);
      __m128i go_rice_tab = _mm_i32gather_epi32(&g_goRiceBits[0][0], go_rice_offset, 4);
      rd_cost_z = _mm256_add_epi64(rd_cost_z, _mm256_cvtepi32_epi64(go_rice_tab));
    }
  } else {
    const int pqAs[4] = {0, 0, 3, 3};
    const int pqBs[4] = {2, 2, 1, 1};
    const int decision_a[4] = {0, 2, 1, 3};
    for (int i = 0; i < 4; i++) {
      const int      state_offset = start + i;
      const int32_t* goRiceTab = g_goRiceBits[state->m_goRicePar[state_offset]];
      const int pqA = pqAs[i];
      const int pqB = pqBs[i];
      int64_t rdCostA = state->m_rdCost[state_offset] + pqDataA->deltaDist[pqA];
      int64_t rdCostB = state->m_rdCost[state_offset] + pqDataA->deltaDist[pqB];
      int64_t rdCostZ = state->m_rdCost[state_offset];
      if (state->m_remRegBins[state_offset] >= 4) {
        if (pqDataA->absLevel[pqA] < 4) {
          rdCostA += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqA]];
        } else {
          const coeff_t value = (pqDataA->absLevel[pqA] - 4) >> 1;
          rdCostA += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqA] - (value << 1)] + goRiceTab[value < RICEMAX ? value : RICEMAX - 1];
        }
        if (pqDataA->absLevel[pqB] < 4) {
          rdCostB += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqB]];
        } else {
          const coeff_t value = (pqDataA->absLevel[pqB] - 4) >> 1;
          rdCostB += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqB] - (value << 1)] + goRiceTab[value < RICEMAX ? value : RICEMAX - 1];
        }
        if (spt == SCAN_ISCSBB) {
          rdCostA += state->m_sigFracBits[state_offset][1];
          rdCostB += state->m_sigFracBits[state_offset][1];
          rdCostZ += state->m_sigFracBits[state_offset][0];
        } else if (spt == SCAN_SOCSBB) {
          rdCostA += state->m_sbbFracBits[state_offset][1] + state->m_sigFracBits[state_offset][1];
          rdCostB += state->m_sbbFracBits[state_offset][1] + state->m_sigFracBits[state_offset][1];
          rdCostZ += state->m_sbbFracBits[state_offset][1] + state->m_sigFracBits[state_offset][0];
        } else if (state->m_numSigSbb[state_offset]) {
          rdCostA += state->m_sigFracBits[state_offset][1];
          rdCostB += state->m_sigFracBits[state_offset][1];
          rdCostZ += state->m_sigFracBits[state_offset][0];
        } else {
          rdCostZ = decisions->rdCost[decision_a[i]];
        }
      } else {
        rdCostA += (1 << SCALE_BITS) + goRiceTab[pqDataA->absLevel[pqA] <= state->m_goRiceZero[state_offset] ? pqDataA->absLevel[pqA] - 1 : (pqDataA->absLevel[pqA] < RICEMAX ? pqDataA->absLevel[pqA] : RICEMAX - 1)];
        rdCostB += (1 << SCALE_BITS) + goRiceTab[pqDataA->absLevel[pqB] <= state->m_goRiceZero[state_offset] ? pqDataA->absLevel[pqB] - 1 : (pqDataA->absLevel[pqB] < RICEMAX ? pqDataA->absLevel[pqB] : RICEMAX - 1)];
        rdCostZ += goRiceTab[state->m_goRiceZero[state_offset]];
      }
      temp_rd_cost_a[i] = rdCostA;
      temp_rd_cost_b[i] = rdCostB;
      temp_rd_cost_z[i] = rdCostZ;
    }
    rd_cost_a = _mm256_loadu_si256((__m256i*)temp_rd_cost_a);
    rd_cost_b = _mm256_loadu_si256((__m256i*)temp_rd_cost_b);
    rd_cost_z = _mm256_loadu_si256((__m256i*)temp_rd_cost_z);
  }
  rd_cost_a = _mm256_permute4x64_epi64(rd_cost_a, 216);
  rd_cost_b = _mm256_permute4x64_epi64(rd_cost_b, 141);
  rd_cost_z = _mm256_permute4x64_epi64(rd_cost_z, 216);
  __m256i rd_cost_decision = _mm256_load_si256((__m256i*)decisions->rdCost);

  __m256i decision_abs_coeff = _mm256_load_si256((__m256i*)decisions->absLevel);
  __m256i decision_prev_state = _mm256_load_si256((__m256i*)decisions->prevId);
  __m256i decision_data = _mm256_permute2x128_si256(decision_abs_coeff, decision_prev_state, 0x20);
  __m256i mask = _mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0);
  decision_data = _mm256_permutevar8x32_epi32(decision_data, mask);

  __m256i a_data = _mm256_set_epi32(3, pqDataA->absLevel[3], 1, pqDataA->absLevel[0], 2, pqDataA->absLevel[3], 0, pqDataA->absLevel[0]);
  __m256i b_data = _mm256_set_epi32(2, pqDataA->absLevel[1], 0, pqDataA->absLevel[2], 3, pqDataA->absLevel[1], 1, pqDataA->absLevel[2]);
  __m256i z_data = _mm256_set_epi32(3, 0, 1, 0, 2, 0, 0, 0);

  __m256i a_vs_b = _mm256_cmpgt_epi64(rd_cost_a, rd_cost_b);
  __m256i cheaper_first = _mm256_blendv_epi8(rd_cost_a, rd_cost_b, a_vs_b);
  __m256i cheaper_first_data = _mm256_blendv_epi8(a_data, b_data, a_vs_b);

  __m256i z_vs_decision = _mm256_cmpgt_epi64(rd_cost_z, rd_cost_decision);
  __m256i cheaper_second = _mm256_blendv_epi8(rd_cost_z, rd_cost_decision, z_vs_decision);
  __m256i cheaper_second_data = _mm256_blendv_epi8(z_data, decision_data, z_vs_decision);

  __m256i final_decision = _mm256_cmpgt_epi64(cheaper_first, cheaper_second);
  __m256i final_rd_cost = _mm256_blendv_epi8(cheaper_first, cheaper_second, final_decision);
  __m256i final_data = _mm256_blendv_epi8(cheaper_first_data, cheaper_second_data, final_decision);

  _mm256_store_si256((__m256i*)decisions->rdCost, final_rd_cost);
  final_data = _mm256_permutevar8x32_epi32(final_data, _mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0));
  _mm256_storeu2_m128i((__m128i *)decisions->prevId, (__m128i *)decisions->absLevel, final_data);
}


static void checkRdCosts(
  const all_depquant_states * const state,
  const enum ScanPosType            spt,
  const PQData *                    pqDataA,
  Decision *                        decisions,
  const int                         decisionA,
  const int                         decisionB,
  const int                         state_offset)
{
  const int pqA = decisionA && decisionB ? 3 : 0;
  const int pqB = decisionA && decisionB ? 1 : 2;
  const int32_t* goRiceTab = g_goRiceBits[state->m_goRicePar[state_offset]];
  int64_t rdCostA = state->m_rdCost[state_offset] + pqDataA->deltaDist[pqA];
  int64_t rdCostB = state->m_rdCost[state_offset] + pqDataA->deltaDist[pqB];
  int64_t rdCostZ = state->m_rdCost[state_offset];
  if (state->m_remRegBins[state_offset] >= 4) {
    if (pqDataA->absLevel[pqA] < 4) {
      rdCostA += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqA]];
    }
    else {
      const coeff_t value = (pqDataA->absLevel[pqA] - 4) >> 1;
      rdCostA += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqA] - (value << 1)] + goRiceTab[value < RICEMAX ? value : RICEMAX - 1];
    }
    if (pqDataA->absLevel[pqB] < 4) {
      rdCostB += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqB]];
    }
    else {
      const coeff_t value = (pqDataA->absLevel[pqB] - 4) >> 1;
      rdCostB += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqB] - (value << 1)] + goRiceTab[value < RICEMAX ? value : RICEMAX - 1];
    }
    if (spt == SCAN_ISCSBB) {
      rdCostA += state->m_sigFracBits[state_offset][1];
      rdCostB += state->m_sigFracBits[state_offset][1];
      rdCostZ += state->m_sigFracBits[state_offset][0];
    }
    else if (spt == SCAN_SOCSBB) {
      rdCostA += state->m_sbbFracBits[state_offset][1] + state->m_sigFracBits[state_offset][1];
      rdCostB += state->m_sbbFracBits[state_offset][1] + state->m_sigFracBits[state_offset][1];
      rdCostZ += state->m_sbbFracBits[state_offset][1] + state->m_sigFracBits[state_offset][0];
    }
    else if (state->m_numSigSbb[state_offset]) {
      rdCostA += state->m_sigFracBits[state_offset][1];
      rdCostB += state->m_sigFracBits[state_offset][1];
      rdCostZ += state->m_sigFracBits[state_offset][0];
    }
    else {
      rdCostZ = decisions->rdCost[decisionA];
    }
  }
  else {
    rdCostA += (1 << SCALE_BITS) + goRiceTab[pqDataA->absLevel[pqA] <= state->m_goRiceZero[state_offset]
                                      ? pqDataA->absLevel[pqA] - 1
                                      : (pqDataA->absLevel[pqA] < RICEMAX ? pqDataA->absLevel[pqA] : RICEMAX - 1)];
    rdCostB += (1 << SCALE_BITS) + goRiceTab[pqDataA->absLevel[pqB] <= state->m_goRiceZero[state_offset]
                                      ? pqDataA->absLevel[pqB] - 1
                                      : (pqDataA->absLevel[pqB] < RICEMAX ? pqDataA->absLevel[pqB] : RICEMAX - 1)];
    rdCostZ += goRiceTab[state->m_goRiceZero[state_offset]];
  }
  if (rdCostA < decisions->rdCost[decisionA]) {
    decisions->rdCost[decisionA] = rdCostA;
    decisions->absLevel[decisionA] = pqDataA->absLevel[pqA];
    decisions->prevId[decisionA] = state->m_stateId[state_offset];
  }
  if (rdCostZ < decisions->rdCost[decisionA]) {
    decisions->rdCost[decisionA] = rdCostZ;
    decisions->absLevel[decisionA] = 0;
    decisions->prevId[decisionA] = state->m_stateId[state_offset];
  }
  if (rdCostB < decisions->rdCost[decisionB]) {
    decisions->rdCost[decisionB] = rdCostB;
    decisions->absLevel[decisionB] = pqDataA->absLevel[pqB];
    decisions->prevId[decisionB] = state->m_stateId[state_offset];
  }
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


static const Decision startDec = { .rdCost = {INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2},
  .absLevel = {-1, -1, -1, -1, 0, 0, 0, 0}, .prevId = {-2, -2, -2, -2, 4, 5, 6, 7} };


static void xDecide(
  all_depquant_states* const all_states,
  depquant_state* const      m_startState,
  quant_block *              qp,
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
  check_rd_costs_avx2(all_states, spt, &pqData, decisions, prev_offset);
  //checkRdCosts(all_states, spt, &pqData, decisions, 0, 2, prev_offset + 0);
  //checkRdCosts(all_states, spt, &pqData, decisions, 2, 0, prev_offset + 1);
  //checkRdCosts(all_states, spt, &pqData, decisions, 1, 3, prev_offset + 2);
  //checkRdCosts(all_states, spt, &pqData, decisions, 3, 1, prev_offset + 3);
  if (spt == SCAN_EOCSBB) {
    checkRdCostSkipSbb(all_states, decisions, 0, skip_offset);
    checkRdCostSkipSbb(all_states, decisions, 1, skip_offset);
    checkRdCostSkipSbb(all_states, decisions, 2, skip_offset);
    checkRdCostSkipSbb(all_states, decisions, 3, skip_offset);
  }

  checkRdCostStart(m_startState, lastOffset, &pqData, decisions, 0);
  checkRdCostStart(m_startState, lastOffset, &pqData, decisions, 2);
}


static INLINE unsigned templateAbsCompare(coeff_t sum)
{
  int rangeIdx = 0;
  if (sum < g_riceT[0]) {
    rangeIdx = 0;
  }
  else if (sum < g_riceT[1]) {
    rangeIdx = 1;
  }
  else if (sum < g_riceT[2]) {
    rangeIdx = 2;
  }
  else if (sum < g_riceT[3]) {
    rangeIdx = 3;
  }
  else {
    rangeIdx = 4;
  }
  return g_riceShift[rangeIdx];
}

static INLINE void update_common_context(
  context_store*   ctxs,
  common_context * cc,
  const uint32_t   scan_pos,
  const uint32_t   cg_pos,
  const uint32_t   width_in_sbb,
  const uint32_t   height_in_sbb,
  const uint32_t   next_sbb_right,
  const uint32_t   next_sbb_below,
  const int        prev_state,
  const int        curr_state)
{
  const uint32_t numSbb = width_in_sbb * height_in_sbb;
  uint8_t* sbbFlags = cc->m_allSbbCtx[cc->m_curr_sbb_ctx_offset + (curr_state & 3)].sbbFlags;
  uint8_t* levels = cc->m_allSbbCtx[cc->m_curr_sbb_ctx_offset + (curr_state & 3)].levels;
  size_t setCpSize = cc->m_nbInfo[scan_pos - 1].maxDist * sizeof(uint8_t);
  if (prev_state != -1 && ctxs->m_allStates.m_refSbbCtxId[prev_state] >= 0) {
    memcpy(sbbFlags, cc->m_allSbbCtx[cc->m_prev_sbb_ctx_offset + ctxs->m_allStates.m_refSbbCtxId[prev_state]].sbbFlags, numSbb * sizeof(uint8_t));
    memcpy(levels + scan_pos, cc->m_allSbbCtx[cc->m_prev_sbb_ctx_offset + ctxs->m_allStates.m_refSbbCtxId[prev_state]].levels + scan_pos, setCpSize);
  }
  else {
    memset(sbbFlags, 0, numSbb * sizeof(uint8_t));
    memset(levels + scan_pos, 0, setCpSize);
  }
  sbbFlags[cg_pos] = !!ctxs->m_allStates.m_numSigSbb[curr_state];
  memcpy(levels + scan_pos, ctxs->m_allStates.m_absLevelsAndCtxInit[curr_state], 16 * sizeof(uint8_t));

  const int       sigNSbb = ((next_sbb_right ? sbbFlags[next_sbb_right] : false) || (next_sbb_below ? sbbFlags[next_sbb_below] : false) ? 1 : 0);
  ctxs->m_allStates.m_numSigSbb[curr_state] = 0;
  if (prev_state != -1) {
    ctxs->m_allStates.m_remRegBins[curr_state] = ctxs->m_allStates.m_remRegBins[prev_state];
  }
  else {
    int ctxBinSampleRatio = 28;
    // (scanInfo.chType == COLOR_Y) ? MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_LUMA : MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_CHROMA;
    ctxs->m_allStates.m_remRegBins[curr_state] = (ctxs->m_allStates.effWidth * ctxs->m_allStates.effHeight * ctxBinSampleRatio) / 16;
  }
  ctxs->m_allStates.m_goRicePar[curr_state] = 0;
  ctxs->m_allStates.m_refSbbCtxId[curr_state] = curr_state & 3;
  ctxs->m_allStates.m_sbbFracBits[curr_state][0] = cc->m_sbbFlagBits[sigNSbb][0];
  ctxs->m_allStates.m_sbbFracBits[curr_state][1] = cc->m_sbbFlagBits[sigNSbb][1];

  uint16_t *templateCtxInit = ctxs->m_allStates.m_absLevelsAndCtxInit[curr_state] + 8;
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
  memset(ctxs->m_allStates.m_absLevelsAndCtxInit[curr_state], 0, 16 * sizeof(uint8_t));
}

static INLINE void updateStateEOS(
  context_store*  ctxs,
  const uint32_t  scan_pos,
  const uint32_t  cg_pos,
  const uint32_t  sigCtxOffsetNext,
  const uint32_t  gtxCtxOffsetNext,
  const uint32_t  width_in_sbb,
  const uint32_t  height_in_sbb,
  const uint32_t  next_sbb_right,
  const uint32_t  next_sbb_below,
  const Decision* decisions,
  int             decision_id);

static void update_state_eos_avx2(context_store* ctxs, const uint32_t scan_pos, const uint32_t cg_pos,
                                  const uint32_t sigCtxOffsetNext, const uint32_t gtxCtxOffsetNext,
                                  const uint32_t width_in_sbb, const uint32_t height_in_sbb,
                                  const uint32_t next_sbb_right, const uint32_t next_sbb_below,
                                  const Decision* decisions)
{
  all_depquant_states* state = &ctxs->m_allStates;
  bool all_above_minus_two = true;
  bool all_between_zero_and_three = true;
  bool all_above_four = true;

  
  int state_offset = ctxs->m_curr_state_offset;
  __m256i rd_cost = _mm256_load_si256((__m256i const*)decisions->rdCost);
  _mm256_store_si256((__m256i *)& ctxs->m_allStates.m_rdCost[state_offset], rd_cost);
  for (int i = 0; i < 4; ++i) {
    all_above_minus_two &= decisions->prevId[i] > -2;
    all_between_zero_and_three &= decisions->prevId[i] >= 0 && decisions->prevId[i] < 4;
    all_above_four &= decisions->prevId[i] >= 4;
  }
  if (all_above_minus_two) {
    bool all_have_previous_state = true;
    __m128i prev_state;
    __m128i prev_state_no_offset;
    __m128i abs_level = _mm_load_si128((const __m128i*)decisions->absLevel);
    if (all_above_four) {
      prev_state = _mm_set1_epi32(ctxs->m_skip_state_offset);
      prev_state_no_offset = _mm_sub_epi32(_mm_load_si128((const __m128i*)decisions->prevId), _mm_set1_epi32(4));
      prev_state = _mm_add_epi32(
        prev_state,
            prev_state_no_offset
      );
      memset(&state->m_numSigSbb[state_offset], 0, 4);
      for (int i = 0; i < 4; ++i) {
        memset(state->m_absLevelsAndCtxInit[state_offset + i], 0, 16 * sizeof(uint8_t));    
      }
    } else if (all_between_zero_and_three) {
      prev_state_no_offset = _mm_set1_epi32(ctxs->m_prev_state_offset);
      prev_state = _mm_add_epi32(
        prev_state_no_offset,
        _mm_load_si128((const __m128i*)decisions->prevId)
      );
      __m128i control = _mm_setr_epi8(0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
      __m128i prev_state_with_ff_high_bytes = _mm_or_epi32(prev_state, _mm_set1_epi32(0xffffff00));
      __m128i num_sig_sbb = _mm_load_si128((const __m128i*)state->m_numSigSbb);
      num_sig_sbb = _mm_shuffle_epi8(num_sig_sbb, prev_state_with_ff_high_bytes);
      num_sig_sbb = _mm_add_epi32(
        num_sig_sbb,
        _mm_min_epi32(abs_level, _mm_set1_epi32(1))
      );

      num_sig_sbb = _mm_shuffle_epi8(num_sig_sbb, control);
      int num_sig_sbb_s = _mm_extract_epi32(num_sig_sbb, 0);
      memcpy(&state->m_numSigSbb[state_offset], &num_sig_sbb_s, 4);

      int32_t prev_state_scalar[4];
      _mm_storeu_si128((__m128i*)prev_state_scalar, prev_state);
      for (int i = 0; i < 4; ++i) {
        memcpy(state->m_absLevelsAndCtxInit[state_offset + i], state->m_absLevelsAndCtxInit[prev_state_scalar[i]], 16 * sizeof(uint8_t));
      }
    } else {
      int prev_state_s[4] = {-1, -1, -1, -1};
      for (int i = 0; i < 4; ++i) {
        const int decision_id = i;
        const int curr_state_offset = state_offset + i;
        if (decisions->prevId[decision_id] >= 4) {
          prev_state_s[i] = ctxs->m_skip_state_offset + (decisions->prevId[decision_id] - 4);
          state->m_numSigSbb[curr_state_offset] = 0;
          memset(state->m_absLevelsAndCtxInit[curr_state_offset], 0, 16 * sizeof(uint8_t));
        } else if (decisions->prevId[decision_id] >= 0) {
          prev_state_s[i] = ctxs->m_prev_state_offset + decisions->prevId[decision_id];
          state->m_numSigSbb[curr_state_offset] = state->m_numSigSbb[prev_state_s[i]] + !!decisions->absLevel[decision_id];
          memcpy(state->m_absLevelsAndCtxInit[curr_state_offset], state->m_absLevelsAndCtxInit[prev_state_s[i]], 16 * sizeof(uint8_t));
        } else {
          state->m_numSigSbb[curr_state_offset] = 1;
          memset(state->m_absLevelsAndCtxInit[curr_state_offset], 0, 16 * sizeof(uint8_t));
          all_have_previous_state = false;
        }
      }
      prev_state = _mm_loadu_si128((__m128i const*)prev_state_s);
    }
    uint32_t level_offset = scan_pos & 15;
    __m128i  max_abs = _mm_min_epi32(abs_level, _mm_set1_epi32(32));
    uint32_t max_abs_s[4];
    _mm_storeu_si128((__m128i*)max_abs_s, max_abs);
    for (int i = 0; i < 4; ++i) {
      uint8_t* levels = (uint8_t*)state->m_absLevelsAndCtxInit[state_offset + i];
      levels[level_offset] = max_abs_s[i];
    }

    // Update common context
    __m128i last;
    {
      const uint32_t numSbb = width_in_sbb * height_in_sbb;
      common_context* cc = &ctxs->m_common_context;
      size_t         setCpSize = cc->m_nbInfo[scan_pos - 1].maxDist * sizeof(uint8_t);
      int previous_state_array[4];
      _mm_storeu_si128((__m128i*)previous_state_array, prev_state);
      for (int curr_state = 0; curr_state < 4; ++curr_state) {
        uint8_t* sbbFlags = cc->m_allSbbCtx[cc->m_curr_sbb_ctx_offset + (curr_state)].sbbFlags;
        uint8_t* levels = cc->m_allSbbCtx[cc->m_curr_sbb_ctx_offset + (curr_state)].levels;
        const int p_state = previous_state_array[curr_state];
        if (p_state != -1 && ctxs->m_allStates.m_refSbbCtxId[p_state] >= 0) {
          const int prev_sbb = cc->m_prev_sbb_ctx_offset + ctxs->m_allStates.m_refSbbCtxId[p_state];
          memcpy(sbbFlags, cc->m_allSbbCtx[prev_sbb].sbbFlags, numSbb * sizeof(uint8_t));
          memcpy(levels + scan_pos, cc->m_allSbbCtx[prev_sbb].levels + scan_pos, setCpSize);
        } else {
          memset(sbbFlags, 0, numSbb * sizeof(uint8_t));
          memset(levels + scan_pos, 0, setCpSize);
        }
        sbbFlags[cg_pos] = !!ctxs->m_allStates.m_numSigSbb[curr_state + state_offset];
        memcpy(levels + scan_pos, ctxs->m_allStates.m_absLevelsAndCtxInit[curr_state + state_offset], 16 * sizeof(uint8_t));
      }

      __m128i sbb_offsets = _mm_set_epi32(3 * numSbb, 2 * numSbb, 1 * numSbb, 0);
      __m128i next_sbb_right_m = _mm_set1_epi32(next_sbb_right);
      __m128i sbb_offsets_right = _mm_add_epi32(sbb_offsets, next_sbb_right_m);
      __m128i sbb_right = next_sbb_right ? _mm_i32gather_epi32((const int *)cc->m_allSbbCtx[cc->m_curr_sbb_ctx_offset].sbbFlags, sbb_offsets_right, 1) : _mm_set1_epi32(0);

      __m128i sbb_offsets_below = _mm_add_epi32(sbb_offsets, _mm_set1_epi32(next_sbb_below));
      __m128i sbb_below = next_sbb_below ? _mm_i32gather_epi32((const int *)cc->m_allSbbCtx[cc->m_curr_sbb_ctx_offset].sbbFlags, sbb_offsets_below, 1) : _mm_set1_epi32(0);

      __m128i sig_sbb = _mm_or_epi32(sbb_right, sbb_below);
      sig_sbb         = _mm_and_si128(sig_sbb, _mm_set1_epi32(0xff));
      sig_sbb = _mm_min_epi32(sig_sbb, _mm_set1_epi32(1));
      __m256i sbb_frac_bits = _mm256_i32gather_epi64((int64_t *)cc->m_sbbFlagBits[0], sig_sbb, 8);
      _mm256_store_si256((__m256i*)state->m_sbbFracBits[state_offset], sbb_frac_bits);

      memset(&state->m_numSigSbb[state_offset], 0, 4);
      memset(&state->m_goRicePar[state_offset], 0, 4);

      uint8_t states[4] = {0, 1, 2, 3};
      memcpy(&state->m_refSbbCtxId[state_offset], states, 4);
      if (all_have_previous_state) {
        __m128i rem_reg_bins = _mm_i32gather_epi32(state->m_remRegBins, prev_state, 4);
        _mm_store_si128((__m128i*) & state->m_remRegBins[state_offset], rem_reg_bins);
      } else {
        const int temp = (state->effWidth * state->effHeight * 28) / 16;
        for (int i = 0; i < 4; ++i) {
          if (previous_state_array[i] != -1) {
            state->m_remRegBins[i + state_offset] = state->m_remRegBins[previous_state_array[i]];
          } else {
            state->m_remRegBins[i + state_offset] = temp;
          }
        }
      }
      
      const int        scanBeg = scan_pos - 16;
      const NbInfoOut* nbOut = cc->m_nbInfo + scanBeg;
      const uint8_t*   absLevels = cc->m_allSbbCtx[cc->m_curr_sbb_ctx_offset].levels + scanBeg;

      __m128i          levels_offsets = _mm_set_epi32(cc->num_coeff * 3, cc->num_coeff * 2, cc->num_coeff * 1, 0);
      __m128i          first_byte = _mm_set1_epi32(0xff);
      __m128i          ones = _mm_set1_epi32(1);
      __m128i         fours = _mm_set1_epi32(4);
      __m256i          all[4];
      uint64_t         temp[4];
      const __m256i v_shuffle = _mm256_set_epi8(15, 14,  7,  6, 13, 12,  5,  4, 11, 10,  3,  2,  9,  8,  1,  0,
                                                31, 30, 23, 22, 29, 28, 21, 20, 27, 26, 19, 18, 25, 24, 17, 16);

      for (int id = 0; id < 16; id++, nbOut++) {
        if (nbOut->num == 0) {
          temp[id % 4] = 0;
          if (id % 4 == 3) {
            all[id / 4] = _mm256_loadu_si256((__m256i const*)temp);
            all[id / 4] = _mm256_shuffle_epi8(all[id / 4], v_shuffle);
          }
          continue;
        }
        __m128i sum_abs = _mm_set1_epi32(0);
        __m128i sum_abs_1 = _mm_set1_epi32(0);
        __m128i sum_num = _mm_set1_epi32(0);
        switch (nbOut->num) {
        case 5:
          {
            __m128i offset = _mm_add_epi32(levels_offsets, _mm_set1_epi32(nbOut->outPos[4]));
            __m128i t = _mm_i32gather_epi32((const int *)absLevels, offset, 1);
            t = _mm_and_si128(t, first_byte);
            sum_abs = _mm_add_epi32(sum_abs, t);
            sum_num   = _mm_add_epi32(sum_num, _mm_min_epi32(t, ones));
            __m128i min_t = _mm_min_epi32(
              t,
              _mm_add_epi32(
                fours,
                _mm_and_si128(t, ones)
              )
            );
            sum_abs_1 = _mm_add_epi32(sum_abs_1, min_t);
          }
        case 4: {
            __m128i offset = _mm_add_epi32(levels_offsets, _mm_set1_epi32(nbOut->outPos[3]));
            __m128i t     = _mm_i32gather_epi32((const int*)absLevels, offset, 1);
            t = _mm_and_si128(t, first_byte);
            sum_abs = _mm_add_epi32(sum_abs, t);
            sum_num   = _mm_add_epi32(sum_num, _mm_min_epi32(t, ones));
            __m128i min_t = _mm_min_epi32(
              t,
              _mm_add_epi32(
                fours,
                _mm_and_si128(t, ones)));
            sum_abs_1 = _mm_add_epi32(sum_abs_1, min_t);
        }
        case 3: {
            __m128i offset = _mm_add_epi32(levels_offsets, _mm_set1_epi32(nbOut->outPos[2]));
            __m128i t     = _mm_i32gather_epi32((const int*)absLevels, offset, 1);
            t = _mm_and_si128(t, first_byte);
            sum_abs = _mm_add_epi32(sum_abs, t);
            sum_num   = _mm_add_epi32(sum_num, _mm_min_epi32(t, ones));
            __m128i min_t = _mm_min_epi32(
              t,
              _mm_add_epi32(
                fours,
                _mm_and_si128(t, ones)));
            sum_abs_1 = _mm_add_epi32(sum_abs_1, min_t);
        }
        case 2: {
            __m128i offset = _mm_add_epi32(levels_offsets, _mm_set1_epi32(nbOut->outPos[1]));
            __m128i t     = _mm_i32gather_epi32((const int*)absLevels, offset, 1);
            t = _mm_and_si128(t, first_byte);
            sum_abs = _mm_add_epi32(sum_abs, t);
            sum_num   = _mm_add_epi32(sum_num, _mm_min_epi32(t, ones));
            __m128i min_t = _mm_min_epi32(
              t,
              _mm_add_epi32(
                fours,
                _mm_and_si128(t, ones)));
            sum_abs_1 = _mm_add_epi32(sum_abs_1, min_t);
        }
        case 1: {
            __m128i offset = _mm_add_epi32(levels_offsets, _mm_set1_epi32(nbOut->outPos[0]));
            __m128i t = _mm_i32gather_epi32((const int*)absLevels, offset, 1);
            t = _mm_and_si128(t, first_byte);
            sum_abs = _mm_add_epi32(sum_abs, t);
            sum_num = _mm_add_epi32(sum_num, _mm_min_epi32(t, ones));
            __m128i min_t = _mm_min_epi32(
              t,
              _mm_add_epi32(
                fours,
                _mm_and_si128(t, ones)));
            sum_abs_1 = _mm_add_epi32(sum_abs_1, min_t);
        }
            break;
        default:
          assert(0);
        }
        sum_abs_1 = _mm_slli_epi32(sum_abs_1, 3);
        sum_abs = _mm_slli_epi32(_mm_min_epi32(_mm_set1_epi32(127), sum_abs), 8);
        __m128i template_ctx_init = _mm_add_epi32(sum_num, sum_abs);
        template_ctx_init = _mm_add_epi32(template_ctx_init, sum_abs_1);
        __m128i shuffle_mask = _mm_setr_epi8(0, 1, 4, 5, 8, 9, 12, 13, 0, 0, 0, 0, 0, 0, 0, 0);
        __m128i shuffled_template_ctx_init = _mm_shuffle_epi8(template_ctx_init, shuffle_mask);
        temp[id % 4] = _mm_extract_epi64(shuffled_template_ctx_init, 0);
        if (id % 4 == 3) {
          all[id / 4] = _mm256_loadu_si256((__m256i const*)temp);
          all[id / 4] = _mm256_shuffle_epi8(all[id / 4], v_shuffle);
          last = template_ctx_init;
        }
      }

      __m256i* v_src_tmp = all;

      __m256i v_tmp[4];
      v_tmp[0] = _mm256_permute2x128_si256(v_src_tmp[0], v_src_tmp[1], 0x20);
      v_tmp[1] = _mm256_permute2x128_si256(v_src_tmp[0], v_src_tmp[1], 0x31);
      v_tmp[2] = _mm256_permute2x128_si256(v_src_tmp[2], v_src_tmp[3], 0x20);
      v_tmp[3] = _mm256_permute2x128_si256(v_src_tmp[2], v_src_tmp[3], 0x31);

      __m256i v_tmp16_lo[2];
      __m256i v_tmp16_hi[2];
      v_tmp16_lo[0] = _mm256_unpacklo_epi32(v_tmp[0], v_tmp[1]);
      v_tmp16_lo[1] = _mm256_unpacklo_epi32(v_tmp[2], v_tmp[3]);
      v_tmp16_hi[0] = _mm256_unpackhi_epi32(v_tmp[0], v_tmp[1]);
      v_tmp16_hi[1] = _mm256_unpackhi_epi32(v_tmp[2], v_tmp[3]);

      v_tmp[0] = _mm256_permute4x64_epi64(v_tmp16_lo[0], _MM_SHUFFLE(3, 1, 2, 0));
      v_tmp[1] = _mm256_permute4x64_epi64(v_tmp16_lo[1], _MM_SHUFFLE(3, 1, 2, 0));
      v_tmp[2] = _mm256_permute4x64_epi64(v_tmp16_hi[0], _MM_SHUFFLE(3, 1, 2, 0));
      v_tmp[3] = _mm256_permute4x64_epi64(v_tmp16_hi[1], _MM_SHUFFLE(3, 1, 2, 0));

      _mm256_store_si256((__m256i*)(state->m_absLevelsAndCtxInit[state_offset] + 8),  _mm256_permute2x128_si256(v_tmp[0], v_tmp[1], 0x20));
      _mm256_store_si256((__m256i*)(state->m_absLevelsAndCtxInit[state_offset + 1] + 8),  _mm256_permute2x128_si256(v_tmp[0], v_tmp[1], 0x31));
      _mm256_store_si256((__m256i*)(state->m_absLevelsAndCtxInit[state_offset + 2] + 8),  _mm256_permute2x128_si256(v_tmp[2], v_tmp[3], 0x20));
      _mm256_store_si256((__m256i*)(state->m_absLevelsAndCtxInit[state_offset + 3] + 8),  _mm256_permute2x128_si256(v_tmp[2], v_tmp[3], 0x31));

      for (int i = 0; i < 4; ++i) {
        memset(state->m_absLevelsAndCtxInit[state_offset + i], 0, 16);
      }
    }

    __m128i sum_num = _mm_and_si128(last, _mm_set1_epi32(7));
    __m128i sum_abs1 = _mm_and_si128(
      _mm_srli_epi32(last, 3),
      _mm_set1_epi32(31));

    __m128i sum_abs_min = _mm_min_epi32(
      _mm_set1_epi32(3),
      _mm_srli_epi32(
        _mm_add_epi32(sum_abs1, _mm_set1_epi32(1)),
        1));

    __m128i offsets = _mm_set_epi32(12 * 3, 12 * 2, 12 * 1, 12 * 0);
    offsets = _mm_add_epi32(offsets, _mm_set1_epi32(sigCtxOffsetNext));
    offsets         = _mm_add_epi32(offsets, sum_abs_min);
    __m256i sig_frac_bits = _mm256_i32gather_epi64((const int64_t *)&state->m_sigFracBitsArray[state_offset][0][0], offsets, 8);
    _mm256_store_si256((__m256i*)&state->m_sigFracBits[state_offset][0], sig_frac_bits);


    __m128i sum_gt1 = _mm_sub_epi32(sum_abs1, sum_num);
    __m128i min_gt1 = _mm_min_epi32(sum_gt1, _mm_set1_epi32(4));
    uint32_t sum_gt1_s[4];
    _mm_storeu_si128((__m128i*)sum_gt1_s, min_gt1);
    for (int i = 0; i < 4; ++i) {
      memcpy(state->m_coeffFracBits[state_offset + i], state->m_gtxFracBitsArray[sum_gt1_s[i] + gtxCtxOffsetNext], sizeof(state->m_coeffFracBits[0]));
    }
  }
  else {
    for (int i = 0; i < 4; i++) {
      updateStateEOS(
        ctxs,
        scan_pos,
        cg_pos,
        sigCtxOffsetNext,
        gtxCtxOffsetNext,
        width_in_sbb,
        height_in_sbb,
        next_sbb_right,
        next_sbb_below,
        decisions,
        i);
    }
  }
}


static INLINE void updateStateEOS(
  context_store*   ctxs,
  const uint32_t   scan_pos,
  const uint32_t   cg_pos,
  const uint32_t   sigCtxOffsetNext,
  const uint32_t   gtxCtxOffsetNext,
  const uint32_t   width_in_sbb,
  const uint32_t   height_in_sbb,
  const uint32_t   next_sbb_right,
  const uint32_t   next_sbb_below,
  const Decision * decisions,
  int              decision_id)
{
  all_depquant_states* state = &ctxs->m_allStates;
  int curr_state_offset = ctxs->m_curr_state_offset + decision_id;
  state->m_rdCost[curr_state_offset] = decisions->rdCost[decision_id];
  if (decisions->prevId[decision_id] > -2) {
    int prvState = -1;
    if (decisions->prevId[decision_id] >= 4) {
      prvState = ctxs->m_skip_state_offset + (decisions->prevId[decision_id] - 4);
      state->m_numSigSbb[curr_state_offset] = 0;
      memset(state->m_absLevelsAndCtxInit[curr_state_offset], 0, 16 * sizeof(uint8_t));
    }
    else if (decisions->prevId[decision_id] >= 0) {
      prvState = ctxs->m_prev_state_offset + decisions->prevId[decision_id];
      state->m_numSigSbb[curr_state_offset] = state->m_numSigSbb[prvState] + !!decisions->absLevel[decision_id];
      memcpy(state->m_absLevelsAndCtxInit[curr_state_offset], state->m_absLevelsAndCtxInit[prvState], 16 * sizeof(uint8_t));
    }
    else {
      state->m_numSigSbb[curr_state_offset] = 1;
      memset(state->m_absLevelsAndCtxInit[curr_state_offset], 0, 16 * sizeof(uint8_t));
    }
    uint8_t* temp = (uint8_t*)(&state->m_absLevelsAndCtxInit[curr_state_offset][scan_pos & 15]);
    *temp = (uint8_t)MIN(32, decisions->absLevel[decision_id]);

    update_common_context(ctxs, state->m_commonCtx, scan_pos, cg_pos, width_in_sbb, height_in_sbb, next_sbb_right,
                          next_sbb_below, prvState, ctxs->m_curr_state_offset + decision_id);

    coeff_t tinit = state->m_absLevelsAndCtxInit[curr_state_offset][8 + ((scan_pos - 1) & 15)];
    coeff_t sumNum = tinit & 7;
    coeff_t sumAbs1 = (tinit >> 3) & 31;
    coeff_t sumGt1 = sumAbs1 - sumNum;
    state->m_sigFracBits[curr_state_offset][0] = state->m_sigFracBitsArray[curr_state_offset][sigCtxOffsetNext + MIN((sumAbs1 + 1) >> 1, 3)][0];
    state->m_sigFracBits[curr_state_offset][1] = state->m_sigFracBitsArray[curr_state_offset][sigCtxOffsetNext + MIN((sumAbs1 + 1) >> 1, 3)][1];

    memcpy(state->m_coeffFracBits[curr_state_offset],
           state->m_gtxFracBitsArray[gtxCtxOffsetNext + (sumGt1 < 4 ? sumGt1 : 4)], sizeof(state->m_coeffFracBits[0]));
  }
}
static INLINE void updateState(
  context_store*  ctxs,
  int             numIPos,
  const uint32_t  scan_pos,
  const Decision* decisions,
  const uint32_t  sigCtxOffsetNext,
  const uint32_t  gtxCtxOffsetNext,
  const NbInfoSbb next_nb_info_ssb,
  const int       baseLevel,
  const bool      extRiceFlag,
  int             decision_id);

static INLINE void update_states_avx2(
  context_store*  ctxs,
  int             numIPos,
  const uint32_t  scan_pos,
  const Decision* decisions,
  const uint32_t  sigCtxOffsetNext,
  const uint32_t  gtxCtxOffsetNext,
  const NbInfoSbb next_nb_info_ssb,
  const int       baseLevel,
  const bool      extRiceFlag)
{
  all_depquant_states* state = &ctxs->m_allStates;

  bool all_non_negative = true;
  bool all_above_minus_two = true;
  bool all_minus_one = true;
  for (int i = 0; i < 4; ++i) {
    all_non_negative &= decisions->prevId[i] >= 0;
    all_above_minus_two &= decisions->prevId[i] > -2;
    all_minus_one &= decisions->prevId[i] == -1;
  }
  int state_offset = ctxs->m_curr_state_offset;
  __m256i rd_cost = _mm256_load_si256((__m256i const*)decisions->rdCost);
  _mm256_store_si256((__m256i *)& ctxs->m_allStates.m_rdCost[state_offset], rd_cost);
  if (all_above_minus_two) {

    bool    rem_reg_all_gte_4 = true;
    bool    rem_reg_all_lt4 = true;

    __m128i abs_level = _mm_load_si128((__m128i const*)decisions->absLevel);
    if (all_non_negative) {
      __m128i prv_states  = _mm_load_si128((__m128i const*)decisions->prevId);
      __m128i prev_offset = _mm_set1_epi32(ctxs->m_prev_state_offset);
      prv_states = _mm_add_epi32(prv_states, prev_offset);
      __m128i control = _mm_setr_epi8(0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
      __m128i shuffled_prev_states = _mm_shuffle_epi8(prv_states, control);
      
      __m128i sig_sbb   = _mm_load_si128((__m128i const*)state->m_numSigSbb);
      sig_sbb = _mm_shuffle_epi8(sig_sbb, shuffled_prev_states);
      __m128i has_coeff = _mm_min_epi32(abs_level, _mm_set1_epi32(1));
      has_coeff         = _mm_shuffle_epi8(has_coeff, control);
      sig_sbb           = _mm_or_epi32(sig_sbb, has_coeff);
      int sig_sbb_i = _mm_extract_epi32(sig_sbb, 0);
      memcpy(&state->m_numSigSbb[state_offset], &sig_sbb_i, 4);
      
      __m128i ref_sbb_ctx_idx = _mm_load_si128((__m128i const*)state->m_refSbbCtxId);
      ref_sbb_ctx_idx = _mm_shuffle_epi8(ref_sbb_ctx_idx, shuffled_prev_states);
      int ref_sbb_ctx = _mm_extract_epi32(ref_sbb_ctx_idx, 0);
      memcpy(&state->m_refSbbCtxId[state_offset], &ref_sbb_ctx, 4);
      
      __m128i go_rice_par = _mm_load_si128((__m128i const*)state->m_goRicePar);
      go_rice_par = _mm_shuffle_epi8(go_rice_par, shuffled_prev_states);
      int go_rice_par_i = _mm_extract_epi32(go_rice_par, 0);
      memcpy(&state->m_goRicePar[state_offset], &go_rice_par_i, 4);

      
      __m256i sbb_frac_bits = _mm256_i32gather_epi64((const int64_t *)state->m_sbbFracBits[0], prv_states, 8);
      _mm256_store_si256((__m256i*)&state->m_sbbFracBits[state_offset][0], sbb_frac_bits);

      __m128i rem_reg_bins = _mm_i32gather_epi32(state->m_remRegBins, prv_states, 4);
      __m128i ones = _mm_set1_epi32(1);
      rem_reg_bins = _mm_sub_epi32(rem_reg_bins, ones);

      __m128i reg_bins_sub = _mm_set1_epi32(0);
      __m128i abs_level_smaller_than_two = _mm_cmplt_epi32(abs_level, _mm_set1_epi32(2));
      __m128i secondary = _mm_blendv_epi8(_mm_set1_epi32(3), abs_level, abs_level_smaller_than_two);

      __m128i rem_reg_bins_smaller_than_four = _mm_cmplt_epi32(rem_reg_bins, _mm_set1_epi32(4));
      reg_bins_sub = _mm_blendv_epi8(secondary, reg_bins_sub, rem_reg_bins_smaller_than_four);
      rem_reg_bins = _mm_sub_epi32(rem_reg_bins, reg_bins_sub);
      _mm_store_si128((__m128i*)&state->m_remRegBins[state_offset], rem_reg_bins);

      __m128i mask = _mm_cmpgt_epi32(rem_reg_bins, _mm_set1_epi32(3)); 
      int     bit_mask = _mm_movemask_epi8(mask);           
      rem_reg_all_gte_4 = (bit_mask == 0xFFFF);
      mask = _mm_cmplt_epi32(rem_reg_bins, _mm_set1_epi32(4));
      bit_mask = _mm_movemask_epi8(mask); 
      rem_reg_all_lt4 = (bit_mask == 0xFFFF);

      int32_t prv_states_scalar[4];
      _mm_storeu_si128((__m128i*)prv_states_scalar, prv_states);
      for (int i = 0; i < 4; ++i) {
        memcpy(state->m_absLevelsAndCtxInit[state_offset + i], state->m_absLevelsAndCtxInit[prv_states_scalar[i]], 48 * sizeof(uint8_t));        
      }
    }
    else if (all_minus_one) {
      memset(&state->m_numSigSbb[state_offset], 1, 4);
      memset(&state->m_refSbbCtxId[state_offset], -1, 4);

      const int a = (state->effWidth * state->effHeight * 28) / 16;

      __m128i   rem_reg_bins = _mm_set1_epi32(a);
      __m128i   sub = _mm_blendv_epi8(
        _mm_set1_epi32(3),
        abs_level,
        _mm_cmplt_epi32(abs_level, _mm_set1_epi32(2))
      );
      rem_reg_bins = _mm_sub_epi32(rem_reg_bins, sub);
      _mm_store_si128((__m128i*) & state->m_remRegBins[state_offset], rem_reg_bins);

      __m128i mask = _mm_cmpgt_epi32(rem_reg_bins, _mm_set1_epi32(3));
      int     bit_mask = _mm_movemask_epi8(mask);
      rem_reg_all_gte_4 = (bit_mask == 0xFFFF);
      mask = _mm_cmplt_epi32(rem_reg_bins, _mm_set1_epi32(4));
      bit_mask = _mm_movemask_epi8(mask);
      rem_reg_all_lt4 = (bit_mask == 0xFFFF);
      
      memset(state->m_absLevelsAndCtxInit[state_offset], 0, 48 * sizeof(uint8_t) * 4);
      
    }
    else {
      for (int i = 0; i< 4; ++i) {
        const int decision_id = i;
        const int state_id = state_offset + i;
        if (decisions->prevId[decision_id] >= 0) {
          const int prvState = ctxs->m_prev_state_offset + decisions->prevId[decision_id];
          state->m_numSigSbb[state_id] = (state->m_numSigSbb[prvState]) || !!decisions->absLevel[decision_id];
          state->m_refSbbCtxId[state_id] = state->m_refSbbCtxId[prvState];
          state->m_sbbFracBits[state_id][0] = state->m_sbbFracBits[prvState][0];
          state->m_sbbFracBits[state_id][1] = state->m_sbbFracBits[prvState][1];
          state->m_remRegBins[state_id] = state->m_remRegBins[prvState] - 1;
          state->m_goRicePar[state_id] = state->m_goRicePar[prvState];
          if (state->m_remRegBins[state_id] >= 4) {
            state->m_remRegBins[state_id] -= (decisions->absLevel[decision_id] < 2 ? (unsigned)decisions->absLevel[decision_id] : 3);
          }
          memcpy(state->m_absLevelsAndCtxInit[state_id], state->m_absLevelsAndCtxInit[prvState], 48 * sizeof(uint8_t));
        } else {
          state->m_numSigSbb[state_id] = 1;
          state->m_refSbbCtxId[state_id] = -1;
          int ctxBinSampleRatio = 28;
          //(scanInfo.chType == CHANNEL_TYPE_LUMA) ? MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_LUMA : MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_CHROMA;
          state->m_remRegBins[state_id] = (state->effWidth * state->effHeight * ctxBinSampleRatio) / 16 - (decisions->absLevel[decision_id] < 2 ? (unsigned)decisions->absLevel[decision_id] : 3);
          memset(state->m_absLevelsAndCtxInit[state_id], 0, 48 * sizeof(uint8_t));
        }
        rem_reg_all_gte_4 &= state->m_remRegBins[state_id] >= 4;
        rem_reg_all_lt4 &= state->m_remRegBins[state_id] < 4;
      }
    }
    uint32_t level_offset = scan_pos & 15;
    __m128i   max_abs = _mm_min_epi32(abs_level, _mm_set1_epi32(32));
    uint32_t max_abs_s[4];
    _mm_storeu_si128((__m128i*)max_abs_s, max_abs);
    for (int i = 0; i < 4; ++i) {
      uint8_t* levels = (uint8_t*)state->m_absLevelsAndCtxInit[state_offset + i];
      levels[level_offset] = max_abs_s[i];
    }
    state->all_gte_four = rem_reg_all_gte_4;
    state->all_lt_four = rem_reg_all_lt4;
    if (rem_reg_all_gte_4) {
      const __m128i  first_two_bytes = _mm_set1_epi32(0xffff);
      const __m128i  first_byte = _mm_set1_epi32(0xff);
      const __m128i  ones = _mm_set1_epi32(1);
      const uint32_t tinit_offset = MIN(level_offset - 1u, 15u) + 8;
      const __m128i levels_start_offsets = _mm_set_epi32(48 * 3, 48 * 2, 48 * 1, 48 * 0);
      const __m128i ctx_start_offsets = _mm_srli_epi32(levels_start_offsets, 1);
      __m128i        tinit = _mm_i32gather_epi32(
        (int *)state->m_absLevelsAndCtxInit[state_offset],
        _mm_add_epi32(ctx_start_offsets, _mm_set1_epi32(tinit_offset)),
        2);
      tinit = _mm_and_si128(tinit, first_two_bytes);
      __m128i sum_abs1 = _mm_and_si128(_mm_srli_epi32(tinit, 3), _mm_set1_epi32(31));
      __m128i sum_num = _mm_and_si128(tinit, _mm_set1_epi32(7));

      uint8_t* levels = (uint8_t*)state->m_absLevelsAndCtxInit[state_offset];
      switch (numIPos) {
      case 5:
        {
          __m128i t = _mm_i32gather_epi32(
            (int *)levels,
            _mm_add_epi32(levels_start_offsets, _mm_set1_epi32(next_nb_info_ssb.inPos[4])),
            1);
          t = _mm_and_si128(t, first_byte);
          __m128i min_arg = _mm_min_epi32(
            _mm_add_epi32(_mm_set1_epi32(4), _mm_and_si128(t, ones)),
            t
          );
          sum_abs1 = _mm_add_epi32(
            sum_abs1,
            min_arg
          );
          sum_num = _mm_add_epi32(
            sum_num,
            _mm_min_epi32(_mm_and_si128(t, first_byte), ones));
        }
      case 4:
        {
          __m128i t = _mm_i32gather_epi32(
            (int*)levels,
            _mm_add_epi32(levels_start_offsets, _mm_set1_epi32(next_nb_info_ssb.inPos[3])),
            1);
          t = _mm_and_si128(t, first_byte);
          __m128i min_arg = _mm_min_epi32(
            _mm_add_epi32(_mm_set1_epi32(4), _mm_and_si128(t, ones)),
            t
          );
          sum_abs1 = _mm_add_epi32(
            sum_abs1,
            min_arg
          );
          sum_num = _mm_add_epi32(
            sum_num,
            _mm_min_epi32(_mm_and_si128(t, first_byte), ones));
        }
      case 3:
        {
          __m128i t = _mm_i32gather_epi32(
            (int*)levels,
            _mm_add_epi32(levels_start_offsets, _mm_set1_epi32(next_nb_info_ssb.inPos[2])),
            1);
          t = _mm_and_si128(t, first_byte);
          __m128i min_arg = _mm_min_epi32(
            _mm_add_epi32(_mm_set1_epi32(4), _mm_and_si128(t, ones)),
            t
          );
          sum_abs1 = _mm_add_epi32(
            sum_abs1,
            min_arg
          );
          sum_num = _mm_add_epi32(
            sum_num,
            _mm_min_epi32(_mm_and_si128(t, first_byte), ones));
        }
      case 2:
        {
          __m128i t = _mm_i32gather_epi32(
            (int*)levels,
            _mm_add_epi32(levels_start_offsets, _mm_set1_epi32(next_nb_info_ssb.inPos[1])),
            1);
          t = _mm_and_si128(t, first_byte);
        __m128i min_arg = _mm_min_epi32(
              _mm_add_epi32(_mm_set1_epi32(4), _mm_and_si128(t, ones)),
              t
            );
          sum_abs1 = _mm_add_epi32(
            sum_abs1,
            min_arg
          );
          sum_num = _mm_add_epi32(
            sum_num,
            _mm_min_epi32(_mm_and_si128(t, first_byte), ones));
        }
      case 1: {
          __m128i t = _mm_i32gather_epi32(
            (int*)levels,
            _mm_add_epi32(levels_start_offsets, _mm_set1_epi32(next_nb_info_ssb.inPos[0])),
            1);
          t = _mm_and_si128(t, first_byte);
          __m128i min_arg = _mm_min_epi32(
            _mm_add_epi32(_mm_set1_epi32(4), _mm_and_si128(t, ones)),
            t
          );
          sum_abs1 = _mm_add_epi32(
            sum_abs1,
            min_arg
            );
          sum_num = _mm_add_epi32(
            sum_num,
            _mm_min_epi32(_mm_and_si128(t, first_byte), ones));
        } break;
      default:
          assert(0);
      }
      __m128i sum_gt1 = _mm_sub_epi32(sum_abs1, sum_num);
      __m128i  offsets = _mm_set_epi32(12 * 3, 12 * 2, 12 * 1, 12 * 0);
      offsets = _mm_add_epi32(offsets, _mm_set1_epi32(sigCtxOffsetNext));
      __m128i temp = _mm_min_epi32(
        _mm_srli_epi32(_mm_add_epi32(sum_abs1, ones), 1),
        _mm_set1_epi32(3));
      offsets = _mm_add_epi32(offsets, temp);
      __m256i sig_frac_bits = _mm256_i32gather_epi64((const int64_t *)state->m_sigFracBitsArray[state_offset][0], offsets, 8);
      _mm256_store_si256((__m256i*)&state->m_sigFracBits[state_offset][0], sig_frac_bits);

      sum_gt1 = _mm_min_epi32(sum_gt1, _mm_set1_epi32(4));
      sum_gt1 = _mm_add_epi32(sum_gt1, _mm_set1_epi32(gtxCtxOffsetNext));
      uint32_t sum_gt1_s[4];
      _mm_storeu_si128((__m128i*)sum_gt1_s, sum_gt1);
      for (int i = 0; i < 4; ++i) {
        memcpy(state->m_coeffFracBits[state_offset + i], state->m_gtxFracBitsArray[sum_gt1_s[i]], sizeof(state->m_coeffFracBits[0]));
      }

      __m128i sum_abs = _mm_srli_epi32(tinit, 8);
      sum_abs = _mm_min_epi32(sum_abs, _mm_set1_epi32(32));
      switch (numIPos) {
        case 5:
          {
          __m128i t = _mm_i32gather_epi32(
            (int*)levels,
            _mm_add_epi32(levels_start_offsets, _mm_set1_epi32(next_nb_info_ssb.inPos[4])),
            1);
          sum_abs = _mm_add_epi32(t, sum_abs);
          }
        case 4:
          {
          __m128i t = _mm_i32gather_epi32(
            (int*)levels,
            _mm_add_epi32(levels_start_offsets, _mm_set1_epi32(next_nb_info_ssb.inPos[3])),
            1);
          sum_abs = _mm_add_epi32(t, sum_abs);
          }
        case 3:
          {
          __m128i t = _mm_i32gather_epi32(
            (int*)levels,
            _mm_add_epi32(levels_start_offsets, _mm_set1_epi32(next_nb_info_ssb.inPos[2])),
            1);
          sum_abs = _mm_add_epi32(t, sum_abs);
          }
        case 2:
          {
          __m128i t = _mm_i32gather_epi32(
            (int*)levels,
            _mm_add_epi32(levels_start_offsets, _mm_set1_epi32(next_nb_info_ssb.inPos[1])),
            1);
          sum_abs = _mm_add_epi32(t, sum_abs);
          }
        case 1:
          {
          __m128i t = _mm_i32gather_epi32(
            (int*)levels,
            _mm_add_epi32(levels_start_offsets, _mm_set1_epi32(next_nb_info_ssb.inPos[0])),
            1);
          sum_abs = _mm_add_epi32(t, sum_abs);
          } break;
        default:
          assert(0);
      }
      sum_abs = _mm_and_si128(sum_abs, first_byte);
      if (extRiceFlag) {
        assert(0 && "Not implemented for avx2");
      } else {
        __m128i sum_all = _mm_max_epi32(
          _mm_min_epi32(
            _mm_set1_epi32(31),
            _mm_sub_epi32(sum_abs, _mm_set1_epi32(20))),
          _mm_set1_epi32(0));
        __m128i temp = _mm_i32gather_epi32(g_goRiceParsCoeff, sum_all, 4);
        __m128i control = _mm_setr_epi8(0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
        __m128i go_rice_par = _mm_shuffle_epi8(temp, control);
        int     go_rice_par_i = _mm_extract_epi32(go_rice_par, 0);
        memcpy(&state->m_goRicePar[state_offset], &go_rice_par_i, 4);
      }
    }

    else if (rem_reg_all_lt4) {
      uint8_t*       levels = (uint8_t*)state->m_absLevelsAndCtxInit[state_offset];
      const __m128i  last_two_bytes = _mm_set1_epi32(0xffff);
      const __m128i  last_byte = _mm_set1_epi32(0xff);
      const uint32_t tinit_offset = MIN(level_offset - 1u, 15u) + 8;
      const __m128i levels_start_offsets = _mm_set_epi32(48 * 3, 48 * 2, 48 * 1, 48 * 0);
      const __m128i ctx_start_offsets = _mm_srli_epi32(levels_start_offsets, 1);
      __m128i       tinit = _mm_i32gather_epi32(
        (int*)state->m_absLevelsAndCtxInit[state_offset],
        _mm_add_epi32(ctx_start_offsets, _mm_set1_epi32(tinit_offset)),
        2);
      tinit = _mm_and_si128(tinit, last_two_bytes);
      __m128i sum_abs = _mm_srli_epi32(tinit, 8);
      switch (numIPos) {
        case 5: {
          __m128i t = _mm_i32gather_epi32(
            (int*)levels,
            _mm_add_epi32(levels_start_offsets, _mm_set1_epi32(next_nb_info_ssb.inPos[4])),
            1);
          t = _mm_and_si128(t, last_byte);
          sum_abs = _mm_add_epi32(sum_abs, t);
        }
        case 4: {
          __m128i t = _mm_i32gather_epi32(
            (int*)levels,
            _mm_add_epi32(levels_start_offsets, _mm_set1_epi32(next_nb_info_ssb.inPos[3])),
            1);
          t = _mm_and_si128(t, last_byte);
          sum_abs = _mm_add_epi32(sum_abs, t);
        }
        case 3: {
          __m128i t = _mm_i32gather_epi32(
            (int*)levels,
            _mm_add_epi32(levels_start_offsets, _mm_set1_epi32(next_nb_info_ssb.inPos[2])),
            1);
          t = _mm_and_si128(t, last_byte);
          sum_abs = _mm_add_epi32(sum_abs, t);
        }
        case 2: {
          __m128i t = _mm_i32gather_epi32(
            (int*)levels,
            _mm_add_epi32(levels_start_offsets, _mm_set1_epi32(next_nb_info_ssb.inPos[1])),
            1);
          t = _mm_and_si128(t, last_byte);
          sum_abs = _mm_add_epi32(sum_abs, t);
        }
        case 1: {
          __m128i t = _mm_i32gather_epi32(
            (int*)levels,
            _mm_add_epi32(levels_start_offsets, _mm_set1_epi32(next_nb_info_ssb.inPos[0])),
            1);
          t = _mm_and_si128(t, last_byte);
          sum_abs = _mm_add_epi32(sum_abs, t);
        } break;
        default:
          assert(0);
      }
      if (extRiceFlag) {
        assert(0 && "Not implemented for avx2");
      } else {
        __m128i sum_all = _mm_min_epi32(_mm_set1_epi32(31), sum_abs);
        __m128i temp = _mm_i32gather_epi32(g_goRiceParsCoeff, sum_all, 4);
        __m128i control = _mm_setr_epi8(0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
        __m128i go_rice_par = _mm_shuffle_epi8(temp, control);
        int     go_rice_par_i = _mm_extract_epi32(go_rice_par, 0);
        memcpy(&state->m_goRicePar[state_offset], &go_rice_par_i, 4);

        
        for (int i = 0; i < 4; ++i) {
          state->m_goRiceZero[state_offset + i] = (i < 2 ? 1 : 2) << state->m_goRicePar[state_offset + i];
          
        }

      }

    }
    else {
      for (int i = 0; i < 4; ++i) {
        const int state_id = state_offset + i;
        uint8_t*  levels = (uint8_t*)(state->m_absLevelsAndCtxInit[state_id]);
        if (state->m_remRegBins[state_id] >= 4) {
          coeff_t tinit = state->m_absLevelsAndCtxInit[state_id][8 + ((scan_pos - 1) & 15)];
          coeff_t sumAbs1 = (tinit >> 3) & 31;
          coeff_t sumNum = tinit & 7;
#define UPDATE(k)                                  \
  {                                                \
    coeff_t t = levels[next_nb_info_ssb.inPos[k]]; \
    sumAbs1 += MIN(4 + (t & 1), t);                \
    sumNum += !!t;                                 \
  }
          switch (numIPos) {
            case 5: UPDATE(4);
            case 4: UPDATE(3);
            case 3: UPDATE(2);
            case 2: UPDATE(1);
            case 1: UPDATE(0); break;
            default: assert(0);
          }
#undef UPDATE
          coeff_t sumGt1 = sumAbs1 - sumNum;
          state->m_sigFracBits[state_id][0] = state->m_sigFracBitsArray[state_id][sigCtxOffsetNext + MIN((sumAbs1 + 1) >> 1, 3)][0];
          state->m_sigFracBits[state_id][1] = state->m_sigFracBitsArray[state_id][sigCtxOffsetNext + MIN((sumAbs1 + 1) >> 1, 3)][1];
          memcpy(state->m_coeffFracBits[state_id], state->m_gtxFracBitsArray[gtxCtxOffsetNext + (sumGt1 < 4 ? sumGt1 : 4)], sizeof(state->m_coeffFracBits[0]));


          coeff_t sumAbs = state->m_absLevelsAndCtxInit[state_id][8 + ((scan_pos - 1) & 15)] >> 8;
#define UPDATE(k)                                  \
  {                                                \
    coeff_t t = levels[next_nb_info_ssb.inPos[k]]; \
    sumAbs += t;                                   \
  }
          switch (numIPos) {
            case 5: UPDATE(4);
            case 4: UPDATE(3);
            case 3: UPDATE(2);
            case 2: UPDATE(1);
            case 1: UPDATE(0); break;
            default: assert(0);
          }
#undef UPDATE
          if (extRiceFlag) {
            unsigned currentShift = templateAbsCompare(sumAbs);
            sumAbs = sumAbs >> currentShift;
            int sumAll = MAX(MIN(31, (int)sumAbs - (int)baseLevel), 0);
            state->m_goRicePar[state_id] = g_goRiceParsCoeff[sumAll];
            state->m_goRicePar[state_id] += currentShift;
          } else {
            int sumAll = MAX(MIN(31, (int)sumAbs - 4 * 5), 0);
            state->m_goRicePar[state_id] = g_goRiceParsCoeff[sumAll];
          }
        } else {
          coeff_t sumAbs = (state->m_absLevelsAndCtxInit[state_id][8 + ((scan_pos - 1) & 15)]) >> 8;
#define UPDATE(k)                                  \
  {                                                \
    coeff_t t = levels[next_nb_info_ssb.inPos[k]]; \
    sumAbs += t;                                   \
  }
          switch (numIPos) {
            case 5: UPDATE(4);
            case 4: UPDATE(3);
            case 3: UPDATE(2);
            case 2: UPDATE(1);
            case 1: UPDATE(0); break;
            default: assert(0);
          }
#undef UPDATE
          if (extRiceFlag) {
            unsigned currentShift = templateAbsCompare(sumAbs);
            sumAbs = sumAbs >> currentShift;
            sumAbs = MIN(31, sumAbs);
            state->m_goRicePar[state_id] = g_goRiceParsCoeff[sumAbs];
            state->m_goRicePar[state_id] += currentShift;
          } else {
            sumAbs = MIN(31, sumAbs);
            state->m_goRicePar[state_id] = g_goRiceParsCoeff[sumAbs];
          }
          state->m_goRiceZero[state_id] = ((state_id & 3) < 2 ? 1 : 2) << state->m_goRicePar[state_id];
        }
      }
    }
  } else {
    for (int i = 0; i < 4; ++i) {
      state->all_gte_four = true;
      state->all_lt_four = true;
      updateState(
        ctxs,
        numIPos,
        scan_pos,
        decisions,
        sigCtxOffsetNext,
        gtxCtxOffsetNext,
        next_nb_info_ssb,
        baseLevel,
        extRiceFlag,
        i);
    }
  }
}


static INLINE void updateState(
  context_store * ctxs,
  int             numIPos,
  const uint32_t  scan_pos,
  const Decision* decisions,
  const uint32_t  sigCtxOffsetNext,
  const uint32_t  gtxCtxOffsetNext,
  const NbInfoSbb next_nb_info_ssb,
  const int       baseLevel,
  const bool      extRiceFlag,
  int             decision_id) {
  all_depquant_states* state = &ctxs->m_allStates;
  int state_id = ctxs->m_curr_state_offset + decision_id;
  // state->m_rdCost[state_id] = decisions->rdCost[decision_id];
  if (decisions->prevId[decision_id] > -2) {
    if (decisions->prevId[decision_id] >= 0) {
      const int prvState = ctxs->m_prev_state_offset + decisions->prevId[decision_id];
      state->m_numSigSbb[state_id] = (state->m_numSigSbb[prvState]) || !!decisions->absLevel[decision_id];
      state->m_refSbbCtxId[state_id] = state->m_refSbbCtxId[prvState];
      state->m_sbbFracBits[state_id][0] = state->m_sbbFracBits[prvState][0];
      state->m_sbbFracBits[state_id][1] = state->m_sbbFracBits[prvState][1];
      state->m_remRegBins[state_id] = state->m_remRegBins[prvState] - 1;
      state->m_goRicePar[state_id] = state->m_goRicePar[prvState];
      if (state->m_remRegBins[state_id] >= 4) {
        state->m_remRegBins[state_id] -= (decisions->absLevel[decision_id] < 2
                                            ? (unsigned)decisions->absLevel[decision_id]
                                            : 3);
      }
      memcpy(state->m_absLevelsAndCtxInit[state_id], state->m_absLevelsAndCtxInit[prvState], 48 * sizeof(uint8_t));
    }
    else {
      state->m_numSigSbb[state_id] = 1;
      state->m_refSbbCtxId[state_id] = -1;
      int ctxBinSampleRatio = 28;
      //(scanInfo.chType == CHANNEL_TYPE_LUMA) ? MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_LUMA : MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_CHROMA;
      state->m_remRegBins[state_id] = (state->effWidth * state->effHeight * ctxBinSampleRatio) / 16 - (
        decisions->absLevel[decision_id] < 2 ? (unsigned)decisions->absLevel[decision_id] : 3);
      memset(state->m_absLevelsAndCtxInit[state_id], 0, 48 * sizeof(uint8_t));
    }
    state->all_gte_four &= state->m_remRegBins[state_id] >= 4;
    state->all_lt_four &= state->m_remRegBins[state_id] < 4;
    uint8_t* levels = (uint8_t*)(state->m_absLevelsAndCtxInit[state_id]);
    levels[scan_pos & 15] = (uint8_t)MIN(32, decisions->absLevel[decision_id]);

    if (state->m_remRegBins[state_id] >= 4) {
      coeff_t tinit = state->m_absLevelsAndCtxInit[state_id][8 + ((scan_pos - 1) & 15)];
      coeff_t sumAbs1 = (tinit >> 3) & 31;
      coeff_t sumNum = tinit & 7;
#define UPDATE(k) {coeff_t t=levels[next_nb_info_ssb.inPos[k]]; sumAbs1+=MIN(4+(t&1),t); sumNum+=!!t; }
      switch (numIPos) {
        case 5: UPDATE(4);
        case 4: UPDATE(3);
        case 3: UPDATE(2);
        case 2: UPDATE(1);
        case 1: UPDATE(0); break;
        default: assert(0);
      }
#undef UPDATE
      coeff_t sumGt1 = sumAbs1 - sumNum;
      state->m_sigFracBits[state_id][0] = state->m_sigFracBitsArray[state_id][sigCtxOffsetNext + MIN(
        (sumAbs1 + 1) >> 1, 3)][0];
      state->m_sigFracBits[state_id][1] = state->m_sigFracBitsArray[state_id][sigCtxOffsetNext + MIN(
        (sumAbs1 + 1) >> 1, 3)][1];
      memcpy(state->m_coeffFracBits[state_id], state->m_gtxFracBitsArray[gtxCtxOffsetNext + (sumGt1 < 4 ? sumGt1 : 4)],
             sizeof(state->m_coeffFracBits[0]));


      coeff_t sumAbs = state->m_absLevelsAndCtxInit[state_id][8 + ((scan_pos - 1) & 15)] >> 8;
#define UPDATE(k) {coeff_t t=levels[next_nb_info_ssb.inPos[k]]; sumAbs+=t; }
      switch (numIPos) {
        case 5: UPDATE(4);
        case 4: UPDATE(3);
        case 3: UPDATE(2);
        case 2: UPDATE(1);
        case 1: UPDATE(0); break;
        default: assert(0);
      }
#undef UPDATE
      if (extRiceFlag) {
        unsigned currentShift = templateAbsCompare(sumAbs);
        sumAbs = sumAbs >> currentShift;
        int sumAll = MAX(MIN(31, (int)sumAbs - (int)baseLevel), 0);
        state->m_goRicePar[state_id] = g_goRiceParsCoeff[sumAll];
        state->m_goRicePar[state_id] += currentShift;
      }
      else {
        int sumAll = MAX(MIN(31, (int)sumAbs - 4 * 5), 0);
        state->m_goRicePar[state_id] = g_goRiceParsCoeff[sumAll];
      }
    }
    else {
      coeff_t sumAbs = (state->m_absLevelsAndCtxInit[state_id][8 + ((scan_pos - 1) & 15)]) >> 8;
#define UPDATE(k) {coeff_t t=levels[next_nb_info_ssb.inPos[k]]; sumAbs+=t; }
      switch (numIPos) {
        case 5: UPDATE(4);
        case 4: UPDATE(3);
        case 3: UPDATE(2);
        case 2: UPDATE(1);
        case 1: UPDATE(0); break;
        default: assert(0);
      }
#undef UPDATE
      if (extRiceFlag) {
        unsigned currentShift = templateAbsCompare(sumAbs);
        sumAbs = sumAbs >> currentShift;
        sumAbs = MIN(31, sumAbs);
        state->m_goRicePar[state_id] = g_goRiceParsCoeff[sumAbs];
        state->m_goRicePar[state_id] += currentShift;
      }
      else {
        sumAbs = MIN(31, sumAbs);
        state->m_goRicePar[state_id] = g_goRiceParsCoeff[sumAbs];
      }
      state->m_goRiceZero[state_id] = ((state_id & 3) < 2 ? 1 : 2) << state->m_goRicePar[state_id];
    }
  }
  else {
    state->all_gte_four &= state->m_remRegBins[state_id] >= 4;
    state->all_lt_four &= state->m_remRegBins[state_id] < 4;
  }
}

static bool same[13];
static void xDecideAndUpdate(
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

  if (scan_pos) {
    if (!(scan_pos & 15)) {
      SWAP(ctxs->m_common_context.m_curr_sbb_ctx_offset, ctxs->m_common_context.m_prev_sbb_ctx_offset, int);
      update_state_eos_avx2(ctxs, scan_pos, scan_info->cg_pos, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], width_in_sbb, height_in_sbb, scan_info->next_sbb_right, scan_info->next_sbb_below, decisions);
      //updateStateEOS(ctxs, scan_pos, scan_info->cg_pos, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], width_in_sbb, height_in_sbb, scan_info->next_sbb_right, scan_info->next_sbb_below, decisions, 0);
      //updateStateEOS(ctxs, scan_pos, scan_info->cg_pos, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], width_in_sbb, height_in_sbb, scan_info->next_sbb_right, scan_info->next_sbb_below, decisions, 1);
      //updateStateEOS(ctxs, scan_pos, scan_info->cg_pos, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], width_in_sbb, height_in_sbb, scan_info->next_sbb_right, scan_info->next_sbb_below, decisions, 2);
      //updateStateEOS(ctxs, scan_pos, scan_info->cg_pos, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], width_in_sbb, height_in_sbb, scan_info->next_sbb_right, scan_info->next_sbb_below, decisions, 3);
      memcpy(decisions->prevId + 4, decisions->prevId, 4 * sizeof(int32_t));
      memcpy(decisions->absLevel + 4, decisions->absLevel, 4 * sizeof(int32_t));
      memcpy(decisions->rdCost + 4, decisions->rdCost, 4 * sizeof(int64_t));
    } else if (!zeroOut) {
      update_states_avx2(ctxs, next_nb_info_ssb.num, scan_pos, decisions, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], next_nb_info_ssb, 4, false);
    /*  updateState(ctxs, next_nb_info_ssb.num, scan_pos, decisions, sigCtxOffsetNext, gtxCtxOffsetNext, next_nb_info_ssb, 4, false, 0);
      updateState(ctxs, next_nb_info_ssb.num, scan_pos, decisions, sigCtxOffsetNext, gtxCtxOffsetNext, next_nb_info_ssb, 4, false, 1);
      updateState(ctxs, next_nb_info_ssb.num, scan_pos, decisions, sigCtxOffsetNext, gtxCtxOffsetNext, next_nb_info_ssb, 4, false, 2);
      updateState(ctxs, next_nb_info_ssb.num, scan_pos, decisions, sigCtxOffsetNext, gtxCtxOffsetNext, next_nb_info_ssb, 4, false, 3);*/
    }

    if (spt == SCAN_SOCSBB) {
      SWAP(ctxs->m_skip_state_offset, ctxs->m_prev_state_offset, int);
    }
  }
}


int uvg_dep_quant(
  const encoder_state_t* const state,
  const cu_info_t* const cur_tu,
  const int width,
  const int height,
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
  dep_quant_context.m_curr_state_offset = 0;
  dep_quant_context.m_prev_state_offset = 4;
  dep_quant_context.m_skip_state_offset = 8;
   
  const uint32_t  lfnstIdx = tree_type != UVG_CHROMA_T  || compID == COLOR_Y ?
                               cur_tu->lfnst_idx :
                               cur_tu->cr_lfnst_idx;
  
  const int       numCoeff = width * height;

  memset(coeff_out, 0x00, width * height * sizeof(coeff_t));
  *absSum                    = 0;

  const bool      is_mts   = compID == COLOR_Y && cur_tu->tr_idx > MTS_SKIP;
  const bool      is_ts    = (cur_tu->tr_skip >> compID) & 1;

  const uint32_t  log2_tr_width  = uvg_g_convert_to_log2[width];
  const uint32_t  log2_tr_height = uvg_g_convert_to_log2[height];
  const uint32_t* const scan     = uvg_get_scan_order_table(SCAN_GROUP_4X4,0,log2_tr_width,log2_tr_height);
  const uint32_t* const cg_scan     = uvg_get_scan_order_table(SCAN_GROUP_UNGROUPED,0,log2_tr_width,log2_tr_height);

  int32_t qp_scaled = uvg_get_scaled_qp(compID, state->qp, (encoder->bitdepth - 8) * 6, encoder->qp_map[0]);
  qp_scaled = is_ts ? MAX(qp_scaled, 4 + 6 * MIN_QP_PRIME_TS) : qp_scaled;
  bool needs_block_size_trafo_scale = !is_ts && ((log2_tr_height + log2_tr_width) % 2 == 1);
  needs_block_size_trafo_scale |= 0; // Non log2 block size

  const int32_t scalinglist_type = (cur_tu->type == CU_INTRA ? 0 : 3) + (int8_t)compID;
  const int32_t *q_coeff = encoder->scaling_list.quant_coeff[log2_tr_width][log2_tr_height][scalinglist_type][qp_scaled % 6];

  if (compID != COLOR_Y) {
    dep_quant_context.m_quant = (quant_block*)& state->quant_blocks[2];
  } else if (cur_tu->type == CU_INTRA && cur_tu->intra.isp_mode != ISP_MODE_NO_ISP) {
    dep_quant_context.m_quant = (quant_block*)&state->quant_blocks[1];    
  } else {
    dep_quant_context.m_quant = (quant_block*)&state->quant_blocks[0];   
  }
  //TODO: no idea when it is safe not to reinit for inter
  if (dep_quant_context.m_quant->needs_init || cur_tu->type == CU_INTER) {
    init_quant_block(state, dep_quant_context.m_quant, cur_tu, log2_tr_width, log2_tr_height, compID, needs_block_size_trafo_scale, -1);
  }
  
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
     (state->encoder_control->cfg.mts && 0 /*sbt used by block*/ &&
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
  const int32_t default_quant_coeff = dep_quant_context.m_quant->m_QScale;
  const int32_t thres               = dep_quant_context.m_quant->m_thresLast;
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
  rate_estimator_t* rate_estimator = (rate_estimator_t *)(compID == COLOR_Y && cur_tu->type == CU_INTRA && cur_tu->intra.isp_mode != ISP_MODE_NO_ISP ?
    &state->rate_estimator[3] : &state->rate_estimator[compID]);
  if(rate_estimator->needs_init || cur_tu->type == CU_INTER) {
    init_rate_esimator(rate_estimator, &state->search_cabac, compID);
    xSetLastCoeffOffset(state, cur_tu, width, height, rate_estimator, compID);
    rate_estimator->needs_init = false;
  }

  reset_common_context(&dep_quant_context.m_common_context, rate_estimator, (width * height) >> 4, numCoeff);
  dep_quant_context.m_common_context.m_nbInfo = encoder->m_scanId2NbInfoOutArray[log2_tr_width][log2_tr_height];
  

  int effectHeight = MIN(32, effHeight);
  int effectWidth = MIN(32, effWidth);
  for (int k = 0; k < 12; k++) {
    dep_quant_context.m_allStates.m_rdCost[k] = INT64_MAX >> 1;
    dep_quant_context.m_allStates.m_numSigSbb[k] = 0;
    dep_quant_context.m_allStates.m_remRegBins[k] = 4; // just large enough for last scan pos
    dep_quant_context.m_allStates.m_refSbbCtxId[k] = -1;
    dep_quant_context.m_allStates.m_sigFracBits[k][0] = rate_estimator->m_sigFracBits[0][0][0];
    dep_quant_context.m_allStates.m_sigFracBits[k][1] = rate_estimator->m_sigFracBits[0][0][1];
    memcpy(dep_quant_context.m_allStates.m_coeffFracBits[k], rate_estimator->m_gtxFracBits[0], sizeof(dep_quant_context.m_allStates.m_coeffFracBits[k]));
    dep_quant_context.m_allStates.m_goRicePar[k] = 0;
    dep_quant_context.m_allStates.m_goRiceZero[k] = 0;

    dep_quant_context.m_allStates.m_sbbFracBits[k][0] = 0;
    dep_quant_context.m_allStates.m_sbbFracBits[k][1] = 0;

    dep_quant_context.m_allStates.m_stateId[k] = k & 3;
    for (int i = 0; i < (compID == COLOR_Y ? 12 : 8); ++i) {
      memcpy(dep_quant_context.m_allStates.m_sigFracBitsArray[k][i], rate_estimator->m_sigFracBits[(k & 3 ? (k & 3) - 1 : 0)][i], sizeof(uint32_t) * 2);
    }
  }

  dep_quant_context.m_allStates.effHeight = effectHeight;
  dep_quant_context.m_allStates.effWidth = effectWidth;
  dep_quant_context.m_allStates.all_gte_four = true;
  dep_quant_context.m_allStates.all_lt_four = false;
  dep_quant_context.m_allStates.m_commonCtx = &dep_quant_context.m_common_context;
  for (int i = 0; i < (compID == COLOR_Y ? 21 : 11); ++i) {
    memcpy(dep_quant_context.m_allStates.m_gtxFracBitsArray[i], rate_estimator->m_gtxFracBits[i], sizeof(int32_t) * 6);
  }

  depquant_state_init(&dep_quant_context.m_startState, rate_estimator->m_sigFracBits[0][0], rate_estimator->m_gtxFracBits[0]);
  dep_quant_context.m_startState.effHeight = effectHeight;
  dep_quant_context.m_startState.effWidth = effectWidth;
  dep_quant_context.m_startState.m_stateId = 0;
  dep_quant_context.m_startState.m_commonCtx = &dep_quant_context.m_common_context;
  for (int i = 0; i < (compID == COLOR_Y ? 12 : 8); ++i) {
    dep_quant_context.m_startState.m_sigFracBitsArray[i] = rate_estimator->m_sigFracBits[0][i];
  }
  for (int i = 0; i < (compID == COLOR_Y ? 21 : 11); ++i) {
    dep_quant_context.m_startState.m_gtxFracBitsArray[i] = rate_estimator->m_gtxFracBits[i];
  }

  const uint32_t height_in_sbb = MAX(height >> 2, 1);
  const uint32_t width_in_sbb = MAX(width >> 2, 1);

  //===== populate trellis =====
  for (int scanIdx = firstTestPos; scanIdx >= 0; scanIdx--) {
    uint32_t blkpos = scan[scanIdx];
    struct dep_quant_scan_info* scan_info = &encoder->scan_info[log2_tr_width][log2_tr_height][scanIdx];

    context_store* ctxs = &dep_quant_context;
    if (enableScalingLists) {
      init_quant_block(state, dep_quant_context.m_quant, cur_tu, log2_tr_width, log2_tr_height, compID, needs_block_size_trafo_scale, q_coeff[blkpos]);

      xDecideAndUpdate(
        rate_estimator,
        ctxs,
        scan_info,
        abs(srcCoeff[blkpos]),
        scanIdx,
        width_in_sbb,
        height_in_sbb,
        encoder->m_scanId2NbInfoSbbArray[log2_tr_width][log2_tr_height][scanIdx ? scanIdx - 1 : 0],
        (zeroOut && (scan_info->pos_x >= effWidth || scan_info->pos_y >= effHeight)),
        q_coeff[blkpos],
        width,
        height,
        compID != 0
        ); //tu.cu->slice->getReverseLastSigCoeffFlag());
    }
    else {
      xDecideAndUpdate(
        rate_estimator,
        ctxs,
        scan_info,
        abs(srcCoeff[blkpos]),
        scanIdx,
        width_in_sbb,
        height_in_sbb,
        encoder->m_scanId2NbInfoSbbArray[log2_tr_width][log2_tr_height][scanIdx ? scanIdx - 1 : 0],
        (zeroOut && (scan_info->pos_x >= effWidth || scan_info->pos_y >= effHeight)),
        default_quant_coeff,
        width,
        height,
        compID != 0); //tu.cu->slice->getReverseLastSigCoeffFlag());
    }
    if(0){
      printf("%d\n", scanIdx);
      for (int i = 0; i < 4; i++) {
        printf("%lld %hu %d\n", ctxs->m_trellis[scanIdx].rdCost[i], ctxs->m_trellis[scanIdx].absLevel[i], ctxs->m_trellis[scanIdx].prevId[i]);
      }
      printf("\n");
    }
  }

  //===== find best path =====
  int prev_id    = -1;
  int64_t  minPathCost = 0;  
  for (int8_t stateId = 0; stateId < 4; stateId++) {
    int64_t pathCost = dep_quant_context.m_trellis[0].rdCost[stateId];
    if (pathCost < minPathCost) {
      prev_id = stateId;
      minPathCost     = pathCost;
    }
  }

  //===== backward scanning =====
  int scanIdx = 0;
  context_store* ctxs = &dep_quant_context;
  //  printf("%d\n", scanIdx);
  //for (int i = 0; i < 4; i++) {
  //  printf("%lld %hu %d\n", ctxs->m_trellis[scanIdx].rdCost[i], ctxs->m_trellis[scanIdx].absLevel[i], ctxs->m_trellis[scanIdx].prevId[i]);
  //}
  //printf("\n");
  for (; prev_id >= 0; scanIdx++) {
    Decision temp       = dep_quant_context.m_trellis[scanIdx];
    int32_t blkpos = scan[scanIdx];
    coeff_out[blkpos] = (srcCoeff[blkpos] < 0 ? -temp.absLevel[prev_id] : temp.absLevel[prev_id]);
    *absSum += temp.absLevel[prev_id];
    prev_id = temp.prevId[prev_id];
  }
  return *absSum;
}


void uvg_dep_quant_dequant(
  const encoder_state_t* const state,
  const int block_type,
  const int width,
  const int height,
  const color_t compID,
  coeff_t* quant_coeff,
  coeff_t * coeff, 
  bool enableScalingLists)
{
  const encoder_control_t* const encoder = state->encoder_control;

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
  const int         qpDQ = uvg_get_scaled_qp(compID, state->qp, (encoder->bitdepth - 8) * 6, encoder->qp_map[0]) + 1;
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
  int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)(compID);

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
