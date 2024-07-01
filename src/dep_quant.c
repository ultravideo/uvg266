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

#include "strategies-depquant.h"
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
  for (int k = 0; k < 2; k++, next_sbb_memory += numSbb * 4llu, next_level_memory += num_coeff * 4llu) {
    ctx->m_allSbbCtx[k].sbbFlags = next_sbb_memory;
    ctx->m_allSbbCtx[k].levels = next_level_memory;
  }
  ctx->m_curr_sbb_ctx_offset = 0;
  ctx->m_prev_sbb_ctx_offset = 1;
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


void uvg_dep_quant_check_rd_costs(
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
  const int      curr_state_without_offset = curr_state & 3;
  uint8_t*       sbbFlags = cc->m_allSbbCtx[cc->m_curr_sbb_ctx_offset].sbbFlags;
  uint8_t*       levels = cc->m_allSbbCtx[cc->m_curr_sbb_ctx_offset].levels;
  size_t         setCpSize = cc->m_nbInfo[scan_pos - 1].maxDist * sizeof(uint8_t);
  if (prev_state != -1) {
    const int8_t prev_sbb_state = ctxs->m_allStates.m_refSbbCtxId[prev_state];
    for (int i = 0; i < numSbb; ++i) {
      sbbFlags[i * 4 + curr_state_without_offset] = cc->m_allSbbCtx[cc->m_prev_sbb_ctx_offset].sbbFlags[i * 4 + prev_sbb_state];
    }
    for (int i = 16; i < setCpSize; ++i) {
      levels[scan_pos * 4 + i * 4 + curr_state_without_offset] = cc->m_allSbbCtx[cc->m_prev_sbb_ctx_offset].sbbFlags[scan_pos * 4 + i * 4 + prev_sbb_state];
    }
  }
  else {
    for (int i = 0; i < numSbb; ++i) {
      sbbFlags[i * 4 + curr_state_without_offset] = 0;
    }
    for (int i = 16; i < setCpSize; ++i) {
      levels[scan_pos * 4 + i * 4 + curr_state_without_offset] = 0;
    }
  }
  sbbFlags[cg_pos * 4 + curr_state_without_offset] = !!ctxs->m_allStates.m_numSigSbb[curr_state];
  for (int i = 0; i < 16; ++i) {
    levels[scan_pos * 4 + i * 4 + curr_state_without_offset] = ctxs->m_allStates.m_absLevels[curr_state / 4][i * 4 + curr_state_without_offset];
  }

  const int sigNSbb = ((next_sbb_right ? sbbFlags[next_sbb_right * 4 + curr_state_without_offset] : false) 
                       || (next_sbb_below ? sbbFlags[next_sbb_below* 4 + curr_state_without_offset] : false) ? 1 : 0);
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
  ctxs->m_allStates.m_refSbbCtxId[curr_state] = curr_state_without_offset;
  ctxs->m_allStates.m_sbbFracBits[curr_state][0] = cc->m_sbbFlagBits[sigNSbb][0];
  ctxs->m_allStates.m_sbbFracBits[curr_state][1] = cc->m_sbbFlagBits[sigNSbb][1];

  uint16_t *templateCtxInit = ctxs->m_allStates.m_ctxInit[ctxs->m_curr_state_offset >> 2];
  const int scanBeg = scan_pos - 16;
  const NbInfoOut* nbOut = cc->m_nbInfo + scanBeg;
  const uint8_t* absLevels = levels + scanBeg * 4;
  for (int id = 0; id < 16; id++, nbOut++) {
    if (nbOut->num) {
      coeff_t sumAbs = 0, sumAbs1 = 0, sumNum = 0;
#define UPDATE(k) {coeff_t t=absLevels[nbOut->outPos[k] * 4 + curr_state_without_offset]; sumAbs+=t; sumAbs1+=MIN(4+(t&1),t); sumNum+=!!t; }
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
      templateCtxInit[curr_state_without_offset + id * 4] = (uint16_t)(sumNum) + ((uint16_t)(sumAbs1 << 3)) + (uint16_t)(MIN(127, sumAbs) << 8);
    }
    else {
      templateCtxInit[curr_state_without_offset + id * 4] = 0;
    }
  }
  for (int i = curr_state_without_offset; i < 64; i += 4) {
    ctxs->m_allStates.m_absLevels[curr_state >> 2][i] = 0;
  }
}



void uvg_dep_quant_update_state_eos(
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
      for (int i = decision_id; i < 64;  i += 4) {
        state->m_absLevels[ctxs->m_curr_state_offset / 4][i] = 0;
      }
    }
    else if (decisions->prevId[decision_id] >= 0) {
      prvState = ctxs->m_prev_state_offset + decisions->prevId[decision_id];
      state->m_numSigSbb[curr_state_offset] = state->m_numSigSbb[prvState] || !!decisions->absLevel[decision_id];
      for (int i = 0; i < 64;  i += 4) {
        state->m_absLevels[ctxs->m_curr_state_offset / 4][i + decision_id] =
          state->m_absLevels[ctxs->m_prev_state_offset / 4][i + decisions->prevId[decision_id]];
      }
    }
    else {
      state->m_numSigSbb[curr_state_offset] = 1;
      for (int i = decision_id; i < 64; i += 4) {
        state->m_absLevels[ctxs->m_curr_state_offset / 4][i] = 0;
      }
    }
    uint8_t* temp = &state->m_absLevels[ctxs->m_curr_state_offset / 4][(scan_pos & 15) * 4 + decision_id];
    *temp = (uint8_t)MIN(255, decisions->absLevel[decision_id]);

    update_common_context(ctxs, state->m_commonCtx, scan_pos, cg_pos, width_in_sbb, height_in_sbb, next_sbb_right,
                          next_sbb_below, prvState, ctxs->m_curr_state_offset + decision_id);

    coeff_t tinit = state->m_ctxInit[ctxs->m_curr_state_offset  >> 2][((scan_pos - 1) & 15) * 4 + decision_id];
    coeff_t sumNum = tinit & 7;
    coeff_t sumAbs1 = (tinit >> 3) & 31;
    coeff_t sumGt1 = sumAbs1 - sumNum;
    state->m_sigFracBits[curr_state_offset][0] = state->m_sigFracBitsArray[curr_state_offset][sigCtxOffsetNext + MIN((sumAbs1 + 1) >> 1, 3)][0];
    state->m_sigFracBits[curr_state_offset][1] = state->m_sigFracBitsArray[curr_state_offset][sigCtxOffsetNext + MIN((sumAbs1 + 1) >> 1, 3)][1];

    memcpy(state->m_coeffFracBits[curr_state_offset],
           state->m_gtxFracBitsArray[gtxCtxOffsetNext + (sumGt1 < 4 ? sumGt1 : 4)], sizeof(state->m_coeffFracBits[0]));
  }
}


void uvg_dep_quant_update_state(
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
  all_depquant_states* state    = &ctxs->m_allStates;
  int                  state_id = ctxs->m_curr_state_offset + decision_id;
  state->m_rdCost[state_id]     = decisions->rdCost[decision_id];
  int32_t prev_id_no_offset     = decisions->prevId[decision_id];
  if (prev_id_no_offset > -2) {
    if (prev_id_no_offset >= 0) {
      const int prvState = ctxs->m_prev_state_offset + prev_id_no_offset;
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
      for (int i = 0; i < 64; i += 4) {
        state->m_ctxInit[ctxs->m_curr_state_offset  >> 2][decision_id + i] = state->m_ctxInit[ctxs->m_prev_state_offset  >> 2][prev_id_no_offset + i];
      }
      for (int i = 0; i < 64; i += 4) {
        state->m_absLevels[ctxs->m_curr_state_offset  >> 2][decision_id + i] = state->m_absLevels[ctxs->m_prev_state_offset  >> 2][prev_id_no_offset + i];
      }
    }
    else {
      state->m_numSigSbb[state_id] = 1;
      state->m_refSbbCtxId[state_id] = -1;
      int ctxBinSampleRatio = 28;
      //(scanInfo.chType == CHANNEL_TYPE_LUMA) ? MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_LUMA : MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_CHROMA;
      state->m_remRegBins[state_id] = (state->effWidth * state->effHeight * ctxBinSampleRatio) / 16 - (
        decisions->absLevel[decision_id] < 2 ? (unsigned)decisions->absLevel[decision_id] : 3);
      for (int i = decision_id; i < 64; i += 4) {
        state->m_absLevels[ctxs->m_curr_state_offset >> 2][i] = 0;
      }
      for (int i = decision_id; i < 64; i += 4) {
        state->m_ctxInit[ctxs->m_curr_state_offset >> 2][i] = 0;
      }
    }
    state->all_gte_four &= state->m_remRegBins[state_id] >= 4;
    state->all_lt_four &= state->m_remRegBins[state_id] < 4;
    uint8_t* levels = state->m_absLevels[ctxs->m_curr_state_offset >> 2];
    levels[(scan_pos & 15) * 4 + decision_id] = (uint8_t)MIN(32, decisions->absLevel[decision_id]);

    if (state->m_remRegBins[state_id] >= 4) {
      coeff_t tinit = state->m_ctxInit[ctxs->m_curr_state_offset  >> 2][((scan_pos - 1) & 15) * 4 + decision_id];
      coeff_t sumAbs1 = (tinit >> 3) & 31;
      coeff_t sumNum = tinit & 7;
#define UPDATE(k) {coeff_t t=levels[next_nb_info_ssb.inPos[k] * 4 + decision_id]; sumAbs1+=MIN(4+(t&1),t); sumNum+=!!t; }
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


      coeff_t sumAbs = state->m_ctxInit[ctxs->m_curr_state_offset  >> 2][((scan_pos - 1) & 15) * 4 + decision_id] >> 8;
#define UPDATE(k) {coeff_t t=levels[next_nb_info_ssb.inPos[k] * 4 + decision_id]; sumAbs+=t; }
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
      coeff_t sumAbs = state->m_ctxInit[ctxs->m_curr_state_offset  >> 2][((scan_pos - 1) & 15) * 4 + decision_id] >> 8;
#define UPDATE(k) {coeff_t t=levels[next_nb_info_ssb.inPos[k] * 4 + decision_id]; sumAbs+=t; }
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
  uvg_find_first_non_zero_coeff(
    srcCoeff,
    enableScalingLists,
    &dep_quant_context,
    scan,
    q_coeff,
    &firstTestPos,
    width, 
    height);
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
  } else if (compID == COLOR_U && state->encoder_control->cfg.jccr) {
    xSetLastCoeffOffset(state, cur_tu, width, height, rate_estimator, compID);    
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

  const int      default_quant_coeff = dep_quant_context.m_quant->m_QScale;
  //===== populate trellis =====
  for (int scanIdx = firstTestPos; scanIdx >= 0; scanIdx--) {
    uint32_t blkpos = scan[scanIdx];
    struct dep_quant_scan_info* scan_info = &encoder->scan_info[log2_tr_width][log2_tr_height][scanIdx];

    context_store* ctxs = &dep_quant_context;
    if (enableScalingLists) {
      init_quant_block(state, dep_quant_context.m_quant, cur_tu, log2_tr_width, log2_tr_height, compID, needs_block_size_trafo_scale, q_coeff[blkpos]);

      uvg_dep_quant_decide_and_update(
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
      uvg_dep_quant_decide_and_update(
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
      int  qIdx = (level * 2) + (level > 0 ? -(state >> 1) : (state >> 1));
      int64_t  nomTCoeff = ((int64_t)qIdx * (int64_t)invQScale + add) >> ((shift < 0) ? 0 : shift);
      coeff[rasterPos] = (coeff_t)CLIP(minTCoeff, maxTCoeff, nomTCoeff);
    }
    state = (32040 >> ((state << 2) + ((level & 1) << 1))) & 3;   // the 16-bit value "32040" represent the state transition table
  }
}
