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

#ifndef DEP_QUANT_H_
#define DEP_QUANT_H_

#include "cu.h"
#include "global.h"

#define SM_NUM_CTX_SETS_SIG   3
#define SM_NUM_CTX_SETS_GTX   2
#define SM_MAX_NUM_SIG_SBB_CTX 2
#define SM_MAX_NUM_SIG_CTX    12
#define SM_MAX_NUM_GTX_CTX    21
#define SCALE_BITS         15
#define RICEMAX            32

typedef struct encoder_control_t encoder_control_t;

enum ScanPosType { SCAN_ISCSBB = 0, SCAN_SOCSBB = 1, SCAN_EOCSBB = 2 };

struct dep_quant_scan_info
{
  uint8_t sig_ctx_offset[2];
  uint8_t gtx_ctx_offset[2];
  uint16_t cg_pos;
  uint16_t  pos_y;
  uint16_t  pos_x;
  uint8_t next_sbb_right;
  uint8_t next_sbb_below;
};

typedef struct
{
  int     m_QShift;
  int64_t m_QAdd;
  int64_t m_QScale;
  int64_t m_maxQIdx;
  int64_t m_thresLast;
  int64_t m_thresSSbb;
  // distortion normalization
  int     m_DistShift;
  int64_t m_DistAdd;
  int64_t m_DistStepAdd;
  int64_t m_DistOrgFact;
  bool    needs_init;
} quant_block;

typedef struct
{
  int32_t  m_lastBitsX[TR_MAX_WIDTH];
  int32_t  m_lastBitsY[TR_MAX_WIDTH];
  uint32_t m_sigSbbFracBits[SM_MAX_NUM_SIG_SBB_CTX][2];
  uint32_t m_sigFracBits[SM_NUM_CTX_SETS_SIG][SM_MAX_NUM_SIG_CTX][2];
  int32_t  m_gtxFracBits[SM_MAX_NUM_GTX_CTX][6];
  bool     needs_init;
} rate_estimator_t;


typedef struct
{
  uint8_t num;
  uint8_t inPos[5];
} NbInfoSbb;

typedef struct
{
  uint16_t maxDist;
  uint16_t num;
  uint16_t outPos[5];
} NbInfoOut;

typedef struct {
  int32_t absLevel[4];
  int64_t deltaDist[4];
} PQData;

typedef struct {
  int64_t ALIGNED(32) rdCost[8];
  int32_t ALIGNED(32) absLevel[8];
  int32_t ALIGNED(32) prevId[8];
  uint8_t zero_out;
} Decision;


typedef struct {
  uint8_t* sbbFlags;
  uint8_t* levels;
} SbbCtx;

typedef struct {
  const NbInfoOut* m_nbInfo;
  uint32_t         m_sbbFlagBits[2][2];
  SbbCtx           m_allSbbCtx[2];
  int              m_curr_sbb_ctx_offset;
  int              m_prev_sbb_ctx_offset;
  uint8_t          sbb_memory[8 * 1024];
  uint8_t          level_memory[8 * TR_MAX_WIDTH * TR_MAX_WIDTH];
} common_context;


typedef struct {
  int64_t  m_rdCost;
  uint16_t m_absLevelsAndCtxInit[24]; // 16x8bit for abs levels + 16x16bit for ctx init id
  int8_t          m_numSigSbb;
  int             m_remRegBins;
  int8_t          m_refSbbCtxId;
  uint32_t        m_sbbFracBits[2];
  uint32_t        m_sigFracBits[2];
  int32_t         m_coeffFracBits[6];
  int8_t          m_goRicePar;
  int8_t          m_goRiceZero;
  int8_t          m_stateId;
  uint32_t*       m_sigFracBitsArray[12];
  int32_t*        m_gtxFracBitsArray[21];
  common_context* m_commonCtx;

  unsigned        effWidth;
  unsigned        effHeight;
} depquant_state;
typedef struct {
  int64_t  ALIGNED(32) m_rdCost[12];
  uint8_t  ALIGNED(32) m_absLevels[3][16 * 4]; 
  uint16_t ALIGNED(32) m_ctxInit[3][16 * 4]; 
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

  unsigned        effWidth;
  unsigned        effHeight;

  bool            all_gte_four;
  bool            all_lt_four;
} all_depquant_states;

typedef struct {
  common_context      m_common_context;
  all_depquant_states m_allStates;
  int                 m_curr_state_offset;
  int                 m_prev_state_offset;
  int                 m_skip_state_offset;
  depquant_state      m_startState;
  quant_block*        m_quant;
  Decision            m_trellis[TR_MAX_WIDTH * TR_MAX_WIDTH];
} context_store;


int uvg_init_nb_info(encoder_control_t* encoder);
void uvg_dealloc_nb_info(encoder_control_t* encoder);


void uvg_dep_quant_dequant(
  const encoder_state_t* const state,
  const int block_type,
  const int width,
  const int height,
  const color_t compID,
  coeff_t* quant_coeff,
  coeff_t* coeff,
  bool enableScalingLists);

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
  const bool enableScalingLists);


void uvg_dep_quant_update_state(
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


void uvg_dep_quant_update_state_eos(
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

void uvg_dep_quant_check_rd_costs(
  const all_depquant_states* const state,
  const enum ScanPosType           spt,
  const PQData*                    pqDataA,
  Decision*                        decisions,
  const int                        decisionA,
  const int                        decisionB,
  const int                        state_offset);
#endif
