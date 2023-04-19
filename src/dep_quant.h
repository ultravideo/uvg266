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
#endif
