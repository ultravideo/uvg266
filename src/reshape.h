#pragma once
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

/**
* \ingroup Reconstruction
* \file
* LMCS reshape.
*/

#include "checkpoint.h"
#include "cu.h"
#include "encoder.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "uvg266.h"
#include "videoframe.h"
#include "image.h"
#include "nal.h"

typedef struct lmcs_seq_info
{
  double binVar[32];
  double binHist[32];
  double normVar[32];
  int    nonZeroCnt;
  double weightVar;
  double weightNorm;
  double minBinVar;
  double maxBinVar;
  double meanBinVar;
  double ratioStdU;
  double ratioStdV;
} lmcs_seq_info;

typedef struct SliceReshapeInfo {
  bool      sliceReshaperEnableFlag;
  bool      sliceReshaperModelPresentFlag;
  unsigned  enableChromaAdj;
  uint32_t  reshaperModelMinBinIdx;
  uint32_t  reshaperModelMaxBinIdx;
  int       reshaperModelBinCWDelta[PIC_CODE_CW_BINS];
  int       maxNbitsNeededDeltaCW;
  int       chrResScalingOffset;
} SliceReshapeInfo;

typedef struct ReshapeCW
{
  uint32_t  binCW[3];
  int       updateCtrl;
  int       adpOption;
  uint32_t  initialCW;
  int rspPicSize;
  int rspFps;
  int rspBaseQP;
  int rspTid;
  int rspSliceQP;
  int rspFpsToIp;
} ReshapeCW;

typedef struct lmcs_aps {
  SliceReshapeInfo        m_sliceReshapeInfo;
  bool                    m_CTUFlag;
  bool                    m_recReshaped;
  uvg_pixel               m_invLUT[1024];
  uvg_pixel               m_fwdLUT[1024];
  int32_t                 m_chromaAdjHelpLUT[PIC_CODE_CW_BINS];
  uint16_t                m_binCW[PIC_ANALYZE_CW_BINS];
  uint16_t                m_initCW;
  bool                    m_reshape;
  uvg_pixel               m_reshapePivot[PIC_CODE_CW_BINS + 1];
  uvg_pixel               m_inputPivot[PIC_CODE_CW_BINS + 1];
  int32_t                 m_fwdScaleCoef[PIC_CODE_CW_BINS];
  int32_t                 m_invScaleCoef[PIC_CODE_CW_BINS];
  int                     m_lumaBD;
  int                     m_reshapeLUTSize;
  int                     m_chromaScale;
  int                     m_vpduX;
  int                     m_vpduY;

  bool                    m_exceedSTD;
  uint32_t                m_binImportance[PIC_ANALYZE_CW_BINS];
  int                     m_tcase;
  int                     m_rateAdpMode;
  bool                    m_useAdpCW;
  uint16_t                m_initCWAnalyze;
  ReshapeCW               m_reshapeCW;
  uvg_pixel               m_cwLumaWeight[PIC_CODE_CW_BINS];
  double                  m_chromaWeight;
  int                     m_chromaAdj;
  int                     m_binNum;
  lmcs_seq_info           m_srcSeqStats;
  lmcs_seq_info           m_rspSeqStats;

} lmcs_aps;


void uvg_free_lmcs_aps(lmcs_aps* aps);

void uvg_init_lmcs_seq_stats(lmcs_seq_info* stats, int32_t m_binNum);

void uvg_init_lmcs_aps(lmcs_aps* aps, int picWidth, int picHeight, uint32_t maxCUWidth, uint32_t maxCUHeight, int bitDepth);

void uvg_calc_seq_stats(struct encoder_state_t* const state, const videoframe_t* frame, lmcs_seq_info* stats, lmcs_aps* aps);

void uvg_lmcs_preanalyzer(struct encoder_state_t* const state, const videoframe_t* frame, lmcs_aps* aps, uint32_t signalType);

void uvg_construct_reshaper_lmcs(lmcs_aps* aps);

void uvg_encode_lmcs_adaptive_parameter_set(encoder_state_t* const state);

int uvg_calculate_lmcs_chroma_adj_vpdu_nei(encoder_state_t* const state, lmcs_aps* aps, int x, int y);