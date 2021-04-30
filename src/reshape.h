#pragma once
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2021 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * Kvazaar is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
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
#include "kvazaar.h"
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
  uint32_t* binCW;
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
  kvz_pixel*              m_invLUT;
  kvz_pixel*              m_fwdLUT;
  int32_t                 m_chromaAdjHelpLUT[PIC_CODE_CW_BINS];
  uint16_t                m_binCW[PIC_ANALYZE_CW_BINS];
  uint16_t                m_initCW;
  bool                    m_reshape;
  kvz_pixel               m_reshapePivot[PIC_CODE_CW_BINS + 1];
  kvz_pixel               m_inputPivot[PIC_CODE_CW_BINS + 1];
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
  kvz_pixel               m_cwLumaWeight[PIC_CODE_CW_BINS];
  double                  m_chromaWeight;
  int                     m_chromaAdj;
  int                     m_binNum;
  lmcs_seq_info           m_srcSeqStats;
  lmcs_seq_info           m_rspSeqStats;

} lmcs_aps;


void kvz_free_lmcs_aps(lmcs_aps* aps);

void kvz_init_lmcs_seq_stats(lmcs_seq_info* stats, int32_t m_binNum);

void kvz_init_lmcs_aps(lmcs_aps* aps, int picWidth, int picHeight, uint32_t maxCUWidth, uint32_t maxCUHeight, int bitDepth);

void kvz_calc_seq_stats(struct encoder_state_t* const state, const videoframe_t* frame, lmcs_seq_info* stats, lmcs_aps* aps);

void kvz_lmcs_preanalyzer(struct encoder_state_t* const state, const videoframe_t* frame, lmcs_seq_info* stats, lmcs_aps* aps, uint32_t signalType);