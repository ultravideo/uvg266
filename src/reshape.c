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

#include "reshape.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "cabac.h"
#include "rdo.h"
#include "strategies/strategies-sao.h"
#include "kvz_math.h"

void kvz_free_lmcs_aps(lmcs_aps* aps)
{
  FREE_POINTER(aps->m_invLUT);
  FREE_POINTER(aps->m_fwdLUT);

}

void kvz_init_lmcs_seq_stats(lmcs_seq_info* stats, int32_t m_binNum)
{
  for (int i = 0; i < m_binNum; i++)
  {
    stats->binVar[i] = 0.0;
    stats->binHist[i] = 0.0;
    stats->normVar[i] = 0.0;
  }
  stats->nonZeroCnt = 0;
  stats->weightVar = 0.0;
  stats->weightNorm = 0.0;
  stats->minBinVar = 0.0;
  stats->maxBinVar = 0.0;
  stats->meanBinVar = 0.0;
  stats->ratioStdU = 0.0;
  stats->ratioStdV = 0.0;
}

void kvz_init_lmcs_aps(lmcs_aps* aps, int picWidth, int picHeight, uint32_t maxCUWidth, uint32_t maxCUHeight, int bitDepth)
{
  aps->m_lumaBD = bitDepth;
  aps->m_reshapeLUTSize = 1 << aps->m_lumaBD;
  aps->m_initCWAnalyze = aps->m_reshapeLUTSize / PIC_ANALYZE_CW_BINS;
  aps->m_initCW = aps->m_reshapeLUTSize / PIC_CODE_CW_BINS;

  aps->m_fwdLUT = calloc(1, sizeof(kvz_pixel) * aps->m_reshapeLUTSize);

  aps->m_invLUT = calloc(1, sizeof(kvz_pixel) * aps->m_reshapeLUTSize);

  memset(aps->m_binCW, 0, sizeof(uint16_t) * PIC_ANALYZE_CW_BINS);
  memset(aps->m_binImportance, 0, sizeof(uint32_t) * PIC_ANALYZE_CW_BINS);
  memset(aps->m_reshapePivot, 0, sizeof(kvz_pixel) * PIC_CODE_CW_BINS + 1);
  memset(aps->m_inputPivot, 0, sizeof(kvz_pixel) * PIC_CODE_CW_BINS + 1);

  for (int i = 0; i < PIC_CODE_CW_BINS; i++) {
    aps->m_fwdScaleCoef[i] = 1 << FP_PREC;
  }

  for (int i = 0; i < PIC_CODE_CW_BINS; i++) {
    aps->m_invScaleCoef[i] = 1 << FP_PREC;
  }

  for (int i = 0; i < PIC_CODE_CW_BINS; i++) {
    aps->m_chromaAdjHelpLUT[i] = 1 << CSCALE_FP_PREC;
  }

  aps->m_sliceReshapeInfo.sliceReshaperEnableFlag = true;
  aps->m_sliceReshapeInfo.enableChromaAdj = true;
  aps->m_sliceReshapeInfo.sliceReshaperModelPresentFlag = true;
  aps->m_sliceReshapeInfo.reshaperModelMinBinIdx = 0;
  aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx = PIC_CODE_CW_BINS - 1;
  memset(aps->m_sliceReshapeInfo.reshaperModelBinCWDelta, 0, (PIC_CODE_CW_BINS) * sizeof(int));
  aps->m_sliceReshapeInfo.chrResScalingOffset = 0;

  aps->m_binNum = PIC_CODE_CW_BINS;
  kvz_init_lmcs_seq_stats(&aps->m_srcSeqStats, aps->m_binNum);
  kvz_init_lmcs_seq_stats(&aps->m_rspSeqStats, aps->m_binNum);
}