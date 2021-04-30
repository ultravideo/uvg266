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

/**
-Perform picture analysis for SDR
\param   pcPic describe pointer of current coding picture
\param   sliceType describe the slice type
\param   reshapeCW describe some input info
From VTM 12.1
*/

void kvz_calc_seq_stats(struct encoder_state_t* const state, const videoframe_t* frame, lmcs_seq_info* stats, lmcs_aps* aps)
{
  const encoder_control_t* const encoder = state->encoder_control;

  int32_t m_binNum = PIC_CODE_CW_BINS;
  kvz_pixel* picY = &frame->source->y[CU_TO_PIXEL(0, 0, 0, frame->source->stride)];
  const int width = frame->source->width;
  const int height = frame->source->height;
  const int stride = frame->source->stride;
  uint32_t winLens = (m_binNum == PIC_CODE_CW_BINS) ? (MIN(height, width) / 240) : 2;
  winLens = winLens > 0 ? winLens : 1;

  int64_t tempSq = 0;
  int64_t topSum = 0, topSumSq = 0;
  int64_t leftSum = 0, leftSumSq = 0;
  int64_t* leftColSum = malloc(sizeof(int64_t)*width);
  int64_t* leftColSumSq = malloc(sizeof(int64_t) * width);
  int64_t* topRowSum = malloc(sizeof(int64_t) * height);
  int64_t* topRowSumSq = malloc(sizeof(int64_t) * height);
  int64_t* topColSum = malloc(sizeof(int64_t) * width);
  int64_t* topColSumSq = malloc(sizeof(int64_t) * width);
  uint32_t* binCnt = malloc(sizeof(uint32_t)*m_binNum);
  memset(leftColSum, 0, width * sizeof(int64_t));
  memset(leftColSumSq, 0, width * sizeof(int64_t));
  memset(topRowSum, 0, height * sizeof(int64_t));
  memset(topRowSumSq, 0, height * sizeof(int64_t));
  memset(topColSum, 0, width * sizeof(int64_t));
  memset(topColSumSq, 0, width * sizeof(int64_t));
  memset(binCnt, 0, m_binNum * sizeof(uint32_t));

  kvz_init_lmcs_seq_stats(stats,m_binNum);
  for (uint32_t y = 0; y < height; y++)
  {
    for (uint32_t x = 0; x < width; x++)
    {
      const kvz_pixel pxlY = picY[x];
      int64_t sum = 0, sumSq = 0;
      uint32_t numPixInPart = 0;
      uint32_t y1 = MAX((int)(y - winLens), 0);
      uint32_t y2 = MIN((int)(y + winLens), (height - 1));
      uint32_t x1 = MAX((int)(x - winLens), 0);
      uint32_t x2 = MIN((int)(x + winLens), (width - 1));
      uint32_t bx = 0, by = 0;
      const kvz_pixel* pWinY = &picY[0];
      numPixInPart = (x2 - x1 + 1) * (y2 - y1 + 1);

      if (x == 0 && y == 0)
      {
        for (by = y1; by <= y2; by++)
        {
          for (bx = x1; bx <= x2; bx++)
          {
            tempSq = (int64_t)pWinY[bx] * (int64_t)pWinY[bx];
            leftSum += pWinY[bx];
            leftSumSq += tempSq;
            leftColSum[bx] += pWinY[bx];
            leftColSumSq[bx] += tempSq;
            topColSum[bx] += pWinY[bx];
            topColSumSq[bx] += tempSq;
            topRowSum[by] += pWinY[bx];
            topRowSumSq[by] += tempSq;
          }
          pWinY += stride;
        }
        topSum = leftSum;
        topSumSq = leftSumSq;
        sum = leftSum;
        sumSq = leftSumSq;
      }
      else if (x == 0 && y > 0)
      {
        if (y < height - winLens)
        {
          pWinY += winLens * stride;
          topRowSum[y + winLens] = 0;
          topRowSumSq[y + winLens] = 0;
          for (bx = x1; bx <= x2; bx++)
          {
            topRowSum[y + winLens] += pWinY[bx];
            topRowSumSq[y + winLens] += (int64_t)pWinY[bx] * (int64_t)pWinY[bx];
          }
          topSum += topRowSum[y + winLens];
          topSumSq += topRowSumSq[y + winLens];
        }
        if (y > winLens)
        {
          topSum -= topRowSum[y - 1 - winLens];
          topSumSq -= topRowSumSq[y - 1 - winLens];
        }
        memset(leftColSum, 0, width * sizeof(int64_t));
        memset(leftColSumSq, 0, width * sizeof(int64_t));
        pWinY = &picY[0];
        pWinY -= (y <= winLens ? y : winLens) * stride;
        for (by = y1; by <= y2; by++)
        {
          for (bx = x1; bx <= x2; bx++)
          {
            leftColSum[bx] += pWinY[bx];
            leftColSumSq[bx] += (int64_t)pWinY[bx] * (int64_t)pWinY[bx];
          }
          pWinY += stride;
        }
        leftSum = topSum;
        leftSumSq = topSumSq;
        sum = topSum;
        sumSq = topSumSq;
      }
      else if (x > 0)
      {
        if (x < width - winLens)
        {
          pWinY -= (y <= winLens ? y : winLens) * stride;
          if (y == 0)
          {
            leftColSum[x + winLens] = 0;
            leftColSumSq[x + winLens] = 0;
            for (by = y1; by <= y2; by++)
            {
              leftColSum[x + winLens] += pWinY[x + winLens];
              leftColSumSq[x + winLens] += (int64_t)pWinY[x + winLens] * (int64_t)pWinY[x + winLens];
              pWinY += stride;
            }
          }
          else
          {
            leftColSum[x + winLens] = topColSum[x + winLens];
            leftColSumSq[x + winLens] = topColSumSq[x + winLens];
            if (y < height - winLens)
            {
              pWinY = &picY[0];
              pWinY += winLens * stride;
              leftColSum[x + winLens] += pWinY[x + winLens];
              leftColSumSq[x + winLens] += (int64_t)pWinY[x + winLens] * (int64_t)pWinY[x + winLens];
            }
            if (y > winLens)
            {
              pWinY = &picY[0];
              pWinY -= (winLens + 1) * stride;
              leftColSum[x + winLens] -= pWinY[x + winLens];
              leftColSumSq[x + winLens] -= (int64_t)pWinY[x + winLens] * (int64_t)pWinY[x + winLens];
            }
          }
          topColSum[x + winLens] = leftColSum[x + winLens];
          topColSumSq[x + winLens] = leftColSumSq[x + winLens];
          leftSum += leftColSum[x + winLens];
          leftSumSq += leftColSumSq[x + winLens];
        }
        if (x > winLens)
        {
          leftSum -= leftColSum[x - 1 - winLens];
          leftSumSq -= leftColSumSq[x - 1 - winLens];
        }
        sum = leftSum;
        sumSq = leftSumSq;
      }

      double average = (double)(sum) / numPixInPart;
      double variance = (double)(sumSq) / numPixInPart - average * average;
      int binLen = aps->m_reshapeLUTSize / m_binNum;
      uint32_t binIdx = (uint32_t)(pxlY / binLen);
      if (aps->m_lumaBD > 10)
      {
        average = average / (double)(1 << (aps->m_lumaBD - 10));
        variance = variance / (double)(1 << (2 * aps->m_lumaBD - 20));
      }
      else if (aps->m_lumaBD < 10)
      {
        average = average * (double)(1 << (10 - aps->m_lumaBD));
        variance = variance * (double)(1 << (20 - 2 * aps->m_lumaBD));
      }
      double varLog10 = log10(variance + 1.0);
      stats->binVar[binIdx] += varLog10;
      binCnt[binIdx]++;
    }
    picY += stride;
  }

  for (int b = 0; b < m_binNum; b++)
  {
    stats->binHist[b] = (double)binCnt[b] / (double)(aps->m_reshapeCW.rspPicSize);
    stats->binVar[b] = (binCnt[b] > 0) ? (stats->binVar[b] / binCnt[b]) : 0.0;
  }
  FREE_POINTER(binCnt);
  FREE_POINTER(topColSum);
  FREE_POINTER(topColSumSq);
  FREE_POINTER(topRowSum);
  FREE_POINTER(topRowSumSq);
  FREE_POINTER(leftColSum);
  FREE_POINTER(leftColSumSq);

  stats->minBinVar = 5.0;
  stats->maxBinVar = 0.0;
  stats->meanBinVar = 0.0;
  stats->nonZeroCnt = 0;
  for (int b = 0; b < m_binNum; b++)
  {
    if (stats->binHist[b] > 0.001)
    {
      stats->nonZeroCnt++;
      stats->meanBinVar += stats->binVar[b];
      if (stats->binVar[b] > stats->maxBinVar) { stats->maxBinVar = stats->binVar[b]; }
      if (stats->binVar[b] < stats->minBinVar) { stats->minBinVar = stats->binVar[b]; }
    }
  }
  stats->meanBinVar /= (double)stats->nonZeroCnt;
  for (int b = 0; b < m_binNum; b++)
  {
    if (stats->meanBinVar > 0.0)
    {
      stats->normVar[b] = stats->binVar[b] / stats->meanBinVar;
    }
    stats->weightVar += stats->binHist[b] * stats->binVar[b];
    stats->weightNorm += stats->binHist[b] * stats->normVar[b];
  }

  picY = &frame->source->y[CU_TO_PIXEL(0, 0, 0, frame->source->stride)];
  double avgY = 0.0;
  double varY = 0.0;
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      avgY += picY[x];
      varY += (double)picY[x] * (double)picY[x];
    }
    picY += stride;
  }
  avgY = avgY / (width * height);
  varY = varY / (width * height) - avgY * avgY;

  if (encoder->chroma_format != KVZ_CSP_400)
  {
    // ToDo: Handle other than YUV 4:2:0
    assert(encoder->chroma_format == KVZ_CSP_420);

    kvz_pixel* picU = &frame->source->u[CU_TO_PIXEL(0, 0, 0, frame->source->stride/2)];
    kvz_pixel* picV = &frame->source->v[CU_TO_PIXEL(0, 0, 0, frame->source->stride / 2)];
    const int widthC = frame->source->width/2;
    const int heightC = frame->source->height/2;
    const int strideC = frame->source->stride/2;


    double avgU = 0.0, avgV = 0.0;
    double varU = 0.0, varV = 0.0;
    for (int y = 0; y < heightC; y++)
    {
      for (int x = 0; x < widthC; x++)
      {
        avgU += picU[x];
        avgV += picV[x];
        varU += (int64_t)picU[x] * (int64_t)picU[x];
        varV += (int64_t)picV[x] * (int64_t)picV[x];
      }
      picU += strideC;
      picV += strideC;
    }
    avgU = avgU / (widthC * heightC);
    avgV = avgV / (widthC * heightC);
    varU = varU / (widthC * heightC) - avgU * avgU;
    varV = varV / (widthC * heightC) - avgV * avgV;
    if (varY > 0)
    {
      stats->ratioStdU = sqrt(varU) / sqrt(varY);
      stats->ratioStdV = sqrt(varV) / sqrt(varY);
    }
  }
}