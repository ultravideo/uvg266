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
  //FREE_POINTER(aps->m_invLUT);
  //FREE_POINTER(aps->m_fwdLUT);

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

  //aps->m_fwdLUT = calloc(1, sizeof(kvz_pixel) * aps->m_reshapeLUTSize);
  //aps->m_invLUT = calloc(1, sizeof(kvz_pixel) * aps->m_reshapeLUTSize);
  assert(bitDepth <= 10); // ToDo: support bit depth larger than 10

  memset(aps->m_fwdLUT, 0, sizeof(kvz_pixel) * aps->m_reshapeLUTSize);
  memset(aps->m_invLUT, 0, sizeof(kvz_pixel) * aps->m_reshapeLUTSize);
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
  uint32_t winLens = (aps->m_binNum == PIC_CODE_CW_BINS) ? (MIN(height, width) / 240) : 2;
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

static void swap_int(int* xp, int* yp) { int temp = *xp;  *xp = *yp;  *yp = temp; }
static void swap_double(double* xp, double* yp) { double temp = *xp;  *xp = *yp;  *yp = temp; }

// Bubble Sort to  descending order with index
static void bubbleSortDsd(double* array, int* idx, int n)
{
  int i, j;
  bool swapped;
  for (i = 0; i < n - 1; i++)
  {
    swapped = false;
    for (j = 0; j < n - i - 1; j++)
    {
      if (array[j] < array[j + 1])
      {
        swap_double(&array[j], &array[j + 1]);
        swap_int(&idx[j], &idx[j + 1]);
        swapped = true;
      }
    }
    if (swapped == false)
    {
      break;
    }
  }
}

static void cwPerturbation(lmcs_aps* aps, int startBinIdx, int endBinIdx, uint16_t maxCW)
{
  for (int i = 0; i < aps->m_binNum; i++)
  {
    if (i >= startBinIdx && i <= endBinIdx)
    {
      aps->m_binCW[i] = (uint32_t)round((double)maxCW / (endBinIdx - startBinIdx + 1));
    }
    else
    {
      aps->m_binCW[i] = 0;
    }
  }

  double hist = 0.0;
  uint16_t delta1 = 0, delta2 = 0;
  for (int i = 0; i < aps->m_binNum; i++)
  {
    if (aps->m_srcSeqStats.binHist[i] > 0.001)
    {
      hist = aps->m_srcSeqStats.binHist[i] > 0.4 ? 0.4 : aps->m_srcSeqStats.binHist[i];
      delta1 = (uint16_t)(10.0 * hist + 0.5);
      delta2 = (uint16_t)(20.0 * hist + 0.5);
      if (aps->m_srcSeqStats.normVar[i] < 0.8)
      {
        aps->m_binCW[i] = aps->m_binCW[i] + delta2;
      }
      else if (aps->m_srcSeqStats.normVar[i] < 0.9)
      {
        aps->m_binCW[i] = aps->m_binCW[i] + delta1;
      }
      if (aps->m_srcSeqStats.normVar[i] > 1.2)
      {
        aps->m_binCW[i] = aps->m_binCW[i] - delta2;
      }
      else if (aps->m_srcSeqStats.normVar[i] > 1.1)
      {
        aps->m_binCW[i] = aps->m_binCW[i] - delta1;
      }
    }
  }
}

static void cwReduction(lmcs_aps* aps, int startBinIdx, int endBinIdx)
{
  int bdShift = aps->m_lumaBD - 10;
  int totCW = bdShift != 0 ? (bdShift > 0 ? aps->m_reshapeLUTSize / (1 << bdShift) : aps->m_reshapeLUTSize * (1 << (-bdShift))) : aps->m_reshapeLUTSize;
  int maxAllowedCW = totCW - 1, usedCW = 0;
  for (int i = 0; i < aps->m_binNum; i++)
  {
    usedCW += aps->m_binCW[i];
  }
  if (usedCW > maxAllowedCW)
  {
    int deltaCW = usedCW - maxAllowedCW;
    int divCW = deltaCW / (endBinIdx - startBinIdx + 1);
    int modCW = deltaCW - divCW * (endBinIdx - startBinIdx + 1);
    if (divCW > 0)
    {
      for (int i = startBinIdx; i <= endBinIdx; i++)
      {
        aps->m_binCW[i] -= divCW;
      }
    }
    for (int i = startBinIdx; i <= endBinIdx; i++)
    {
      if (modCW == 0)
      {
        break;
      }
      if (aps->m_binCW[i] > 0)
      {
        aps->m_binCW[i]--;
        modCW--;
      }
    }
  }
}

static void deriveReshapeParametersSDR(lmcs_aps* aps, bool* intraAdp, bool* interAdp)
{
  bool   isSkipCase = false;
  bool   isLowCase = false;
  int    firstBinVarLessThanVal1 = 0;
  int    firstBinVarLessThanVal2 = 0;
  int    firstBinVarLessThanVal3 = 0;
  double percBinVarLessThenVal1 = 0.0;
  double percBinVarLessThenVal2 = 0.0;
  double percBinVarLessThenVal3 = 0.0;
  assert(aps->m_binNum > 2); // For static analyzer to figure out we are not out-of-bounds
  assert(aps->m_binNum <= PIC_ANALYZE_CW_BINS);
  //int* binIdxSortDsd = malloc(sizeof(int) * aps->m_binNum);
  //double* binVarSortDsd = malloc(sizeof(double) * aps->m_binNum);
  //double* binVarSortDsdCDF = malloc(sizeof(double) * aps->m_binNum);
  int binIdxSortDsd[PIC_ANALYZE_CW_BINS];
  double binVarSortDsd[PIC_ANALYZE_CW_BINS];
  double binVarSortDsdCDF[PIC_ANALYZE_CW_BINS];

  double ratioWeiVar = 0.0, ratioWeiVarNorm = 0.0;
  int startBinIdx = aps->m_sliceReshapeInfo.reshaperModelMinBinIdx;
  int endBinIdx = aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx;

  for (int b = 0; b < aps->m_binNum; b++)
  {
    binVarSortDsd[b] = aps->m_srcSeqStats.binVar[b];
    binIdxSortDsd[b] = b;
  }
  bubbleSortDsd(binVarSortDsd, binIdxSortDsd, aps->m_binNum);
  binVarSortDsdCDF[0] = aps->m_srcSeqStats.binHist[binIdxSortDsd[0]];
  for (int b = 1; b < aps->m_binNum; b++)
  {
    binVarSortDsdCDF[b] = binVarSortDsdCDF[b - 1] + aps->m_srcSeqStats.binHist[binIdxSortDsd[b]];
  }
  for (int b = 0; b < aps->m_binNum - 1; b++)
  {
    if (binVarSortDsd[b] > 3.4)
    {
      firstBinVarLessThanVal1 = b + 1;
    }
    if (binVarSortDsd[b] > 2.8)
    {
      firstBinVarLessThanVal2 = b + 1;
    }
    if (binVarSortDsd[b] > 2.5)
    {
      firstBinVarLessThanVal3 = b + 1;
    }
  }
  percBinVarLessThenVal1 = binVarSortDsdCDF[firstBinVarLessThanVal1];
  percBinVarLessThenVal2 = binVarSortDsdCDF[firstBinVarLessThanVal2];
  percBinVarLessThenVal3 = binVarSortDsdCDF[firstBinVarLessThanVal3];
  //FREE_POINTER(binIdxSortDsd);
  //FREE_POINTER(binVarSortDsd);
  //FREE_POINTER(binVarSortDsdCDF);

  cwPerturbation(aps, startBinIdx, endBinIdx, (uint16_t)aps->m_reshapeCW.binCW[1]);
  cwReduction(aps, startBinIdx, endBinIdx);
  kvz_init_lmcs_seq_stats(&aps->m_rspSeqStats, aps->m_binNum);

  for (int b = 0; b < aps->m_binNum; b++)
  {
    double scale = (aps->m_binCW[b] > 0) ? ((double)aps->m_binCW[b] / (double)aps->m_initCWAnalyze) : 1.0;
    aps->m_rspSeqStats.binHist[b] = aps->m_srcSeqStats.binHist[b];
    aps->m_rspSeqStats.binVar[b] = aps->m_srcSeqStats.binVar[b] + 2.0 * log10(scale);
  }
  aps->m_rspSeqStats.minBinVar = 5.0;
  aps->m_rspSeqStats.maxBinVar = 0.0;
  aps->m_rspSeqStats.meanBinVar = 0.0;
  aps->m_rspSeqStats.nonZeroCnt = 0;
  for (int b = 0; b < aps->m_binNum; b++)
  {
    if (aps->m_rspSeqStats.binHist[b] > 0.001)
    {
      aps->m_rspSeqStats.nonZeroCnt++;
      aps->m_rspSeqStats.meanBinVar += aps->m_rspSeqStats.binVar[b];
      if (aps->m_rspSeqStats.binVar[b] > aps->m_rspSeqStats.maxBinVar)
      {
        aps->m_rspSeqStats.maxBinVar = aps->m_rspSeqStats.binVar[b];
      }
      if (aps->m_rspSeqStats.binVar[b] < aps->m_rspSeqStats.minBinVar)
      {
        aps->m_rspSeqStats.minBinVar = aps->m_rspSeqStats.binVar[b];
      }
    }
  }
  aps->m_rspSeqStats.meanBinVar /= (double)aps->m_rspSeqStats.nonZeroCnt;
  for (int b = 0; b < aps->m_binNum; b++)
  {
    if (aps->m_rspSeqStats.meanBinVar > 0.0)
    {
      aps->m_rspSeqStats.normVar[b] = aps->m_rspSeqStats.binVar[b] / aps->m_rspSeqStats.meanBinVar;
    }
    aps->m_rspSeqStats.weightVar += aps->m_rspSeqStats.binHist[b] * aps->m_rspSeqStats.binVar[b];
    aps->m_rspSeqStats.weightNorm += aps->m_rspSeqStats.binHist[b] * aps->m_rspSeqStats.normVar[b];
  }
  ratioWeiVar = aps->m_rspSeqStats.weightVar / aps->m_srcSeqStats.weightVar;
  ratioWeiVarNorm = aps->m_rspSeqStats.weightNorm / aps->m_srcSeqStats.weightNorm;

  if ((aps->m_srcSeqStats.binHist[0] + aps->m_srcSeqStats.binHist[aps->m_binNum - 1]) > 0.0001 && aps->m_srcSeqStats.binHist[aps->m_binNum - 2] < 0.001)
  {
    if (percBinVarLessThenVal3 > 0.8 && percBinVarLessThenVal2 > 0.4 && aps->m_srcSeqStats.binVar[aps->m_binNum - 2] > 4.8)
    {
      isSkipCase = true;
    }
    else if (percBinVarLessThenVal3 < 0.1 && percBinVarLessThenVal1 < 0.05 && aps->m_srcSeqStats.binVar[aps->m_binNum - 2] < 4.0)
    {
      isSkipCase = true;
    }
  }
  if (isSkipCase)
  {
    *intraAdp = false;
    *interAdp = false;
    return;
  }

  if (aps->m_reshapeCW.rspPicSize > 5184000)
  {
    isLowCase = true;
  }
  else if (aps->m_srcSeqStats.binVar[1] > 4.0)
  {
    isLowCase = true;
  }
  else if (aps->m_rspSeqStats.meanBinVar > 3.4 && ratioWeiVarNorm > 1.005 && ratioWeiVar > 1.02)
  {
    isLowCase = true;
  }
  else if (aps->m_rspSeqStats.meanBinVar > 3.1 && ratioWeiVarNorm > 1.005 && ratioWeiVar > 1.04)
  {
    isLowCase = true;
  }
  else if (aps->m_rspSeqStats.meanBinVar > 2.8 && ratioWeiVarNorm > 1.01 && ratioWeiVar > 1.04)
  {
    isLowCase = true;
  }

  if (aps->m_reshapeCW.updateCtrl == 0)
  {
    aps->m_reshapeCW.binCW[1] = 1022;
    if (isLowCase)
    {
      *intraAdp = false;
      aps->m_rateAdpMode = 1;
      aps->m_reshapeCW.binCW[1] = 980;
      if (aps->m_srcSeqStats.binHist[aps->m_binNum - 2] > 0.05)
      {
        aps->m_reshapeCW.binCW[1] = 896;
        if (aps->m_srcSeqStats.binVar[aps->m_binNum - 2] < 1.2)
        {
          aps->m_reshapeCW.binCW[1] = 938;
        }
      }
      else if (percBinVarLessThenVal2 < 0.8 && percBinVarLessThenVal3 == 1.0)
      {
        aps->m_rateAdpMode = 1;
        aps->m_reshapeCW.binCW[1] = 938;
      }
    }
    if (aps->m_srcSeqStats.binHist[aps->m_binNum - 2] < 0.001)
    {
      if (aps->m_srcSeqStats.binHist[1] > 0.05 && aps->m_srcSeqStats.binVar[1] > 3.0)
      {
        *intraAdp = true;
        aps->m_rateAdpMode = 1;
        aps->m_reshapeCW.binCW[1] = 784;
      }
      else if (aps->m_srcSeqStats.binHist[1] < 0.006)
      {
        *intraAdp = false;
        aps->m_rateAdpMode = 0;
        aps->m_reshapeCW.binCW[1] = 1008;
      }
      else if (percBinVarLessThenVal3 < 0.5)
      {
        *intraAdp = true;
        aps->m_rateAdpMode = 0;
        aps->m_reshapeCW.binCW[1] = 1022;
      }
    }
    else if ((aps->m_srcSeqStats.maxBinVar > 4.0 && aps->m_rspSeqStats.meanBinVar > 3.2 && percBinVarLessThenVal2 < 0.25) || ratioWeiVar < 1.03)
    {
      *intraAdp = true;
      aps->m_rateAdpMode = 0;
      aps->m_reshapeCW.binCW[1] = 1022;
    }
    if (*intraAdp == true && aps->m_rateAdpMode == 0)
    {
      aps->m_tcase = 9;
    }
  }
  else if (aps->m_reshapeCW.updateCtrl == 1)
  {
    aps->m_reshapeCW.binCW[1] = 952;
    if (isLowCase)
    {
      if (aps->m_reshapeCW.rspPicSize > 5184000)
      {
        aps->m_rateAdpMode = 1;
        aps->m_reshapeCW.binCW[1] = 812;
      }
      if (aps->m_srcSeqStats.binHist[aps->m_binNum - 2] > 0.05)
      {
        aps->m_rateAdpMode = 1;
        aps->m_reshapeCW.binCW[1] = 812;
        if (aps->m_srcSeqStats.binHist[aps->m_binNum - 2] > 0.1 || aps->m_srcSeqStats.binHist[1] > 0.1)
        {
          aps->m_rateAdpMode = 0;
          aps->m_reshapeCW.binCW[1] = 924;
        }
      }
      else if (percBinVarLessThenVal2 < 0.8 && percBinVarLessThenVal3 == 1.0)
      {
        aps->m_rateAdpMode = 1;
        aps->m_reshapeCW.binCW[1] = 896;
      }
      else if (percBinVarLessThenVal2 > 0.98 && aps->m_srcSeqStats.binHist[1] > 0.05)
      {
        aps->m_rateAdpMode = 0;
        aps->m_reshapeCW.binCW[1] = 784;
      }
      else if (percBinVarLessThenVal2 < 0.1)
      {
        aps->m_rateAdpMode = 0;
        aps->m_reshapeCW.binCW[1] = 1022;
      }
    }
    if (aps->m_srcSeqStats.binHist[1] > 0.1 && (aps->m_srcSeqStats.binVar[1] > 1.8 && aps->m_srcSeqStats.binVar[1] < 3.0))
    {
      aps->m_rateAdpMode = 1;
      if (aps->m_srcSeqStats.binVar[aps->m_binNum - 2] > 1.2 && aps->m_srcSeqStats.binVar[aps->m_binNum - 2] < 4.0)
      {
        aps->m_reshapeCW.binCW[1] = 784;
      }
    }
    else if (aps->m_srcSeqStats.binHist[aps->m_binNum - 2] < 0.001)
    {
      if (aps->m_srcSeqStats.binHist[1] > 0.05 && aps->m_srcSeqStats.binVar[1] > 3.0)
      {
        aps->m_rateAdpMode = 1;
        aps->m_reshapeCW.binCW[1] = 784;
      }
      else if (aps->m_srcSeqStats.binHist[1] < 0.006)
      {
        aps->m_rateAdpMode = 0;
        aps->m_reshapeCW.binCW[1] = 980;
      }
      else if (percBinVarLessThenVal3 < 0.5)
      {
        aps->m_rateAdpMode = 0;
        aps->m_reshapeCW.binCW[1] = 924;
      }
    }
    else if ((aps->m_srcSeqStats.maxBinVar > 4.0 && aps->m_rspSeqStats.meanBinVar > 3.2 && percBinVarLessThenVal2 < 0.25) || ratioWeiVar < 1.03)
    {
      aps->m_rateAdpMode = 0;
      aps->m_reshapeCW.binCW[1] = 980;
    }
  }
  else
  {
    aps->m_useAdpCW = true;
    aps->m_reshapeCW.binCW[0] = 36;  aps->m_reshapeCW.binCW[1] = 30;
    if (isLowCase)
    {
      if (aps->m_srcSeqStats.binHist[aps->m_binNum - 2] > 0.05)
      {
        aps->m_useAdpCW = false;
        aps->m_rateAdpMode = 1;
        aps->m_reshapeCW.binCW[1] = 896;
        if (aps->m_srcSeqStats.binHist[1] > 0.005)
        {
          aps->m_rateAdpMode = 0;
        }
      }
      else if (percBinVarLessThenVal2 < 0.8 && percBinVarLessThenVal3 == 1.0)
      {
        aps->m_reshapeCW.binCW[1] = 28;
      }
    }
    if (aps->m_srcSeqStats.binHist[1] > 0.1 && aps->m_srcSeqStats.binVar[1] > 1.8 && aps->m_srcSeqStats.binVar[1] < 3.0)
    {
      aps->m_useAdpCW = false;
      aps->m_rateAdpMode = 1;
      aps->m_reshapeCW.binCW[1] = 952;
    }
    else if (aps->m_srcSeqStats.binHist[1] > 0.05 && aps->m_srcSeqStats.binHist[aps->m_binNum - 2] < 0.001 && aps->m_srcSeqStats.binVar[1] > 3.0)
    {
      aps->m_useAdpCW = false;
      aps->m_rateAdpMode = 1;
      aps->m_reshapeCW.binCW[1] = 784;
    }
    else if (aps->m_srcSeqStats.binHist[1] > 0.05 && aps->m_srcSeqStats.binHist[aps->m_binNum - 2] < 0.005 && aps->m_srcSeqStats.binVar[1] > 1.0 && aps->m_srcSeqStats.binVar[1] < 1.5)
    {
      aps->m_rateAdpMode = 2;
      aps->m_reshapeCW.binCW[0] = 38;
    }
    else if (aps->m_srcSeqStats.binHist[1] < 0.005 && aps->m_srcSeqStats.binHist[aps->m_binNum - 2] > 0.05 && aps->m_srcSeqStats.binVar[aps->m_binNum - 2] > 1.0 && aps->m_srcSeqStats.binVar[aps->m_binNum - 2] < 1.5)
    {
      aps->m_rateAdpMode = 2;
      aps->m_reshapeCW.binCW[0] = 36;
    }
    else if (aps->m_srcSeqStats.binHist[1] > 0.02 && aps->m_srcSeqStats.binHist[aps->m_binNum - 2] > 0.04 && aps->m_srcSeqStats.binVar[1] < 2.0 && aps->m_srcSeqStats.binVar[aps->m_binNum - 2] < 1.5)
    {
      aps->m_rateAdpMode = 2;
      aps->m_reshapeCW.binCW[0] = 34;
    }
    else if ((aps->m_srcSeqStats.binHist[1] > 0.05 && aps->m_srcSeqStats.binHist[aps->m_binNum - 2] > 0.2 && aps->m_srcSeqStats.binVar[1] > 3.0 && aps->m_srcSeqStats.binVar[1] < 4.0) || ratioWeiVar < 1.03)
    {
      aps->m_rateAdpMode = 1;
      aps->m_reshapeCW.binCW[0] = 34;
    }
    else if (aps->m_srcSeqStats.binVar[1] < 4.0 && percBinVarLessThenVal2 == 1.0 && percBinVarLessThenVal3 == 1.0)
    {
      aps->m_rateAdpMode = 0;
      aps->m_reshapeCW.binCW[0] = 34;
    }
    if (aps->m_useAdpCW && !isLowCase)
    {
      aps->m_reshapeCW.binCW[1] = 66 - aps->m_reshapeCW.binCW[0];
    }
  }
}

static void deriveReshapeParameters(double* array, int start, int end, ReshapeCW* respCW, double* alpha, double* beta)
{
  double minVar = 10.0, maxVar = 0.0;
  for (int b = start; b <= end; b++)
  {
    if (array[b] < minVar)
    {
      minVar = array[b];
    }
    if (array[b] > maxVar)
    {
      maxVar = array[b];
    }
  }
  double maxCW = (double)respCW->binCW[0];
  double minCW = (double)respCW->binCW[1];
  *alpha = (minCW - maxCW) / (maxVar - minVar);
  *beta = (maxCW * maxVar - minCW * minVar) / (maxVar - minVar);
}

void kvz_lmcs_preanalyzer(struct encoder_state_t* const state, const videoframe_t* frame, lmcs_aps* aps, uint32_t signalType)
{

  enum kvz_slice_type sliceType = state->frame->slicetype;
  aps->m_sliceReshapeInfo.sliceReshaperModelPresentFlag = true;
  aps->m_sliceReshapeInfo.sliceReshaperEnableFlag = true;

  int modIP = state->frame->poc - state->frame->poc / aps->m_reshapeCW.rspFpsToIp * aps->m_reshapeCW.rspFpsToIp;
  if (sliceType == KVZ_SLICE_I || (aps->m_reshapeCW.updateCtrl == 2 && modIP == 0))
  {
    if (aps->m_sliceReshapeInfo.sliceReshaperModelPresentFlag == true)
    {
      aps->m_binNum = PIC_CODE_CW_BINS;
      int stdMin = 16 << (aps->m_lumaBD - 8);
      int stdMax = 235 << (aps->m_lumaBD - 8);
      int binLen = aps->m_reshapeLUTSize / aps->m_binNum;
      int startBinIdx = stdMin / binLen;
      int endBinIdx = stdMax / binLen;
      aps->m_sliceReshapeInfo.reshaperModelMinBinIdx = startBinIdx;
      aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx = endBinIdx;
      aps->m_initCWAnalyze = aps->m_lumaBD > 10 ? (binLen >> (aps->m_lumaBD - 10)) : aps->m_lumaBD < 10 ? (binLen << (10 - aps->m_lumaBD)) : binLen;
      for (int b = 0; b < aps->m_binNum; b++)
      {
        aps->m_binCW[b] = aps->m_initCWAnalyze;
      }

      aps->m_reshape = true;
      aps->m_useAdpCW = false;
      aps->m_exceedSTD = false;
      aps->m_chromaWeight = 1.0;
      aps->m_sliceReshapeInfo.enableChromaAdj = 1;
      aps->m_rateAdpMode = 0;  aps->m_tcase = 0;
      bool intraAdp = true, interAdp = true;

      kvz_calc_seq_stats(state, frame, &aps->m_srcSeqStats, aps);

      bool isFlat = true;
      for (int b = 0; b < aps->m_binNum; b++)
      {
        if (aps->m_srcSeqStats.binVar[b] > 0)
        {
          isFlat = false;
        }
      }
      if (isFlat)
      {
        intraAdp = false;
        interAdp = false;
      }
      if (aps->m_binNum == PIC_CODE_CW_BINS)
      {
        if ((aps->m_srcSeqStats.binHist[0] + aps->m_srcSeqStats.binHist[aps->m_binNum - 1]) > 0.005)
        {
          aps->m_exceedSTD = true;
        }
        if (aps->m_srcSeqStats.binHist[aps->m_binNum - 1] > 0.0003)
        {
          intraAdp = false;
          interAdp = false;
        }
        if (aps->m_srcSeqStats.binHist[0] > 0.03)
        {
          intraAdp = false;
          interAdp = false;
        }
      }
      else if (aps->m_binNum == PIC_ANALYZE_CW_BINS)
      {
        if ((aps->m_srcSeqStats.binHist[0] + aps->m_srcSeqStats.binHist[1] + aps->m_srcSeqStats.binHist[aps->m_binNum - 2]
          + aps->m_srcSeqStats.binHist[aps->m_binNum - 1])
  > 0.01)
        {
          aps->m_exceedSTD = true;
        }
        if ((aps->m_srcSeqStats.binHist[aps->m_binNum - 2] + aps->m_srcSeqStats.binHist[aps->m_binNum - 1]) > 0.0003)
        {
          intraAdp = false;
          interAdp = false;
        }
        if ((aps->m_srcSeqStats.binHist[0] + aps->m_srcSeqStats.binHist[1]) > 0.03)
        {
          intraAdp = false;
          interAdp = false;
        }
      }
      if (aps->m_exceedSTD)
      {
        for (int i = 0; i < aps->m_binNum; i++)
        {
          if (aps->m_srcSeqStats.binHist[i] > 0 && i < startBinIdx)
          {
            startBinIdx = i;
          }
          if (aps->m_srcSeqStats.binHist[i] > 0 && i > endBinIdx)
          {
            endBinIdx = i;
          }
        }
        aps->m_sliceReshapeInfo.reshaperModelMinBinIdx = startBinIdx;
        aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx = endBinIdx;
      }

      if ((aps->m_srcSeqStats.ratioStdU + aps->m_srcSeqStats.ratioStdV) > 1.5 && aps->m_srcSeqStats.binHist[1] > 0.5)
      {
        intraAdp = false;
        interAdp = false;
      }
      if (aps->m_srcSeqStats.ratioStdU > 0.36 && aps->m_srcSeqStats.ratioStdV > 0.2 && aps->m_reshapeCW.rspPicSize > 5184000)
      {
        aps->m_sliceReshapeInfo.enableChromaAdj = 0; aps->m_chromaWeight = 1.05;
        if ((aps->m_srcSeqStats.ratioStdU + aps->m_srcSeqStats.ratioStdV) < 0.69)
        {
          aps->m_chromaWeight = 0.95;
        }
      }

      if (interAdp)
      {
        if (aps->m_reshapeCW.adpOption)
        {
          aps->m_reshapeCW.binCW[0] = 0; aps->m_reshapeCW.binCW[1] = aps->m_reshapeCW.initialCW;
          aps->m_rateAdpMode = aps->m_reshapeCW.adpOption - 2 * (aps->m_reshapeCW.adpOption / 2);
          if (aps->m_reshapeCW.adpOption == 2)
          {
            aps->m_tcase = 9;
          }
          else if (aps->m_reshapeCW.adpOption > 2)
          {
            intraAdp = false;
          }
        }
        else if (signalType == RESHAPE_SIGNAL_SDR)
        {
          aps->m_reshapeCW.binCW[0] = 0; aps->m_reshapeCW.binCW[1] = 1022;
          deriveReshapeParametersSDR(aps, &intraAdp, &interAdp);
        }
        else if (signalType == RESHAPE_SIGNAL_HLG)
        {
          if (aps->m_reshapeCW.updateCtrl == 0)
          {
            aps->m_rateAdpMode = 0;  aps->m_tcase = 9;
            aps->m_reshapeCW.binCW[1] = 952;
            if (aps->m_srcSeqStats.meanBinVar < 2.5)
            {
              aps->m_reshapeCW.binCW[1] = 840;
            }
          }
          else
          {
            aps->m_useAdpCW = true;
            aps->m_rateAdpMode = 2;
            if (aps->m_binNum == PIC_CODE_CW_BINS)
            {
              aps->m_reshapeCW.binCW[0] = 72;
              aps->m_reshapeCW.binCW[1] = 58;
            }
            else if (aps->m_binNum == PIC_ANALYZE_CW_BINS)
            {
              aps->m_reshapeCW.binCW[0] = 36;
              aps->m_reshapeCW.binCW[1] = 30;
            }
            if (aps->m_srcSeqStats.meanBinVar < 2.5)
            {
              intraAdp = false;
              interAdp = false;
            }
          }
        }
      }

      if (aps->m_rateAdpMode == 2 && aps->m_reshapeCW.rspBaseQP <= 22)
      {
        intraAdp = false;
        interAdp = false;
      }
      aps->m_sliceReshapeInfo.sliceReshaperEnableFlag = intraAdp;
      if (!intraAdp && !interAdp)
      {
        aps->m_sliceReshapeInfo.sliceReshaperModelPresentFlag = false;
        aps->m_reshape = false;
        return;
      }

      if (aps->m_rateAdpMode == 1 && aps->m_reshapeCW.rspBaseQP <= 22)
      {
        for (int i = 0; i < aps->m_binNum; i++)
        {
          if (i >= startBinIdx && i <= endBinIdx)
          {
            aps->m_binCW[i] = aps->m_initCWAnalyze + 2;
          }
          else
          {
            aps->m_binCW[i] = 0;
          }
        }
      }
      else if (aps->m_useAdpCW)
      {
        if (signalType == RESHAPE_SIGNAL_SDR && aps->m_reshapeCW.updateCtrl == 2)
        {
          aps->m_binNum = PIC_ANALYZE_CW_BINS;
          startBinIdx = startBinIdx * 2;
          endBinIdx = endBinIdx * 2 + 1;
          kvz_calc_seq_stats(state, frame, &aps->m_srcSeqStats, aps);
        }
        double alpha = 1.0, beta = 0.0;
        deriveReshapeParameters(aps->m_srcSeqStats.binVar, startBinIdx, endBinIdx, &aps->m_reshapeCW, &alpha, &beta);
        for (int i = 0; i < aps->m_binNum; i++)
        {
          if (i >= startBinIdx && i <= endBinIdx)
          {
            aps->m_binCW[i] = (uint32_t)round(alpha * aps->m_srcSeqStats.binVar[i] + beta);
          }
          else
          {
            aps->m_binCW[i] = 0;
          }
        }
      }
      else
      {
        cwPerturbation(aps, startBinIdx, endBinIdx, (uint16_t)aps->m_reshapeCW.binCW[1]);
      }
      cwReduction(aps, startBinIdx, endBinIdx);
    }
    aps->m_chromaAdj = aps->m_sliceReshapeInfo.enableChromaAdj;
  }
  else // Inter slices
  {
    aps->m_sliceReshapeInfo.sliceReshaperModelPresentFlag = false;
    aps->m_sliceReshapeInfo.enableChromaAdj = aps->m_chromaAdj;
    if (!aps->m_reshape)
    {
      aps->m_sliceReshapeInfo.sliceReshaperEnableFlag = false;
    }
    else
    {
      const int cTid = aps->m_reshapeCW.rspTid;
      bool enableRsp = aps->m_tcase == 5 ? false : (aps->m_tcase < 5 ? (cTid < aps->m_tcase + 1 ? false : true) : (cTid <= 10 - aps->m_tcase ? true : false));
      aps->m_sliceReshapeInfo.sliceReshaperEnableFlag = enableRsp;
    }
  }
}


static void adjust_lmcs_pivot(lmcs_aps* aps)
{
  int bdShift = aps->m_lumaBD - 10;
  int totCW = bdShift != 0 ? (bdShift > 0 ? aps->m_reshapeLUTSize / (1 << bdShift) : aps->m_reshapeLUTSize * (1 << (-bdShift))) : aps->m_reshapeLUTSize;
  int orgCW = totCW / PIC_CODE_CW_BINS;
  int log2SegSize = aps->m_lumaBD - 5;//kvz_math_floor_log2(LMCS_SEG_NUM);

  aps->m_reshapePivot[0] = 0;
  for (int i = 0; i < PIC_CODE_CW_BINS; i++)
  {
    aps->m_reshapePivot[i + 1] = aps->m_reshapePivot[i] + aps->m_binCW[i];
  }
  int segIdxMax = (aps->m_reshapePivot[aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx + 1] >> log2SegSize);
  for (int i = aps->m_sliceReshapeInfo.reshaperModelMinBinIdx; i <= aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx; i++)
  {
    aps->m_reshapePivot[i + 1] = aps->m_reshapePivot[i] + aps->m_binCW[i];
    int segIdxCurr = (aps->m_reshapePivot[i] >> log2SegSize);
    int segIdxNext = (aps->m_reshapePivot[i + 1] >> log2SegSize);

    if ((segIdxCurr == segIdxNext) && (aps->m_reshapePivot[i] != (segIdxCurr << log2SegSize)))
    {
      if (segIdxCurr == segIdxMax)
      {
        aps->m_reshapePivot[i] = aps->m_reshapePivot[aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx + 1];
        for (int j = i; j <= aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx; j++)
        {
          aps->m_reshapePivot[j + 1] = aps->m_reshapePivot[i];
          aps->m_binCW[j] = 0;
        }
        aps->m_binCW[i - 1] = aps->m_reshapePivot[i] - aps->m_reshapePivot[i - 1];
        break;
      }
      else
      {
        int16_t adjustVal = ((segIdxCurr + 1) << log2SegSize) - aps->m_reshapePivot[i + 1];
        aps->m_reshapePivot[i + 1] += adjustVal;
        aps->m_binCW[i] += adjustVal;

        for (int j = i + 1; j <= aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx; j++)
        {
          if (aps->m_binCW[j] < (adjustVal + (orgCW >> 3)))
          {
            adjustVal -= (aps->m_binCW[j] - (orgCW >> 3));
            aps->m_binCW[j] = (orgCW >> 3);
          }
          else
          {
            aps->m_binCW[j] -= adjustVal;
            adjustVal = 0;
          }
          if (adjustVal == 0)
          {
            break;
          }
        }
      }
    }
  }

  for (int i = PIC_CODE_CW_BINS - 1; i >= 0; i--)
  {
    if (aps->m_binCW[i] > 0)
    {
      aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx = i;
      break;
    }
  }
}

static int get_pwl_idx_inv(lmcs_aps* aps,int lumaVal)
{
  int idxS = 0;
  for (idxS = aps->m_sliceReshapeInfo.reshaperModelMinBinIdx; (idxS <= aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx); idxS++)
  {
    if (lumaVal < aps->m_reshapePivot[idxS + 1])     break;
  }
  return MIN(idxS, PIC_CODE_CW_BINS - 1);
}

void kvz_construct_reshaper_lmcs(lmcs_aps* aps)
{
  int bdShift = aps->m_lumaBD - 10;
  int totCW = bdShift != 0 ? (bdShift > 0 ? aps->m_reshapeLUTSize / (1 << bdShift) : aps->m_reshapeLUTSize * (1 << (-bdShift))) : aps->m_reshapeLUTSize;
  int histLenth = totCW / aps->m_binNum;
  int log2HistLenth = kvz_math_floor_log2(histLenth);
  int i;

  if (aps->m_binNum == PIC_ANALYZE_CW_BINS)
  {
    for (int i = 0; i < PIC_CODE_CW_BINS; i++)
    {
      aps->m_binCW[i] = aps->m_binCW[2 * i] + aps->m_binCW[2 * i + 1];
    }
  }
  for (int i = 0; i <= PIC_CODE_CW_BINS; i++)
  {
    aps->m_inputPivot[i] = aps->m_initCW * i;
  }

  aps->m_sliceReshapeInfo.reshaperModelMinBinIdx = 0;
  aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx = PIC_CODE_CW_BINS - 1;
  for (int i = 0; i < PIC_CODE_CW_BINS; i++)
  {
    if (aps->m_binCW[i] > 0)
    {
      aps->m_sliceReshapeInfo.reshaperModelMinBinIdx = i;
      break;
    }
  }
  for (int i = PIC_CODE_CW_BINS - 1; i >= 0; i--)
  {
    if (aps->m_binCW[i] > 0)
    {
      aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx = i;
      break;
    }
  }

  if (bdShift != 0)
  {
    for (int i = 0; i < PIC_ANALYZE_CW_BINS; i++)
    {
      aps->m_binCW[i] = bdShift > 0 ? aps->m_binCW[i] * (1 << bdShift) : aps->m_binCW[i] / (1 << (-bdShift));
    }
  }

  adjust_lmcs_pivot(aps);

  int maxAbsDeltaCW = 0, absDeltaCW = 0, deltaCW = 0;
  for (int i = aps->m_sliceReshapeInfo.reshaperModelMinBinIdx; i <= aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx; i++)
  {
    deltaCW = (int)aps->m_binCW[i] - (int)aps->m_initCW;
    aps->m_sliceReshapeInfo.reshaperModelBinCWDelta[i] = deltaCW;
    absDeltaCW = (deltaCW < 0) ? (-deltaCW) : deltaCW;
    if (absDeltaCW > maxAbsDeltaCW)
    {
      maxAbsDeltaCW = absDeltaCW;
    }
  }
  aps->m_sliceReshapeInfo.maxNbitsNeededDeltaCW = MAX(1, 1 + kvz_math_floor_log2(maxAbsDeltaCW));

  histLenth = aps->m_initCW;
  log2HistLenth = kvz_math_floor_log2(histLenth);

  int sumBins = 0;
  for (i = 0; i < PIC_CODE_CW_BINS; i++) { sumBins += aps->m_binCW[i]; }
  assert(sumBins < aps->m_reshapeLUTSize && "SDR CW assignment is wrong!!");
  for (int i = 0; i < PIC_CODE_CW_BINS; i++)
  {
    aps->m_reshapePivot[i + 1] = aps->m_reshapePivot[i] + aps->m_binCW[i];
    aps->m_fwdScaleCoef[i] = ((int32_t)aps->m_binCW[i] * (1 << FP_PREC) + (1 << (log2HistLenth - 1))) >> log2HistLenth;
    if (aps->m_binCW[i] == 0)
    {
      aps->m_invScaleCoef[i] = 0;
      aps->m_chromaAdjHelpLUT[i] = 1 << CSCALE_FP_PREC;
    }
    else
    {
      aps->m_invScaleCoef[i] = (int32_t)(aps->m_initCW * (1 << FP_PREC) / aps->m_binCW[i]);
      aps->m_chromaAdjHelpLUT[i] = (int32_t)(aps->m_initCW * (1 << FP_PREC) / (aps->m_binCW[i] + aps->m_sliceReshapeInfo.chrResScalingOffset));
    }
  }
  for (int lumaSample = 0; lumaSample < aps->m_reshapeLUTSize; lumaSample++)
  {
    int idxY = lumaSample / aps->m_initCW;
    int tempVal = aps->m_reshapePivot[idxY] + ((aps->m_fwdScaleCoef[idxY] * (lumaSample - aps->m_inputPivot[idxY]) + (1 << (FP_PREC - 1))) >> FP_PREC);
    aps->m_fwdLUT[lumaSample] = CLIP((kvz_pixel)0, (kvz_pixel)((1 << aps->m_lumaBD) - 1), (kvz_pixel)(tempVal));

    int idxYInv = get_pwl_idx_inv(aps,lumaSample);
    int invSample = aps->m_inputPivot[idxYInv] + ((aps->m_invScaleCoef[idxYInv] * (lumaSample - aps->m_reshapePivot[idxYInv]) + (1 << (FP_PREC - 1))) >> FP_PREC);
    aps->m_invLUT[lumaSample] = CLIP((kvz_pixel)0, (kvz_pixel)((1 << aps->m_lumaBD) - 1), (kvz_pixel)(invSample));
  }
  for (i = 0; i < PIC_CODE_CW_BINS; i++)
  {
    int start = i * histLenth;
    int end = (i + 1) * histLenth - 1;
    aps->m_cwLumaWeight[i] = aps->m_fwdLUT[end] - aps->m_fwdLUT[start];
  }
}


static void code_lmcs_aps(encoder_state_t* const state, lmcs_aps* aps)
{
  bitstream_t* const stream = &state->stream;

  WRITE_UE(stream, aps->m_sliceReshapeInfo.reshaperModelMinBinIdx, "lmcs_min_bin_idx");
  WRITE_UE(stream, PIC_CODE_CW_BINS - 1 - aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx, "lmcs_delta_max_bin_idx");
  assert(aps->m_sliceReshapeInfo.maxNbitsNeededDeltaCW > 0);
  WRITE_UE(stream, aps->m_sliceReshapeInfo.maxNbitsNeededDeltaCW - 1, "lmcs_delta_cw_prec_minus1");

  for (int i = aps->m_sliceReshapeInfo.reshaperModelMinBinIdx; i <= aps->m_sliceReshapeInfo.reshaperModelMaxBinIdx; i++)
  {
    int deltaCW = aps->m_sliceReshapeInfo.reshaperModelBinCWDelta[i];
    int signCW = (deltaCW < 0) ? 1 : 0;
    int absCW = (deltaCW < 0) ? (-deltaCW) : deltaCW;
    WRITE_U(stream, absCW, aps->m_sliceReshapeInfo.maxNbitsNeededDeltaCW, "lmcs_delta_abs_cw[ i ]");
    if (absCW > 0)
    {
      WRITE_U(stream, signCW, 1, "lmcs_delta_sign_cw_flag[ i ]");
    }
  }
  // ToDo: LMCS Chroma scaling
  /*  
  int deltaCRS = aps-> ->chromaPresentFlag ? aps->m_sliceReshapeInfo.chrResScalingOffset : 0;
  int signCRS = (deltaCRS < 0) ? 1 : 0;
  int absCRS = (deltaCRS < 0) ? (-deltaCRS) : deltaCRS;
  if (pcAPS->chromaPresentFlag)
  {
    WRITE_CODE(absCRS, 3, "lmcs_delta_abs_crs");
  }
  if (absCRS > 0)
  {
    WRITE_FLAG(signCRS, "lmcs_delta_sign_crs_flag");
  }
  */
}



void kvz_encode_lmcs_adaptive_parameter_set(encoder_state_t* const state)
{
  bitstream_t* const stream = &state->stream;

  if (state->encoder_control->cfg.lmcs_enable) {
    // ToDo: Write LMCS APS NAL
    
    kvz_nal_write(stream, NAL_UNIT_PREFIX_APS, 0, state->frame->first_nal);
    state->frame->first_nal = false;
#ifdef KVZ_DEBUG
  printf("=========== Adaptation Parameter Set  ===========\n");
#endif
   

   WRITE_U(stream, (int)1/*LMCS_APS*/, 3, "aps_params_type");
   WRITE_U(stream, 0, 5, "adaptation_parameter_set_id");
   WRITE_U(stream, state->encoder_control->chroma_format != KVZ_CSP_400, 1, "aps_chroma_present_flag");

   code_lmcs_aps(state, state->slice->lmcs_aps);

   WRITE_U(stream, 0, 1, "aps_extension_flag"); //Implementation when this flag is equal to 1 should be added when it is needed. Currently in the spec we don't have case when this flag is equal to 1
   kvz_bitstream_add_rbsp_trailing_bits(stream);
    
  }
}
