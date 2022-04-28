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

#include "global.h"

#if COMPILE_INTEL_SSE41
#include "uvg266.h"
#if UVG_BIT_DEPTH == 8
#include "strategies/sse41/alf-sse41.h"

#include <immintrin.h>
#include <stdlib.h>

#include "strategyselector.h"

static void alf_derive_classification_blk_sse41(encoder_state_t * const state,
  const int shift,
  const int n_height,
  const int n_width,
  const int blk_pos_x,
  const int blk_pos_y,
  const int blk_dst_x,
  const int blk_dst_y,
  const int vb_ctu_height,
  int vb_pos)
{
  videoframe_t* const frame = state->tile->frame;
  const size_t imgStride = frame->rec->stride;
  const uvg_pixel *  srcExt    = state->tile->frame->rec->y;

  const int imgHExtended = n_height + 4;
  const int imgWExtended = n_width + 4;

  const int posX = blk_pos_x;
  const int posY = blk_pos_y;

  //alf_classifier** classifier = state->tile->frame->alf_info->classifier;

  // 18x40 array
  uint16_t colSums[(CLASSIFICATION_BLK_SIZE + 4) >> 1]
                  [CLASSIFICATION_BLK_SIZE + 8];

  for (int i = 0; i < imgHExtended; i += 2)
  {
    const size_t offset = (i + posY - 3) * imgStride + posX - 3;

    const uvg_pixel*imgY0 = &srcExt[offset];
    const uvg_pixel*imgY1 = &srcExt[offset + imgStride];
    const uvg_pixel*imgY2 = &srcExt[offset + imgStride * 2];
    const uvg_pixel*imgY3 = &srcExt[offset + imgStride * 3];

    // pixel padding for gradient calculation
    int pos      = blk_dst_y - 2 + i;
    int posInCTU = pos & (vb_ctu_height - 1);
    if (pos > 0 && posInCTU == vb_pos - 2)
    {
      imgY3 = imgY2;
    }
    else if (pos > 0 && posInCTU == vb_pos)
    {
      imgY0 = imgY1;
    }

    __m128i prev = _mm_setzero_si128();

    for (int j = 0; j < imgWExtended; j += 8)
    {
      const __m128i x00 = _mm_loadu_si128((const __m128i*) (imgY0 + j));
      const __m128i x01 = _mm_loadu_si128((const __m128i*) (imgY1 + j));
      const __m128i x02 = _mm_loadu_si128((const __m128i*) (imgY2 + j));
      const __m128i x03 = _mm_loadu_si128((const __m128i*) (imgY3 + j));

      const __m128i x04 = _mm_loadu_si128((const __m128i*) (imgY0 + j + 2));
      const __m128i x05 = _mm_loadu_si128((const __m128i*) (imgY1 + j + 2));
      const __m128i x06 = _mm_loadu_si128((const __m128i*) (imgY2 + j + 2));
      const __m128i x07 = _mm_loadu_si128((const __m128i*) (imgY3 + j + 2));

      const __m128i x0 = _mm_unpacklo_epi8(x00, _mm_setzero_si128());
      const __m128i x1 = _mm_unpacklo_epi8(x01, _mm_setzero_si128());
      const __m128i x2 = _mm_unpacklo_epi8(x02, _mm_setzero_si128());
      const __m128i x3 = _mm_unpacklo_epi8(x03, _mm_setzero_si128());

      const __m128i x4 = _mm_unpacklo_epi8(x04, _mm_setzero_si128());
      const __m128i x5 = _mm_unpacklo_epi8(x05, _mm_setzero_si128());
      const __m128i x6 = _mm_unpacklo_epi8(x06, _mm_setzero_si128());
      const __m128i x7 = _mm_unpacklo_epi8(x07, _mm_setzero_si128());

      const __m128i nw = _mm_blend_epi16(x0, x1, 0xaa);
      const __m128i n = _mm_blend_epi16(x0, x5, 0x55);
      const __m128i ne = _mm_blend_epi16(x4, x5, 0xaa);
      const __m128i w = _mm_blend_epi16(x1, x2, 0xaa);
      const __m128i e = _mm_blend_epi16(x5, x6, 0xaa);
      const __m128i sw = _mm_blend_epi16(x2, x3, 0xaa);
      const __m128i s = _mm_blend_epi16(x2, x7, 0x55);
      const __m128i se = _mm_blend_epi16(x6, x7, 0xaa);

      __m128i c = _mm_blend_epi16(x1, x6, 0x55);
      c = _mm_add_epi16(c, c);
      __m128i d = _mm_shuffle_epi8(c, _mm_setr_epi8(2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13));

      const __m128i ver = _mm_abs_epi16(_mm_sub_epi16(c, _mm_add_epi16(n, s)));
      const __m128i hor = _mm_abs_epi16(_mm_sub_epi16(d, _mm_add_epi16(w, e)));
      const __m128i di0 = _mm_abs_epi16(_mm_sub_epi16(d, _mm_add_epi16(nw, se)));
      const __m128i di1 = _mm_abs_epi16(_mm_sub_epi16(d, _mm_add_epi16(ne, sw)));

      const __m128i hv = _mm_hadd_epi16(ver, hor);
      const __m128i di = _mm_hadd_epi16(di0, di1);
      const __m128i all = _mm_hadd_epi16(hv, di);

      const __m128i t = _mm_blend_epi16(all, prev, 0xaa);
      _mm_storeu_si128((__m128i*) & colSums[i >> 1][j], _mm_hadd_epi16(t, all));
      prev = all; 

      if (j + 8 < imgWExtended)
      {
        j += 8;

        const __m128i x0 = _mm_unpackhi_epi8(x00, _mm_setzero_si128());
        const __m128i x1 = _mm_unpackhi_epi8(x01, _mm_setzero_si128());
        const __m128i x2 = _mm_unpackhi_epi8(x02, _mm_setzero_si128());
        const __m128i x3 = _mm_unpackhi_epi8(x03, _mm_setzero_si128());

        const __m128i x4 = _mm_unpackhi_epi8(x04, _mm_setzero_si128());
        const __m128i x5 = _mm_unpackhi_epi8(x05, _mm_setzero_si128());
        const __m128i x6 = _mm_unpackhi_epi8(x06, _mm_setzero_si128());
        const __m128i x7 = _mm_unpackhi_epi8(x07, _mm_setzero_si128());

        const __m128i nw = _mm_blend_epi16(x0, x1, 0xaa);
        const __m128i n = _mm_blend_epi16(x0, x5, 0x55);
        const __m128i ne = _mm_blend_epi16(x4, x5, 0xaa);
        const __m128i w = _mm_blend_epi16(x1, x2, 0xaa);
        const __m128i e = _mm_blend_epi16(x5, x6, 0xaa);
        const __m128i sw = _mm_blend_epi16(x2, x3, 0xaa);
        const __m128i s = _mm_blend_epi16(x2, x7, 0x55);
        const __m128i se = _mm_blend_epi16(x6, x7, 0xaa);

        __m128i c = _mm_blend_epi16(x1, x6, 0x55);
        c = _mm_add_epi16(c, c);
        __m128i d = _mm_shuffle_epi8(c, _mm_setr_epi8(2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13));

        const __m128i ver = _mm_abs_epi16(_mm_sub_epi16(c, _mm_add_epi16(n, s)));
        const __m128i hor = _mm_abs_epi16(_mm_sub_epi16(d, _mm_add_epi16(w, e)));
        const __m128i di0 = _mm_abs_epi16(_mm_sub_epi16(d, _mm_add_epi16(nw, se)));
        const __m128i di1 = _mm_abs_epi16(_mm_sub_epi16(d, _mm_add_epi16(ne, sw)));

        const __m128i hv = _mm_hadd_epi16(ver, hor);
        const __m128i di = _mm_hadd_epi16(di0, di1);
        const __m128i all = _mm_hadd_epi16(hv, di);

        const __m128i t = _mm_blend_epi16(all, prev, 0xaa);
        _mm_storeu_si128((__m128i*) & colSums[i >> 1][j], _mm_hadd_epi16(t, all));
        prev = all;
      }
    }
  }

  for (int i = 0; i < (n_height >> 1); i += 4)
  {
    for (int j = 0; j < n_width; j += 8)
    {
      __m128i x0, x1, x2, x3, x4, x5, x6, x7;

      const uint32_t z = (2 * i + blk_pos_y) & (vb_ctu_height - 1);
      const uint32_t z2 = (2 * i + 4 + blk_pos_y) & (vb_ctu_height - 1);

      x0 = (z == vb_pos) ? _mm_setzero_si128() : _mm_loadu_si128((__m128i *) &colSums[i + 0][j + 4]);
      x1 = _mm_loadu_si128((__m128i *) &colSums[i + 1][j + 4]);
      x2 = _mm_loadu_si128((__m128i *) &colSums[i + 2][j + 4]);
      x3 = (z == vb_pos - 4) ? _mm_setzero_si128() : _mm_loadu_si128((__m128i *) &colSums[i + 3][j + 4]);

      x4 = (z2 == vb_pos) ? _mm_setzero_si128() : _mm_loadu_si128((__m128i *) &colSums[i + 2][j + 4]);
      x5 = _mm_loadu_si128((__m128i *) &colSums[i + 3][j + 4]);
      x6 = _mm_loadu_si128((__m128i *) &colSums[i + 4][j + 4]);
      x7 = (z2 == vb_pos - 4) ? _mm_setzero_si128() : _mm_loadu_si128((__m128i *) &colSums[i + 5][j + 4]);

      __m128i x0l = _mm_cvtepu16_epi32(x0);
      __m128i x0h = _mm_unpackhi_epi16(x0, _mm_setzero_si128());
      __m128i x1l = _mm_cvtepu16_epi32(x1);
      __m128i x1h = _mm_unpackhi_epi16(x1, _mm_setzero_si128());
      __m128i x2l = _mm_cvtepu16_epi32(x2);
      __m128i x2h = _mm_unpackhi_epi16(x2, _mm_setzero_si128());
      __m128i x3l = _mm_cvtepu16_epi32(x3);
      __m128i x3h = _mm_unpackhi_epi16(x3, _mm_setzero_si128());
      __m128i x4l = _mm_cvtepu16_epi32(x4);
      __m128i x4h = _mm_unpackhi_epi16(x4, _mm_setzero_si128());
      __m128i x5l = _mm_cvtepu16_epi32(x5);
      __m128i x5h = _mm_unpackhi_epi16(x5, _mm_setzero_si128());
      __m128i x6l = _mm_cvtepu16_epi32(x6);
      __m128i x6h = _mm_unpackhi_epi16(x6, _mm_setzero_si128());
      __m128i x7l = _mm_cvtepu16_epi32(x7);
      __m128i x7h = _mm_unpackhi_epi16(x7, _mm_setzero_si128());

      x0l = _mm_add_epi32(x0l, x1l);
      x2l = _mm_add_epi32(x2l, x3l);
      x4l = _mm_add_epi32(x4l, x5l);
      x6l = _mm_add_epi32(x6l, x7l);
      x0h = _mm_add_epi32(x0h, x1h);
      x2h = _mm_add_epi32(x2h, x3h);
      x4h = _mm_add_epi32(x4h, x5h);
      x6h = _mm_add_epi32(x6h, x7h);

      x0l = _mm_add_epi32(x0l, x2l);
      x4l = _mm_add_epi32(x4l, x6l);
      x0h = _mm_add_epi32(x0h, x2h);
      x4h = _mm_add_epi32(x4h, x6h);

      x2l = _mm_unpacklo_epi32(x0l, x4l);
      x2h = _mm_unpackhi_epi32(x0l, x4l);
      x6l = _mm_unpacklo_epi32(x0h, x4h);
      x6h = _mm_unpackhi_epi32(x0h, x4h);

      __m128i sumV  = _mm_unpacklo_epi32(x2l, x6l);
      __m128i sumH  = _mm_unpackhi_epi32(x2l, x6l);
      __m128i sumD0 = _mm_unpacklo_epi32(x2h, x6h);
      __m128i sumD1 = _mm_unpackhi_epi32(x2h, x6h);

      //      uint32_t tempAct = sumV + sumH;
      __m128i tempAct = _mm_add_epi32(sumV, sumH);

      //      const uint32_t activity = std::min<uint32_t>(15, tempAct * scale >> shift);
      //      static const uint8_t th[16] = { 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4 };
      //      uint8_t class_idx = th[activity];
      const uint32_t scale  = (z == vb_pos - 4 || z == vb_pos) ? 96 : 64;
      const uint32_t scale2 = (z2 == vb_pos - 4 || z2 == vb_pos) ? 96 : 64;
      __m128i activity = _mm_mullo_epi32(tempAct, _mm_unpacklo_epi64(_mm_set1_epi32(scale), _mm_set1_epi32(scale2)));
      activity         = _mm_srl_epi32(activity, _mm_cvtsi32_si128(shift));
      activity         = _mm_min_epi32(activity, _mm_set1_epi32(15));
      __m128i class_idx = _mm_shuffle_epi8(_mm_setr_epi8(0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4), activity);

      //      if (sumV > sumH)
      //      {
      //        hv1       = sumV;
      //        hv0       = sumH;
      //        dirTempHV = 0;
      //      }
      //      else
      //      {
      //        hv1       = sumH;
      //        hv0       = sumV;
      //        dirTempHV = 1;
      //      }
      __m128i dirTempHVMinus1 = _mm_cmpgt_epi32(sumV, sumH);
      __m128i hv1             = _mm_max_epi32(sumV, sumH);
      __m128i hv0             = _mm_min_epi32(sumV, sumH);

      //      if (sumD0 > sumD1)
      //      {
      //        d1       = sumD0;
      //        d0       = sumD1;
      //        dirTempD = 0;
      //      }
      //      else
      //      {
      //        d1       = sumD1;
      //        d0       = sumD0;
      //        dirTempD = 1;
      //      }
      __m128i dirTempDMinus1 = _mm_cmpgt_epi32(sumD0, sumD1);
      __m128i d1             = _mm_max_epi32(sumD0, sumD1);
      __m128i d0             = _mm_min_epi32(sumD0, sumD1);

      //      int dirIdx;
      //      if (d1 * hv0 > hv1 * d0)
      //      {
      //        hvd1   = d1;
      //        hvd0   = d0;
      //        dirIdx = 0;
      //      }
      //      else
      //      {
      //        hvd1   = hv1;
      //        hvd0   = hv0;
      //        dirIdx = 2;
      //      }
      __m128i a      = _mm_xor_si128(_mm_mullo_epi32(d1, hv0), _mm_set1_epi32(0x80000000));
      __m128i b      = _mm_xor_si128(_mm_mullo_epi32(hv1, d0), _mm_set1_epi32(0x80000000));
      __m128i dirIdx = _mm_cmpgt_epi32(a, b);
      __m128i hvd1   = _mm_blendv_epi8(hv1, d1, dirIdx);
      __m128i hvd0   = _mm_blendv_epi8(hv0, d0, dirIdx);

      //      if (hvd1 * 2 > 9 * hvd0)
      //      {
      //        class_idx += (dirIdx + 2) * 5;
      //      }
      //      else if (hvd1 > 2 * hvd0)
      //      {
      //        class_idx += (dirIdx + 1) * 5;
      //      }
      __m128i strength1 = _mm_cmpgt_epi32(hvd1, _mm_add_epi32(hvd0, hvd0));
      __m128i strength2 = _mm_cmpgt_epi32(_mm_add_epi32(hvd1, hvd1), _mm_add_epi32(hvd0, _mm_slli_epi32(hvd0, 3)));
      __m128i offset    = _mm_and_si128(strength1, _mm_set1_epi32(5));
      class_idx          = _mm_add_epi32(class_idx, offset);
      class_idx          = _mm_add_epi32(class_idx, _mm_and_si128(strength2, _mm_set1_epi32(5)));
      offset            = _mm_andnot_si128(dirIdx, offset);
      offset            = _mm_add_epi32(offset, offset);
      class_idx          = _mm_add_epi32(class_idx, offset);

      //      uint8_t transpose_idx = 2 * dirTempD + dirTempHV;
      __m128i transpose_idx = _mm_set1_epi32(3);
      transpose_idx         = _mm_add_epi32(transpose_idx, dirTempHVMinus1);
      transpose_idx         = _mm_add_epi32(transpose_idx, dirTempDMinus1);
      transpose_idx         = _mm_add_epi32(transpose_idx, dirTempDMinus1);

      int yOffset = 2 * i + blk_pos_y;
      int xOffset = j + blk_pos_x;

      static_assert(sizeof(alf_classifier) == 2, "alf_classifier type must be 16 bits wide");
      __m128i v;
      v = _mm_unpacklo_epi8(class_idx, transpose_idx);
      v = _mm_shuffle_epi8(v, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 8, 9, 8, 9, 8, 9, 8, 9));
      _mm_storeu_si128((__m128i *) (state->tile->frame->alf_info->classifier[yOffset] + xOffset), v);
      _mm_storeu_si128((__m128i *) (state->tile->frame->alf_info->classifier[yOffset + 1] + xOffset), v);
      _mm_storeu_si128((__m128i *) (state->tile->frame->alf_info->classifier[yOffset + 2] + xOffset), v);
      _mm_storeu_si128((__m128i *) (state->tile->frame->alf_info->classifier[yOffset + 3] + xOffset), v);
      v = _mm_unpackhi_epi8(class_idx, transpose_idx);
      v = _mm_shuffle_epi8(v, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 8, 9, 8, 9, 8, 9, 8, 9));
      _mm_storeu_si128((__m128i *) (state->tile->frame->alf_info->classifier[yOffset + 4] + xOffset), v);
      _mm_storeu_si128((__m128i *) (state->tile->frame->alf_info->classifier[yOffset + 5] + xOffset), v);
      _mm_storeu_si128((__m128i *) (state->tile->frame->alf_info->classifier[yOffset + 6] + xOffset), v);
      _mm_storeu_si128((__m128i *) (state->tile->frame->alf_info->classifier[yOffset + 7] + xOffset), v);
    }
  }
}


INLINE static void process2coeffs_5x5(__m128i params[2][3], __m128i *cur, __m128i *accumA, __m128i *accumB, const int i, const uvg_pixel* ptr0, const uvg_pixel* ptr1, const uvg_pixel* ptr2, const uvg_pixel* ptr3) {
  const __m128i val00 = _mm_sub_epi16(_mm_unpacklo_epi8(_mm_loadu_si128((const __m128i*) ptr0), _mm_setzero_si128()), *cur);
  const __m128i val10 = _mm_sub_epi16(_mm_unpacklo_epi8(_mm_loadu_si128((const __m128i*) ptr2), _mm_setzero_si128()), *cur);
  const __m128i val01 = _mm_sub_epi16(_mm_unpacklo_epi8(_mm_loadu_si128((const __m128i*) ptr1), _mm_setzero_si128()), *cur);
  const __m128i val11 = _mm_sub_epi16(_mm_unpacklo_epi8(_mm_loadu_si128((const __m128i*) ptr3), _mm_setzero_si128()), *cur);
  __m128i val01A = _mm_unpacklo_epi16(val00, val10);
  __m128i val01B = _mm_unpackhi_epi16(val00, val10);
  __m128i val01C = _mm_unpacklo_epi16(val01, val11);
  __m128i val01D = _mm_unpackhi_epi16(val01, val11);

  __m128i limit01A = params[1][i];

  val01A = _mm_min_epi16(val01A, limit01A);
  val01B = _mm_min_epi16(val01B, limit01A);
  val01C = _mm_min_epi16(val01C, limit01A);
  val01D = _mm_min_epi16(val01D, limit01A);

  limit01A = _mm_sub_epi16(_mm_setzero_si128(), limit01A);

  val01A = _mm_max_epi16(val01A, limit01A);
  val01B = _mm_max_epi16(val01B, limit01A);
  val01C = _mm_max_epi16(val01C, limit01A);
  val01D = _mm_max_epi16(val01D, limit01A);

  val01A = _mm_add_epi16(val01A, val01C);
  val01B = _mm_add_epi16(val01B, val01D);

  __m128i coeff01A = params[0][i];

  *accumA = _mm_add_epi32(*accumA, _mm_madd_epi16(val01A, coeff01A));
  *accumB = _mm_add_epi32(*accumB, _mm_madd_epi16(val01B, coeff01A));
};


static void alf_filter_5x5_block_sse41(encoder_state_t* const state,
  const uvg_pixel* src_pixels,
  uvg_pixel* dst_pixels,
  const int src_stride,
  const int dst_stride,
  const short* filter_set,
  const int16_t* fClipSet,
  clp_rng clp_rng,
  const int width,
  const int height,
  int x_pos,
  int y_pos,
  int blk_dst_x,
  int blk_dst_y,
  int vb_pos,
  const int vb_ctu_height)
{


  assert((vb_ctu_height & (vb_ctu_height - 1)) == 0 && "vb_ctu_height must be a power of 2");
  
  //alf_component_id compId = COMPONENT_Cb;

  const size_t srcStride = src_stride;
  const size_t dstStride = dst_stride;

  const int SHIFT = state->encoder_control->bitdepth - 1;
  const int ROUND = 1 << (SHIFT - 1);
  const __m128i mmOffset1 = _mm_set1_epi32((1 << ((SHIFT + 3) - 1)) - ROUND);

  const size_t STEP_X = 8;
  const size_t STEP_Y = 4;

  assert(y_pos % STEP_Y == 0 && "Wrong startHeight in filtering");
  assert(x_pos % STEP_X == 0 && "Wrong startWidth in filtering");
  assert(height % STEP_Y == 0 && "Wrong endHeight in filtering");
  assert(width % 4 == 0 && "Wrong endWidth in filtering");

  const uvg_pixel* src = src_pixels + y_pos * srcStride + x_pos;
  uvg_pixel* dst = dst_pixels + blk_dst_y * dstStride + blk_dst_x;



  const __m128i mmOffset = _mm_set1_epi32(ROUND);
  const __m128i mmMin = _mm_set1_epi16(clp_rng.min);
  const __m128i mmMax = _mm_set1_epi16(clp_rng.max);

  __m128i params[2][3];
  __m128i fs = _mm_loadu_si128((__m128i*) filter_set);
  params[0][0] = _mm_shuffle_epi32(fs, 0x00);
  params[0][1] = _mm_shuffle_epi32(fs, 0x55);
  params[0][2] = _mm_shuffle_epi32(fs, 0xaa);
  __m128i fc = _mm_loadu_si128((__m128i*) fClipSet);
  params[1][0] = _mm_shuffle_epi32(fc, 0x00);
  params[1][1] = _mm_shuffle_epi32(fc, 0x55);
  params[1][2] = _mm_shuffle_epi32(fc, 0xaa);

  const __m128i mask = _mm_set_epi8(16, 16, 16, 16, 16, 16, 16, 16, 14, 12, 10, 8, 6, 4, 2, 0);

  for (size_t i = 0; i < height; i += STEP_Y)
  {
    for (size_t j = 0; j < width; j += STEP_X)
    {

      for (size_t ii = 0; ii < STEP_Y; ii++)
      {
        const uvg_pixel* pImg0, * pImg1, * pImg2, * pImg3, * pImg4;

        pImg0 = src + j + ii * srcStride;
        pImg1 = pImg0 + srcStride;
        pImg2 = pImg0 - srcStride;
        pImg3 = pImg1 + srcStride;
        pImg4 = pImg2 - srcStride;

        const int yVb = (blk_dst_y + i + ii) & (vb_ctu_height - 1);
        if (yVb < vb_pos && (yVb >= vb_pos - 2))   // above
        {
          pImg1 = (yVb == vb_pos - 1) ? pImg0 : pImg1;
          pImg3 = (yVb >= vb_pos - 2) ? pImg1 : pImg3;

          pImg2 = (yVb == vb_pos - 1) ? pImg0 : pImg2;
          pImg4 = (yVb >= vb_pos - 2) ? pImg2 : pImg4;
        }
        else if (yVb >= vb_pos && (yVb <= vb_pos + 1))   // bottom
        {
          pImg2 = (yVb == vb_pos) ? pImg0 : pImg2;
          pImg4 = (yVb <= vb_pos + 1) ? pImg2 : pImg4;

          pImg1 = (yVb == vb_pos) ? pImg0 : pImg1;
          pImg3 = (yVb <= vb_pos + 1) ? pImg1 : pImg3;
        }
        __m128i cur = _mm_unpacklo_epi8(_mm_loadu_si128((const __m128i*) pImg0), _mm_setzero_si128());
        
        __m128i accumA = mmOffset;
        __m128i accumB = mmOffset;

        

        process2coeffs_5x5(params, &cur, &accumA, &accumB, 0, pImg3 + 0, pImg4 + 0, pImg1 + 1, pImg2 - 1);
        process2coeffs_5x5(params, &cur, &accumA, &accumB, 1, pImg1 + 0, pImg2 + 0, pImg1 - 1, pImg2 + 1);
        process2coeffs_5x5(params, &cur, &accumA, &accumB, 2, pImg0 + 2, pImg0 - 2, pImg0 + 1, pImg0 - 1);
        bool isNearVBabove = yVb < vb_pos && (yVb >= vb_pos - 1);
        bool isNearVBbelow = yVb >= vb_pos && (yVb <= vb_pos);
        if (!(isNearVBabove || isNearVBbelow))
        {
          accumA = _mm_srai_epi32(accumA, SHIFT);
          accumB = _mm_srai_epi32(accumB, SHIFT);
        }
        else
        {
          accumA = _mm_srai_epi32(_mm_add_epi32(accumA, mmOffset1), SHIFT + 3);
          accumB = _mm_srai_epi32(_mm_add_epi32(accumB, mmOffset1), SHIFT + 3);
        }
        accumA = _mm_packs_epi32(accumA, accumB);
        accumA = _mm_add_epi16(accumA, cur);
        accumA = _mm_min_epi16(mmMax, _mm_max_epi16(accumA, mmMin));
        
        if (j + STEP_X <= width)
        {
          //_mm_storeu_si128((__m128i*) (dst + ii * dstStride + j), accumA);
          _mm_storel_epi64((__m128i*) (dst + ii * dstStride + j), _mm_shuffle_epi8(accumA, mask));
        }
        else
        {
          //_mm_storel_epi64((__m128i*) (dst + ii * dstStride + j), accumA);
          _mm_store_ss((float*) (dst + ii * dstStride + j), _mm_castsi128_ps(_mm_shuffle_epi8(accumA, mask)));
        }
      }

    }

    src += srcStride * STEP_Y;
    dst += dstStride * STEP_Y;
  }
}

#define sh(x) 0x0202 * (x & 7) + 0x0100 + 0x1010 * (x & 8)

static const uint16_t shuffleTab[4][2][8] = {
  {
    { sh(0), sh(1), sh(2), sh(3), sh(4), sh(5), sh(6), sh(7) },
    { sh(8), sh(9), sh(10), sh(11), sh(12), sh(13), sh(14), sh(15) },
  },
  {
    { sh(9), sh(4), sh(10), sh(8), sh(1), sh(5), sh(11), sh(7) },
    { sh(3), sh(0), sh(2), sh(6), sh(12), sh(13), sh(14), sh(15) },
  },
  {
    { sh(0), sh(3), sh(2), sh(1), sh(8), sh(7), sh(6), sh(5) },
    { sh(4), sh(9), sh(10), sh(11), sh(12), sh(13), sh(14), sh(15) },
  },
  {
    { sh(9), sh(8), sh(10), sh(4), sh(3), sh(7), sh(11), sh(5) },
    { sh(1), sh(0), sh(2), sh(6), sh(12), sh(13), sh(14), sh(15) },
  },
};



INLINE static void process2coeffs_7x7(__m128i params[2][2][6], __m128i *cur, __m128i *accumA, __m128i *accumB, const int i, const uvg_pixel* ptr0, const uvg_pixel* ptr1, const uvg_pixel* ptr2, const uvg_pixel* ptr3) {
  const __m128i val00 = _mm_sub_epi16(_mm_unpacklo_epi8(_mm_loadu_si128((const __m128i*) ptr0), _mm_setzero_si128()), *cur);
  const __m128i val10 = _mm_sub_epi16(_mm_unpacklo_epi8(_mm_loadu_si128((const __m128i*) ptr2), _mm_setzero_si128()), *cur);
  const __m128i val01 = _mm_sub_epi16(_mm_unpacklo_epi8(_mm_loadu_si128((const __m128i*) ptr1), _mm_setzero_si128()), *cur);
  const __m128i val11 = _mm_sub_epi16(_mm_unpacklo_epi8(_mm_loadu_si128((const __m128i*) ptr3), _mm_setzero_si128()), *cur);

  __m128i val01A = _mm_unpacklo_epi16(val00, val10);
  __m128i val01B = _mm_unpackhi_epi16(val00, val10);
  __m128i val01C = _mm_unpacklo_epi16(val01, val11);
  __m128i val01D = _mm_unpackhi_epi16(val01, val11);

  __m128i limit01A = params[0][1][i];
  __m128i limit01B = params[1][1][i];

  val01A = _mm_min_epi16(val01A, limit01A);
  val01B = _mm_min_epi16(val01B, limit01B);
  val01C = _mm_min_epi16(val01C, limit01A);
  val01D = _mm_min_epi16(val01D, limit01B);

  limit01A = _mm_sub_epi16(_mm_setzero_si128(), limit01A);
  limit01B = _mm_sub_epi16(_mm_setzero_si128(), limit01B);

  val01A = _mm_max_epi16(val01A, limit01A);
  val01B = _mm_max_epi16(val01B, limit01B);
  val01C = _mm_max_epi16(val01C, limit01A);
  val01D = _mm_max_epi16(val01D, limit01B);

  val01A = _mm_add_epi16(val01A, val01C);
  val01B = _mm_add_epi16(val01B, val01D);

  const __m128i coeff01A = params[0][0][i];
  const __m128i coeff01B = params[1][0][i];

  *accumA = _mm_add_epi32(*accumA, _mm_madd_epi16(val01A, coeff01A));
  *accumB = _mm_add_epi32(*accumB, _mm_madd_epi16(val01B, coeff01B));
};




static void alf_filter_7x7_block_sse41(encoder_state_t* const state,
  const uvg_pixel* src_pixels,
  uvg_pixel* dst_pixels,
  const int src_stride,
  const int dst_stride,
  const short* filter_set,
  const int16_t* fClipSet,
  clp_rng clp_rng,
  const int width,
  const int height,
  int x_pos,
  int y_pos,
  int blk_dst_x,
  int blk_dst_y,
  int vb_pos,
  const int vb_ctu_height)
{
  assert((vb_ctu_height & (vb_ctu_height - 1)) == 0 && "vb_ctu_height must be a power of 2");
  //alf_component_id compId = COMPONENT_Y;


  const size_t srcStride = src_stride;
  const size_t dstStride = dst_stride;

  const int SHIFT = state->encoder_control->bitdepth - 1;
  const int ROUND = 1 << (SHIFT - 1);

  const size_t STEP_X = 8;
  const size_t STEP_Y = 4;

  assert(y_pos % STEP_Y == 0 && "Wrong startHeight in filtering");
  assert(x_pos % STEP_X == 0 && "Wrong startWidth in filtering");
  assert(height % STEP_Y == 0 && "Wrong endHeight in filtering");
  assert(width % STEP_X == 0 && "Wrong endWidth in filtering");

  const uvg_pixel* src = src_pixels + y_pos * srcStride + x_pos;
  uvg_pixel* dst = dst_pixels + blk_dst_y * dstStride + blk_dst_x;

  const __m128i mmOffset = _mm_set1_epi32(ROUND);
  const __m128i mmOffset1 = _mm_set1_epi32((1 << ((SHIFT + 3) - 1)) - ROUND);
  const __m128i mmMin = _mm_set1_epi16(clp_rng.min);
  const __m128i mmMax = _mm_set1_epi16(clp_rng.max);

  const __m128i mask = _mm_set_epi8(16, 16, 16, 16, 16, 16, 16, 16, 14, 12, 10, 8, 6, 4, 2, 0);

  for (size_t i = 0; i < height; i += STEP_Y)
  {
    const alf_classifier* pClass = state->tile->frame->alf_info->classifier[blk_dst_y + i] + blk_dst_x;

    for (size_t j = 0; j < width; j += STEP_X)
    {
      __m128i params[2][2][6];

      for (int k = 0; k < 2; ++k)
      {
        const alf_classifier* cl = &pClass[j + 4 * k];

        const int transpose_idx = cl->transpose_idx;
        const int class_idx = cl->class_idx;

        static_assert(sizeof(*filter_set) == 2, "ALF coeffs must be 16-bit wide");
        static_assert(sizeof(*fClipSet) == 2, "ALF clip values must be 16-bit wide");

        __m128i rawCoeff0, rawCoeff1;
        __m128i rawClip0, rawClip1;

        rawCoeff0 = _mm_loadu_si128((const __m128i*) (filter_set + class_idx * MAX_NUM_ALF_LUMA_COEFF));
        rawCoeff1 = _mm_loadl_epi64((const __m128i*) (filter_set + class_idx * MAX_NUM_ALF_LUMA_COEFF + 8));

        rawClip0 = _mm_loadu_si128((const __m128i*) (fClipSet + class_idx * MAX_NUM_ALF_LUMA_COEFF));
        rawClip1 = _mm_loadl_epi64((const __m128i*) (fClipSet + class_idx * MAX_NUM_ALF_LUMA_COEFF + 8));

        const __m128i s0 = _mm_loadu_si128((const __m128i*) shuffleTab[transpose_idx][0]);
        const __m128i s1 = _mm_xor_si128(s0, _mm_set1_epi8((char)0x80));
        const __m128i s2 = _mm_loadu_si128((const __m128i*) shuffleTab[transpose_idx][1]);
        const __m128i s3 = _mm_xor_si128(s2, _mm_set1_epi8((char)0x80));

        const __m128i rawCoeffLo = _mm_or_si128(_mm_shuffle_epi8(rawCoeff0, s0), _mm_shuffle_epi8(rawCoeff1, s1));
        const __m128i rawCoeffHi = _mm_or_si128(_mm_shuffle_epi8(rawCoeff0, s2), _mm_shuffle_epi8(rawCoeff1, s3));
        const __m128i rawClipLo = _mm_or_si128(_mm_shuffle_epi8(rawClip0, s0), _mm_shuffle_epi8(rawClip1, s1));
        const __m128i rawClipHi = _mm_or_si128(_mm_shuffle_epi8(rawClip0, s2), _mm_shuffle_epi8(rawClip1, s3));

        params[k][0][0] = _mm_shuffle_epi32(rawCoeffLo, 0x00);
        params[k][0][1] = _mm_shuffle_epi32(rawCoeffLo, 0x55);
        params[k][0][2] = _mm_shuffle_epi32(rawCoeffLo, 0xaa);
        params[k][0][3] = _mm_shuffle_epi32(rawCoeffLo, 0xff);
        params[k][0][4] = _mm_shuffle_epi32(rawCoeffHi, 0x00);
        params[k][0][5] = _mm_shuffle_epi32(rawCoeffHi, 0x55);
        params[k][1][0] = _mm_shuffle_epi32(rawClipLo, 0x00);
        params[k][1][1] = _mm_shuffle_epi32(rawClipLo, 0x55);
        params[k][1][2] = _mm_shuffle_epi32(rawClipLo, 0xaa);
        params[k][1][3] = _mm_shuffle_epi32(rawClipLo, 0xff);
        params[k][1][4] = _mm_shuffle_epi32(rawClipHi, 0x00);
        params[k][1][5] = _mm_shuffle_epi32(rawClipHi, 0x55);
      }

      for (size_t ii = 0; ii < STEP_Y; ii++)
      {
        const uvg_pixel* pImg0, * pImg1, * pImg2, * pImg3, * pImg4, * pImg5, * pImg6;

        pImg0 = src + j + ii * srcStride;
        pImg1 = pImg0 + srcStride;
        pImg2 = pImg0 - srcStride;
        pImg3 = pImg1 + srcStride;
        pImg4 = pImg2 - srcStride;
        pImg5 = pImg3 + srcStride;
        pImg6 = pImg4 - srcStride;

        const int yVb = (blk_dst_y + i + ii) & (vb_ctu_height - 1);
        if (yVb < vb_pos && (yVb >= vb_pos - 4))   // above
        {
          pImg1 = (yVb == vb_pos - 1) ? pImg0 : pImg1;
          pImg3 = (yVb >= vb_pos - 2) ? pImg1 : pImg3;
          pImg5 = (yVb >= vb_pos - 3) ? pImg3 : pImg5;

          pImg2 = (yVb == vb_pos - 1) ? pImg0 : pImg2;
          pImg4 = (yVb >= vb_pos - 2) ? pImg2 : pImg4;
          pImg6 = (yVb >= vb_pos - 3) ? pImg4 : pImg6;
        }
        else if (yVb >= vb_pos && (yVb <= vb_pos + 3))   // bottom
        {
          pImg2 = (yVb == vb_pos) ? pImg0 : pImg2;
          pImg4 = (yVb <= vb_pos + 1) ? pImg2 : pImg4;
          pImg6 = (yVb <= vb_pos + 2) ? pImg4 : pImg6;

          pImg1 = (yVb == vb_pos) ? pImg0 : pImg1;
          pImg3 = (yVb <= vb_pos + 1) ? pImg1 : pImg3;
          pImg5 = (yVb <= vb_pos + 2) ? pImg3 : pImg5;
        }
        __m128i cur = _mm_unpacklo_epi8(_mm_loadu_si128((const __m128i*) pImg0), _mm_setzero_si128());

        __m128i accumA = mmOffset;
        __m128i accumB = mmOffset;

        process2coeffs_7x7(params, &cur, &accumA, &accumB, 0, pImg5 + 0, pImg6 + 0, pImg3 + 1, pImg4 - 1);
        process2coeffs_7x7(params, &cur, &accumA, &accumB, 1, pImg3 + 0, pImg4 + 0, pImg3 - 1, pImg4 + 1);
        process2coeffs_7x7(params, &cur, &accumA, &accumB, 2, pImg1 + 2, pImg2 - 2, pImg1 + 1, pImg2 - 1);
        process2coeffs_7x7(params, &cur, &accumA, &accumB, 3, pImg1 + 0, pImg2 + 0, pImg1 - 1, pImg2 + 1);
        process2coeffs_7x7(params, &cur, &accumA, &accumB, 4, pImg1 - 2, pImg2 + 2, pImg0 + 3, pImg0 - 3);
        process2coeffs_7x7(params, &cur, &accumA, &accumB, 5, pImg0 + 2, pImg0 - 2, pImg0 + 1, pImg0 - 1);


        bool isNearVBabove = yVb < vb_pos && (yVb >= vb_pos - 1);
        bool isNearVBbelow = yVb >= vb_pos && (yVb <= vb_pos);
        if (!(isNearVBabove || isNearVBbelow))
        {
          accumA = _mm_srai_epi32(accumA, SHIFT);
          accumB = _mm_srai_epi32(accumB, SHIFT);
        }
        else
        {
          accumA = _mm_srai_epi32(_mm_add_epi32(accumA, mmOffset1), SHIFT + 3);
          accumB = _mm_srai_epi32(_mm_add_epi32(accumB, mmOffset1), SHIFT + 3);
        }
        accumA = _mm_packs_epi32(accumA, accumB);
        accumA = _mm_add_epi16(accumA, cur);
        accumA = _mm_min_epi16(mmMax, _mm_max_epi16(accumA, mmMin));       

        //_mm_storeu_si128((__m128i*) (dst + ii * dstStride + j), accumA);
        _mm_storel_epi64((__m128i*) (dst + ii * dstStride + j), _mm_shuffle_epi8(accumA, mask));
      }
    }

    src += srcStride * STEP_Y;
    dst += dstStride * STEP_Y;
  }
}



#endif // UVG_BIT_DEPTH == 8
#endif //COMPILE_INTEL_SSE41


int uvg_strategy_register_alf_sse41(void* opaque, uint8_t bitdepth) {
  bool success = true;
#if COMPILE_INTEL_SSE41
#if UVG_BIT_DEPTH == 8
  if (bitdepth == 8){
    success &= uvg_strategyselector_register(opaque, "alf_derive_classification_blk", "sse41", 20, &alf_derive_classification_blk_sse41);
    success &= uvg_strategyselector_register(opaque, "alf_filter_5x5_blk", "sse41", 0, &alf_filter_5x5_block_sse41);
    success &= uvg_strategyselector_register(opaque, "alf_filter_7x7_blk", "sse41", 0, &alf_filter_7x7_block_sse41);
  }
#endif // UVG_BIT_DEPTH == 8
#endif
  return success;
}
