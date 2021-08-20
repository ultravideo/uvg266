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

#include "global.h"

#if COMPILE_INTEL_SSE41
#include "kvazaar.h"
#if KVZ_BIT_DEPTH == 8
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
  const kvz_pixel *  srcExt    = state->tile->frame->rec->y;

  const int imgHExtended = n_height + 4;
  const int imgWExtended = n_width + 4;

  const int posX = blk_pos_x;
  const int posY = blk_pos_y;

  alf_classifier** classifier = state->tile->frame->alf_info->classifier;

  // 18x40 array
  uint16_t colSums[(CLASSIFICATION_BLK_SIZE + 4) >> 1]
                  [CLASSIFICATION_BLK_SIZE + 8];

  for (int i = 0; i < imgHExtended; i += 2)
  {
    const size_t offset = (i + posY - 3) * imgStride + posX - 3;

    const kvz_pixel*imgY0 = &srcExt[offset];
    const kvz_pixel*imgY1 = &srcExt[offset + imgStride];
    const kvz_pixel*imgY2 = &srcExt[offset + imgStride * 2];
    const kvz_pixel*imgY3 = &srcExt[offset + imgStride * 3];

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
      //      uint8_t classIdx = th[activity];
      const uint32_t scale  = (z == vb_pos - 4 || z == vb_pos) ? 96 : 64;
      const uint32_t scale2 = (z2 == vb_pos - 4 || z2 == vb_pos) ? 96 : 64;
      __m128i activity = _mm_mullo_epi32(tempAct, _mm_unpacklo_epi64(_mm_set1_epi32(scale), _mm_set1_epi32(scale2)));
      activity         = _mm_srl_epi32(activity, _mm_cvtsi32_si128(shift));
      activity         = _mm_min_epi32(activity, _mm_set1_epi32(15));
      __m128i classIdx = _mm_shuffle_epi8(_mm_setr_epi8(0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4), activity);

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
      //        classIdx += (dirIdx + 2) * 5;
      //      }
      //      else if (hvd1 > 2 * hvd0)
      //      {
      //        classIdx += (dirIdx + 1) * 5;
      //      }
      __m128i strength1 = _mm_cmpgt_epi32(hvd1, _mm_add_epi32(hvd0, hvd0));
      __m128i strength2 = _mm_cmpgt_epi32(_mm_add_epi32(hvd1, hvd1), _mm_add_epi32(hvd0, _mm_slli_epi32(hvd0, 3)));
      __m128i offset    = _mm_and_si128(strength1, _mm_set1_epi32(5));
      classIdx          = _mm_add_epi32(classIdx, offset);
      classIdx          = _mm_add_epi32(classIdx, _mm_and_si128(strength2, _mm_set1_epi32(5)));
      offset            = _mm_andnot_si128(dirIdx, offset);
      offset            = _mm_add_epi32(offset, offset);
      classIdx          = _mm_add_epi32(classIdx, offset);

      //      uint8_t transposeIdx = 2 * dirTempD + dirTempHV;
      __m128i transposeIdx = _mm_set1_epi32(3);
      transposeIdx         = _mm_add_epi32(transposeIdx, dirTempHVMinus1);
      transposeIdx         = _mm_add_epi32(transposeIdx, dirTempDMinus1);
      transposeIdx         = _mm_add_epi32(transposeIdx, dirTempDMinus1);

      int yOffset = 2 * i + blk_pos_y;
      int xOffset = j + blk_pos_x;

      static_assert(sizeof(alf_classifier) == 2, "ALFClassifier type must be 16 bits wide");
      __m128i v;
      v = _mm_unpacklo_epi8(classIdx, transposeIdx);
      v = _mm_shuffle_epi8(v, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 8, 9, 8, 9, 8, 9, 8, 9));
      _mm_storeu_si128((__m128i *) (classifier[yOffset] + xOffset), v);
      _mm_storeu_si128((__m128i *) (classifier[yOffset + 1] + xOffset), v);
      _mm_storeu_si128((__m128i *) (classifier[yOffset + 2] + xOffset), v);
      _mm_storeu_si128((__m128i *) (classifier[yOffset + 3] + xOffset), v);
      v = _mm_unpackhi_epi8(classIdx, transposeIdx);
      v = _mm_shuffle_epi8(v, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 8, 9, 8, 9, 8, 9, 8, 9));
      _mm_storeu_si128((__m128i *) (classifier[yOffset + 4] + xOffset), v);
      _mm_storeu_si128((__m128i *) (classifier[yOffset + 5] + xOffset), v);
      _mm_storeu_si128((__m128i *) (classifier[yOffset + 6] + xOffset), v);
      _mm_storeu_si128((__m128i *) (classifier[yOffset + 7] + xOffset), v);
    }
  }
}

#endif // KVZ_BIT_DEPTH == 8
#endif //COMPILE_INTEL_SSE41


int kvz_strategy_register_alf_sse41(void* opaque, uint8_t bitdepth) {
  bool success = true;
#if COMPILE_INTEL_SSE41
#if KVZ_BIT_DEPTH == 8
  if (bitdepth == 8){
    success &= kvz_strategyselector_register(opaque, "alf_derive_classification_blk", "sse41", 20, &alf_derive_classification_blk_sse41);
  }
#endif // KVZ_BIT_DEPTH == 8
#endif
  return success;
}
