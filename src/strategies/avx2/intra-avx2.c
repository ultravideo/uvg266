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

#include "strategies/avx2/intra-avx2.h"


#if COMPILE_INTEL_AVX2 && defined X86_64
#include "uvg266.h"
#include "cu.h"
#include "tables.h"
#if UVG_BIT_DEPTH == 8

#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "global.h"
#include "intra-avx2.h"
#include "intra_avx2_tables.h"
#include "strategies/avx2/mip_data_avx2.h"
#include "uvg_math.h"

 #include "strategyselector.h"
 #include "strategies/missing-intel-intrinsics.h"


static const int16_t cubic_filter[32][4] =
{
  { 0, 64,  0,  0 },
  { -1, 63,  2,  0 },
  { -2, 62,  4,  0 },
  { -2, 60,  7, -1 },
  { -2, 58, 10, -2 },
  { -3, 57, 12, -2 },
  { -4, 56, 14, -2 },
  { -4, 55, 15, -2 },
  { -4, 54, 16, -2 },
  { -5, 53, 18, -2 },
  { -6, 52, 20, -2 },
  { -6, 49, 24, -3 },
  { -6, 46, 28, -4 },
  { -5, 44, 29, -4 },
  { -4, 42, 30, -4 },
  { -4, 39, 33, -4 },
  { -4, 36, 36, -4 },
  { -4, 33, 39, -4 },
  { -4, 30, 42, -4 },
  { -4, 29, 44, -5 },
  { -4, 28, 46, -6 },
  { -3, 24, 49, -6 },
  { -2, 20, 52, -6 },
  { -2, 18, 53, -5 },
  { -2, 16, 54, -4 },
  { -2, 15, 55, -4 },
  { -2, 14, 56, -4 },
  { -2, 12, 57, -3 },
  { -2, 10, 58, -2 },
  { -1,  7, 60, -2 },
  { 0,  4, 62, -2 },
  { 0,  2, 63, -1 },
};


// Specified in JVET-T2001 8.4.5.2.13 Table 25
// These are the fC interpolation filter coefficients
static const int8_t cubic_filter_8bit_c[32][4] =
{
  { 0, 64,  0,  0 },
  { -1, 63,  2,  0 },
  { -2, 62,  4,  0 },
  { -2, 60,  7, -1 },
  { -2, 58, 10, -2 },
  { -3, 57, 12, -2 },
  { -4, 56, 14, -2 },
  { -4, 55, 15, -2 },
  { -4, 54, 16, -2 },
  { -5, 53, 18, -2 },
  { -6, 52, 20, -2 },
  { -6, 49, 24, -3 },
  { -6, 46, 28, -4 },
  { -5, 44, 29, -4 },
  { -4, 42, 30, -4 },
  { -4, 39, 33, -4 },
  { -4, 36, 36, -4 },
  { -4, 33, 39, -4 },
  { -4, 30, 42, -4 },
  { -4, 29, 44, -5 },
  { -4, 28, 46, -6 },
  { -3, 24, 49, -6 },
  { -2, 20, 52, -6 },
  { -2, 18, 53, -5 },
  { -2, 16, 54, -4 },
  { -2, 15, 55, -4 },
  { -2, 14, 56, -4 },
  { -2, 12, 57, -3 },
  { -2, 10, 58, -2 },
  { -1,  7, 60, -2 },
  { 0,  4, 62, -2 },
  { 0,  2, 63, -1 },
};

// Specified in JVET-T2001 8.4.5.2.13 Table 25
// These are the fG interpolation filter coefficients
static const int8_t cubic_filter_8bit_g[32][4] = 
{
  {16, 32, 16, 0},
  {16, 32, 16, 0},
  {15, 31, 17, 1},
  {15, 31, 17, 1},
  {14, 30, 18, 2},
  {14, 30, 18, 2},
  {13, 29, 19, 3},
  {13, 29, 19, 3},
  {12, 28, 20, 4},
  {12, 28, 20, 4},
  {11, 27, 21, 5},
  {11, 27, 21, 5},
  {10, 26, 22, 6},
  {10, 26, 22, 6},
  { 9, 25, 23, 7},
  { 9, 25, 23, 7},
  { 8, 24, 24, 8},
  { 8, 24, 24, 8},
  { 7, 23, 25, 9},
  { 7, 23, 25, 9},
  { 6, 22, 26, 10},
  { 6, 22, 26, 10},
  { 5, 21, 27, 11},
  { 5, 21, 27, 11},
  { 4, 20, 28, 12},
  { 4, 20, 28, 12},
  { 3, 19, 29, 13},
  { 3, 19, 29, 13},
  { 2, 18, 30, 14},
  { 2, 18, 30, 14},
  { 1, 17, 31, 15},
  { 1, 17, 31, 15}
};


static void angular_pred_w4_ver_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t* delta_int, const int16_t* delta_fract, const int height, const int use_cubic)
{
  const int width = 4;

  const __m256i p_shuf_01 = _mm256_setr_epi8(
    0x00, 0x01, 0x01, 0x02, 0x02, 0x03, 0x03, 0x04,
    0x08, 0x09, 0x09, 0x0a, 0x0a, 0x0b, 0x0b, 0x0c,
    0x00, 0x01, 0x01, 0x02, 0x02, 0x03, 0x03, 0x04,
    0x08, 0x09, 0x09, 0x0a, 0x0a, 0x0b, 0x0b, 0x0c
  );

  const __m256i p_shuf_23 = _mm256_setr_epi8(
    0x02, 0x03, 0x03, 0x04, 0x04, 0x05, 0x05, 0x06,
    0x0a, 0x0b, 0x0b, 0x0c, 0x0c, 0x0d, 0x0d, 0x0e,
    0x02, 0x03, 0x03, 0x04, 0x04, 0x05, 0x05, 0x06,
    0x0a, 0x0b, 0x0b, 0x0c, 0x0c, 0x0d, 0x0d, 0x0e
  );

  const __m256i w_shuf_01 = _mm256_setr_epi8(
    0x00, 0x02, 0x00, 0x02, 0x00, 0x02, 0x00, 0x02,
    0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a,
    0x00, 0x02, 0x00, 0x02, 0x00, 0x02, 0x00, 0x02,
    0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a
  );

  const __m256i w_shuf_23 = _mm256_setr_epi8(
    0x04, 0x06, 0x04, 0x06, 0x04, 0x06, 0x04, 0x06,
    0x0c, 0x0e, 0x0c, 0x0e, 0x0c, 0x0e, 0x0c, 0x0e,
    0x04, 0x06, 0x04, 0x06, 0x04, 0x06, 0x04, 0x06,
    0x0c, 0x0e, 0x0c, 0x0e, 0x0c, 0x0e, 0x0c, 0x0e
  );
  
  int16_t f[4][4] = { { 0 } };

  // For a 4 width block, height must be at least 4. Handle 4 lines at once
  for (int y = 0; y < height; y += 4) {
    if (use_cubic) {
      memcpy(f[0], cubic_filter[delta_fract[y + 0]], 8);
      memcpy(f[1], cubic_filter[delta_fract[y + 1]], 8);
      memcpy(f[2], cubic_filter[delta_fract[y + 2]], 8);
      memcpy(f[3], cubic_filter[delta_fract[y + 3]], 8);
    }
    else {
      for (int yy = 0; yy < 4; ++yy) {
        const int16_t offset = (delta_fract[y + yy] >> 1);
        f[yy][0] = 16 - offset;
        f[yy][1] = 32 - offset;
        f[yy][2] = 16 + offset;
        f[yy][3] = offset;
      }
    }

    // Do 4-tap intra interpolation filtering
    uvg_pixel* p = (uvg_pixel*)ref_main;
    // This solution assumes the delta int values to be 64-bit
    // Cast from 16-bit to 64-bit.
    __m256i vidx = _mm256_setr_epi64x(delta_int[y + 0],
                                      delta_int[y + 1],
                                      delta_int[y + 2],
                                      delta_int[y + 3]);
    __m256i all_weights = _mm256_loadu_si256((__m256i*)f);
    __m256i w01 = _mm256_shuffle_epi8(all_weights, w_shuf_01);
    __m256i w23 = _mm256_shuffle_epi8(all_weights, w_shuf_23);

    for (int_fast32_t x = 0; x + 3 < width; x += 4, p += 4) {

      __m256i vp = _mm256_i64gather_epi64((const long long int*)p, vidx, 1);
      __m256i vp_01 = _mm256_shuffle_epi8(vp, p_shuf_01);
      __m256i vp_23 = _mm256_shuffle_epi8(vp, p_shuf_23);

      __m256i dot_01 = _mm256_maddubs_epi16(vp_01, w01);
      __m256i dot_23 = _mm256_maddubs_epi16(vp_23, w23);
      __m256i sum = _mm256_add_epi16(dot_01, dot_23);
      sum = _mm256_add_epi16(sum, _mm256_set1_epi16(32));
      sum = _mm256_srai_epi16(sum, 6);

      __m128i lo = _mm256_castsi256_si128(sum);
      __m128i hi = _mm256_extracti128_si256(sum, 1);
      __m128i filtered = _mm_packus_epi16(lo, hi);

      _mm_storeu_si128((__m128i*)(dst + (y * width + x)), filtered);
    }
  }
}

static void angular_pred_w8_ver_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t* delta_int, const int16_t* delta_fract, const int height, const int use_cubic)
{
  const int width = 8;

  const __m128i p_shuf_01 = _mm_setr_epi8(
    0x00, 0x01, 0x01, 0x02, 0x02, 0x03, 0x03, 0x04,
    0x04, 0x05, 0x05, 0x06, 0x06, 0x07, 0x07, 0x08
  );

  const __m128i p_shuf_23 = _mm_setr_epi8(
    0x02, 0x03, 0x03, 0x04, 0x04, 0x05, 0x05, 0x06,
    0x06, 0x07, 0x07, 0x08, 0x08, 0x09, 0x09, 0x0a
  );

  const __m256i w_shuf_01 = _mm256_setr_epi8(
    0x00, 0x02, 0x00, 0x02, 0x00, 0x02, 0x00, 0x02,
    0x00, 0x02, 0x00, 0x02, 0x00, 0x02, 0x00, 0x02,
    0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a,
    0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a
  );

  const __m256i w_shuf_23 = _mm256_setr_epi8(
    0x04, 0x06, 0x04, 0x06, 0x04, 0x06, 0x04, 0x06,
    0x04, 0x06, 0x04, 0x06, 0x04, 0x06, 0x04, 0x06,
    0x0c, 0x0e, 0x0c, 0x0e, 0x0c, 0x0e, 0x0c, 0x0e,
    0x0c, 0x0e, 0x0c, 0x0e, 0x0c, 0x0e, 0x0c, 0x0e
  );

  // For a 8 width block, height must be at least 2. Handle 2 lines at once
  for (int y = 0; y < height; y += 2) {
    __m256i all_weights;
    if (use_cubic) {
      int16_t tmp[8];
      memcpy(&tmp[0], cubic_filter[delta_fract[y + 0]], 4 * sizeof(int16_t));
      memcpy(&tmp[4], cubic_filter[delta_fract[y + 1]], 4 * sizeof(int16_t));
      all_weights = _mm256_setr_epi64x(*(int64_t*)&tmp[0], *(int64_t*)&tmp[4], *(int64_t*)&tmp[0], *(int64_t*)&tmp[4]);
    }
    else {
      int16_t tmp[8];
      for (int yy = 0; yy < 2; ++yy) {
        const int16_t offset = (delta_fract[y + yy] >> 1);
        const int idx = yy * 4;
        tmp[idx + 0] = 16 - offset;
        tmp[idx + 1] = 32 - offset;
        tmp[idx + 2] = 16 + offset;
        tmp[idx + 3] = offset;
      }
      all_weights = _mm256_setr_epi64x(*(int64_t*)&tmp[0], *(int64_t*)&tmp[4], *(int64_t*)&tmp[0], *(int64_t*)&tmp[4]);
    }

    // Do 4-tap intra interpolation filtering
    uvg_pixel* p = (uvg_pixel*)ref_main;

    // Weights are 16-bit, but shuffle will cut out the unnecessary bits.
    __m256i w01 = _mm256_shuffle_epi8(all_weights, w_shuf_01);
    __m256i w23 = _mm256_shuffle_epi8(all_weights, w_shuf_23);

    for (int_fast32_t x = 0; x < width; x += 8, p += 8) {
      __m128i vp0 = _mm_loadu_si128((__m128i*)(p + delta_int[y + 0]));
      __m128i vp1 = _mm_loadu_si128((__m128i*)(p + delta_int[y + 1]));

      __m256i vp_01 = _mm256_castsi128_si256(_mm_shuffle_epi8(vp0, p_shuf_01));
      vp_01 = _mm256_inserti128_si256(vp_01, _mm_shuffle_epi8(vp1, p_shuf_01), 1);

      __m256i vp_23 = _mm256_castsi128_si256(_mm_shuffle_epi8(vp0, p_shuf_23));
      vp_23 = _mm256_inserti128_si256(vp_23, _mm_shuffle_epi8(vp1, p_shuf_23), 1);

      __m256i dot_01 = _mm256_maddubs_epi16(vp_01, w01);
      __m256i dot_23 = _mm256_maddubs_epi16(vp_23, w23);
      __m256i sum = _mm256_add_epi16(dot_01, dot_23);
      sum = _mm256_add_epi16(sum, _mm256_set1_epi16(32));
      sum = _mm256_srai_epi16(sum, 6);

      __m128i lo = _mm256_castsi256_si128(sum);
      __m128i hi = _mm256_extracti128_si256(sum, 1);
      __m128i filtered = _mm_packus_epi16(lo, hi);

      _mm_store_si128((__m128i*)(dst + (y * 8)), filtered);
    }
  }
}

static void angular_pred_w16_ver_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t* delta_int, const int16_t* delta_fract, const int width, const int height, const int use_cubic)
{
  const __m256i p_shuf_01 = _mm256_setr_epi8(
    0x00, 0x01, 0x01, 0x02, 0x02, 0x03, 0x03, 0x04,
    0x04, 0x05, 0x05, 0x06, 0x06, 0x07, 0x07, 0x08,
    0x00, 0x01, 0x01, 0x02, 0x02, 0x03, 0x03, 0x04,
    0x04, 0x05, 0x05, 0x06, 0x06, 0x07, 0x07, 0x08
  );

  const __m256i p_shuf_23 = _mm256_setr_epi8(
    0x02, 0x03, 0x03, 0x04, 0x04, 0x05, 0x05, 0x06,
    0x06, 0x07, 0x07, 0x08, 0x08, 0x09, 0x09, 0x0a,
    0x02, 0x03, 0x03, 0x04, 0x04, 0x05, 0x05, 0x06,
    0x06, 0x07, 0x07, 0x08, 0x08, 0x09, 0x09, 0x0a
  );

  const __m256i w_shuf_01 = _mm256_setr_epi8(
    0x00, 0x02, 0x00, 0x02, 0x00, 0x02, 0x00, 0x02,
    0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a,
    0x00, 0x02, 0x00, 0x02, 0x00, 0x02, 0x00, 0x02,
    0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a
  );

  const __m256i w_shuf_23 = _mm256_setr_epi8(
    0x04, 0x06, 0x04, 0x06, 0x04, 0x06, 0x04, 0x06,
    0x0c, 0x0e, 0x0c, 0x0e, 0x0c, 0x0e, 0x0c, 0x0e,
    0x04, 0x06, 0x04, 0x06, 0x04, 0x06, 0x04, 0x06,
    0x0c, 0x0e, 0x0c, 0x0e, 0x0c, 0x0e, 0x0c, 0x0e
  );

  //int16_t f[4][4] = { { 0 } };

  // For a 16 width block, height can be 1.
  for (int y = 0; y < height; ++y) {
    __m256i all_weights;
    if (use_cubic) {
      //memcpy(f[0], cubic_filter[delta_fract[y + 0]], 8);
      //memcpy(f[1], cubic_filter[delta_fract[y + 1]], 8);
      //memcpy(f[2], cubic_filter[delta_fract[y + 2]], 8);
      //memcpy(f[3], cubic_filter[delta_fract[y + 3]], 8);
      //int64_t *tmp = (int64_t*)&delta_fract[y];
      int16_t tmp[4];
      memcpy(&tmp, cubic_filter[delta_fract[y]], 8);
      all_weights = _mm256_set1_epi64x(*(int64_t*)tmp);
    }
    else {
      const int16_t offset = (delta_fract[y] >> 1);
      int16_t tmp[4];
      tmp[0] = 16 - offset;
      tmp[1] = 32 - offset;
      tmp[2] = 16 + offset;
      tmp[3] = offset;
      all_weights = _mm256_set1_epi64x(*(int64_t*)tmp);
    }

    // Do 4-tap intra interpolation filtering
    uvg_pixel* p = (uvg_pixel*)ref_main;
    // This solution assumes the delta int values to be 64-bit
    // Cast from 16-bit to 64-bit.
    __m256i vidx = _mm256_setr_epi64x(delta_int[y] + 0,
                                      delta_int[y] + 4,
                                      delta_int[y] + 8,
                                      delta_int[y] + 12);

    //__m256i all_weights = _mm256_loadu_si256((__m256i*)f);

    __m256i w01 = _mm256_shuffle_epi8(all_weights, w_shuf_01);
    __m256i w23 = _mm256_shuffle_epi8(all_weights, w_shuf_23);

    for (int_fast32_t x = 0; x < width; x += 16, p += 16) {
      __m256i vp = _mm256_loadu_si256((__m256i*)(p + delta_int[y]));

      __m256i tmp = _mm256_permute4x64_epi64(vp, _MM_SHUFFLE(2, 1, 1, 0));

      __m256i vp_01 = _mm256_shuffle_epi8(tmp, p_shuf_01);
      __m256i vp_23 = _mm256_shuffle_epi8(tmp, p_shuf_23);

      __m256i dot_01 = _mm256_maddubs_epi16(vp_01, w01);
      __m256i dot_23 = _mm256_maddubs_epi16(vp_23, w23);
      __m256i sum = _mm256_add_epi16(dot_01, dot_23);
      sum = _mm256_add_epi16(sum, _mm256_set1_epi16(32));
      sum = _mm256_srai_epi16(sum, 6);

      __m128i lo = _mm256_castsi256_si128(sum);
      __m128i hi = _mm256_extracti128_si256(sum, 1);
      __m128i filtered = _mm_packus_epi16(lo, hi);

      _mm_store_si128((__m128i*)dst, filtered);
      dst += 16;
    }
  }
}


static void angular_pred_w4_hor_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t* delta_int, const int16_t* delta_fract, const int height, const int use_cubic)
{
  const int width = 4;

  const __m256i p_shuf_01 = _mm256_setr_epi8(
    0x00, 0x01, 0x08, 0x09, 0x01, 0x02, 0x09, 0x0a,
    0x02, 0x03, 0x0a, 0x0b, 0x03, 0x04, 0x0b, 0x0c,
    0x00, 0x01, 0x08, 0x09, 0x01, 0x02, 0x09, 0x0a,
    0x02, 0x03, 0x0a, 0x0b, 0x03, 0x04, 0x0b, 0x0c
  );

  const __m256i p_shuf_23 = _mm256_setr_epi8(
    0x02, 0x03, 0x0a, 0x0b, 0x03, 0x04, 0x0b, 0x0c,
    0x04, 0x05, 0x0c, 0x0d, 0x05, 0x06, 0x0d, 0x0e,
    0x02, 0x03, 0x0a, 0x0b, 0x03, 0x04, 0x0b, 0x0c,
    0x04, 0x05, 0x0c, 0x0d, 0x05, 0x06, 0x0d, 0x0e
  );

  const __m256i w_shuf_01 = _mm256_setr_epi8(
    0x00, 0x02, 0x08, 0x0a, 0x00, 0x02, 0x08, 0x0a,
    0x00, 0x02, 0x08, 0x0a, 0x00, 0x02, 0x08, 0x0a,
    0x00, 0x02, 0x08, 0x0a, 0x00, 0x02, 0x08, 0x0a,
    0x00, 0x02, 0x08, 0x0a, 0x00, 0x02, 0x08, 0x0a
  );

  const __m256i w_shuf_23 = _mm256_setr_epi8(
    0x04, 0x06, 0x0c, 0x0e, 0x04, 0x06, 0x0c, 0x0e,
    0x04, 0x06, 0x0c, 0x0e, 0x04, 0x06, 0x0c, 0x0e,
    0x04, 0x06, 0x0c, 0x0e, 0x04, 0x06, 0x0c, 0x0e,
    0x04, 0x06, 0x0c, 0x0e, 0x04, 0x06, 0x0c, 0x0e
  );

  const __m128i r_shuffle = _mm_setr_epi8(
    0x00, 0x01, 0x08, 0x09, 0x02, 0x03, 0x0a, 0x0b,
    0x04, 0x05, 0x0c, 0x0d, 0x06, 0x07, 0x0e, 0x0f
  );

  int16_t f[4][4] = { { 0 } };

  // For a 4 width block, height must be at least 4. Handle 4 lines at once
  for (int y = 0; y < height; y += 4) {
    
    // Do 4-tap intra interpolation filtering
    uvg_pixel* p = (uvg_pixel*)(ref_main + y);

    for (int_fast32_t x = 0; x < width; x += 4) {
      if (use_cubic) {
        memcpy(f[0], cubic_filter[delta_fract[x + 0]], 8);
        memcpy(f[1], cubic_filter[delta_fract[x + 1]], 8);
        memcpy(f[2], cubic_filter[delta_fract[x + 2]], 8);
        memcpy(f[3], cubic_filter[delta_fract[x + 3]], 8);
      }
      else {
        for (int xx = 0; xx < 4; ++xx) {
          const int16_t offset = (delta_fract[x + xx] >> 1);
          f[xx][0] = 16 - offset;
          f[xx][1] = 32 - offset;
          f[xx][2] = 16 + offset;
          f[xx][3] = offset;
        }
      }

      // This solution assumes the delta int values to be 64-bit
      // Cast from 16-bit to 64-bit.
      __m256i vidx = _mm256_setr_epi64x(delta_int[x + 0],
                                        delta_int[x + 1],
                                        delta_int[x + 2],
                                        delta_int[x + 3]);
      __m256i all_weights = _mm256_loadu_si256((__m256i*)f);
      __m256i w01 = _mm256_shuffle_epi8(all_weights, w_shuf_01);
      __m256i w23 = _mm256_shuffle_epi8(all_weights, w_shuf_23);

      __m256i vp = _mm256_i64gather_epi64((const long long int*)p, vidx, 1);
      //__m256i vp = _mm256_loadu_si256((__m256i*)(p + delta_int[y]));

      //__m256i tmp = _mm256_permute4x64_epi64(vp, _MM_SHUFFLE(2, 1, 1, 0));

      __m256i vp_01 = _mm256_shuffle_epi8(vp, p_shuf_01);
      __m256i vp_23 = _mm256_shuffle_epi8(vp, p_shuf_23);

      __m256i dot_01 = _mm256_maddubs_epi16(vp_01, w01);
      __m256i dot_23 = _mm256_maddubs_epi16(vp_23, w23);
      __m256i sum = _mm256_add_epi16(dot_01, dot_23);
      sum = _mm256_add_epi16(sum, _mm256_set1_epi16(32));
      sum = _mm256_srai_epi16(sum, 6);

      __m128i lo = _mm256_castsi256_si128(sum);
      __m128i hi = _mm256_extracti128_si256(sum, 1);
      __m128i filtered = _mm_packus_epi16(lo, hi);

      _mm_store_si128((__m128i*)dst, _mm_shuffle_epi8(filtered, r_shuffle));
      dst += 16;
    }
  }
}

static void angular_pred_w8_hor_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t* delta_int, const int16_t* delta_fract, const int height, const int use_cubic)
{
  const int width = 8;

  int8_t f[8][4] = { { 0 } };
  if (use_cubic) {
    memcpy(f[0], cubic_filter_8bit_c[delta_fract[0]], sizeof(int8_t) * 4);
    memcpy(f[1], cubic_filter_8bit_c[delta_fract[1]], sizeof(int8_t) * 4);
    memcpy(f[2], cubic_filter_8bit_c[delta_fract[2]], sizeof(int8_t) * 4);
    memcpy(f[3], cubic_filter_8bit_c[delta_fract[3]], sizeof(int8_t) * 4);
    memcpy(f[4], cubic_filter_8bit_c[delta_fract[4]], sizeof(int8_t) * 4);
    memcpy(f[5], cubic_filter_8bit_c[delta_fract[5]], sizeof(int8_t) * 4);
    memcpy(f[6], cubic_filter_8bit_c[delta_fract[6]], sizeof(int8_t) * 4);
    memcpy(f[7], cubic_filter_8bit_c[delta_fract[7]], sizeof(int8_t) * 4);
  }
  else {
    for (int x = 0; x < 8; ++x) {
      const int8_t offset = (delta_fract[x] >> 1);
      f[x][0] = 16 - offset;
      f[x][1] = 32 - offset;
      f[x][2] = 16 + offset;
      f[x][3] = offset;
    }
  }

  __m128i tmp = _mm_load_si128((__m128i*)delta_int);
  __m256i vidx = _mm256_cvtepi16_epi32(tmp);
  __m256i weights = _mm256_loadu_si256((__m256i*)f);

  for (int y = 0; y < height; y += 2) {

    // Do 4-tap intra interpolation filtering
    uvg_pixel* p = (uvg_pixel*)(ref_main + y);
    __m256i vp0 = _mm256_i32gather_epi32((const int*)(p + 0), vidx, 1);
    __m256i vp1 = _mm256_i32gather_epi32((const int*)(p + 1), vidx, 1);

    __m256i dot_01 = _mm256_maddubs_epi16(vp0, weights);
    __m256i dot_23 = _mm256_maddubs_epi16(vp1, weights);
    __m256i sum = _mm256_hadd_epi16(dot_01, dot_23);
    sum = _mm256_add_epi16(sum, _mm256_set1_epi16(32));
    sum = _mm256_srai_epi16(sum, 6);

    __m128i lo = _mm256_castsi256_si128(sum);
    __m128i hi = _mm256_extracti128_si256(sum, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);
    filtered = _mm_shuffle_epi32(filtered, _MM_SHUFFLE(3, 1, 2, 0));
      
    _mm_store_si128((__m128i*)dst,  filtered);
 
    dst += 16;
  }
}

static void angular_pred_w16_hor_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t* delta_int, const int16_t* delta_fract, const int width, const int height, const int use_cubic)
{
  int8_t f[64][4] = { { 0 } };
  if (use_cubic) {
    for (int x = 0; x < width; ++x) {
      memcpy(f[x], cubic_filter_8bit_c[delta_fract[x]], sizeof(int8_t) * 4);
    }
  }
  else {
    for (int x = 0; x < width; ++x) {
      const int8_t offset = (delta_fract[x] >> 1);
      f[x][0] = 16 - offset;
      f[x][1] = 32 - offset;
      f[x][2] = 16 + offset;
      f[x][3] = offset;
    }
  }

  for (int x = 0; x < width; x += 16) {
    __m128i tmp0 = _mm_load_si128((__m128i*)&delta_int[x]);
    __m128i tmp1 = _mm_load_si128((__m128i*)&delta_int[x + 8]);
    __m256i vidx0 = _mm256_cvtepi16_epi32(tmp0);
    __m256i vidx1 = _mm256_cvtepi16_epi32(tmp1);

    __m256i w0 = _mm256_loadu_si256((__m256i*) & f[x + 0]);
    __m256i w1 = _mm256_loadu_si256((__m256i*) & f[x + 8]);

    // Width 16, handle one row at a time
    for (int y = 0; y < height; ++y) {

      // Do 4-tap intra interpolation filtering
      uvg_pixel* p = (uvg_pixel*)(ref_main + y);
      __m256i vp0 = _mm256_i32gather_epi32((const int*)p, vidx0, 1);
      __m256i vp1 = _mm256_i32gather_epi32((const int*)p, vidx1, 1);

      __m256i dot_01 = _mm256_maddubs_epi16(vp0, w0);
      __m256i dot_23 = _mm256_maddubs_epi16(vp1, w1);
      __m256i sum = _mm256_hadd_epi16(dot_01, dot_23);
      sum = _mm256_add_epi16(sum, _mm256_set1_epi16(32));
      sum = _mm256_srai_epi16(sum, 6);

      __m128i lo = _mm256_castsi256_si128(sum);
      __m128i hi = _mm256_extracti128_si256(sum, 1);
      __m128i filtered = _mm_packus_epi16(lo, hi);
      filtered = _mm_shuffle_epi32(filtered, _MM_SHUFFLE(3, 1, 2, 0));

      _mm_store_si128((__m128i*)(dst + (y * width + x)), filtered);
    }
  }
}


static void angular_pred_generic_linear_filter(uvg_pixel* dst, uvg_pixel* ref, const int width, const int height, const int16_t* delta_int, const int16_t* delta_fract)
{
  // 2-tap filter
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      uvg_pixel ref1 = ref[x + delta_int[y] + 1];
      uvg_pixel ref2 = ref[x + delta_int[y] + 2];
      //dst[y * width + x] = ref1 + ((delta_fract[y] * (ref2 - ref1) + 16) >> 5);
      dst[y * width + x] = ((32 - delta_fract[y]) * ref1 + delta_fract[y] * ref2 + 16) >> 5;
    }
  }
}


// Linear interpolation filter for width 4 has a different call, since it uses premade tables for coefficients
static void angular_pred_linear_filter_w4_ver_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int16_t* delta_int, const int32_t pred_mode)
{
  const int16_t* dint = delta_int;
  const __m128i v16s = _mm_set1_epi16(16);

  const int mode_idx = pred_mode <= 34 ? pred_mode - 2 : 66 - pred_mode;
  const int weight_table_offset = coeff_table_mode_offsets[mode_idx];
  const int vnum = coeff_vector128_num_by_mode[mode_idx];
  const int modulo = vnum - 1;
  int offset_num = 0;

  ALIGNED(16) int16_t shuffle_vector_offsets[8];
  memcpy(shuffle_vector_offsets, &intra_chroma_linear_interpolation_w4_ver_shuffle_vector_offset[mode_idx * 8], sizeof(int16_t) * 8);

  // Height has to be at least 4, handle 4 lines at once
  for (int y = 0; y < height; y += 4) {
    // Load refs from smallest index onwards, shuffle will handle the rest. The smallest index will be at one of these delta int table indices
    const int16_t min_offset = 1 + MIN(dint[0], dint[3]);
    dint += 4;
    // Load enough reff samples to cover four 4 width lines. Shuffles will put the samples in correct places.
    const __m128i vsrc_raw = _mm_loadu_si128((const __m128i*) & ref[min_offset]);
    const int offset = weight_table_offset + (offset_num * 16);
    
    const __m128i vcoeff0 = _mm_load_si128((const __m128i*)&intra_chroma_linear_interpolation_weights_w4_ver[offset]);
    const __m128i vcoeff1 = vnum == 1 ? vcoeff0 : _mm_load_si128((const __m128i*)&intra_chroma_linear_interpolation_weights_w4_ver[offset + 16]);

    const __m128i vshuf0 = _mm_load_si128((const __m128i*) &intra_chroma_linear_interpolation_shuffle_vectors_w4_ver[shuffle_vector_offsets[y >> 2] + 0]);
    const __m128i vshuf1 = _mm_load_si128((const __m128i*) &intra_chroma_linear_interpolation_shuffle_vectors_w4_ver[shuffle_vector_offsets[y >> 2] + 16]);

    __m128i vsrc0 = _mm_shuffle_epi8(vsrc_raw, vshuf0);
    __m128i vsrc1 = _mm_shuffle_epi8(vsrc_raw, vshuf1);
    
    __m128i res0 = _mm_maddubs_epi16(vsrc0, vcoeff0);
    __m128i res1 = _mm_maddubs_epi16(vsrc1, vcoeff1);
    res0 = _mm_add_epi16(res0, v16s);
    res1 = _mm_add_epi16(res1, v16s);
    res0 = _mm_srai_epi16(res0, 5);
    res1 = _mm_srai_epi16(res1, 5);

    _mm_store_si128((__m128i*)dst, _mm_packus_epi16(res0, res1));
    dst += 16;
    offset_num += 2;
    // This resets the offset number to 0 when it reaches the end of the table. Only works on powers of 2.
    offset_num &= modulo;
  }
}


static void angular_pred_linear_filter_w8_ver_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int16_t* delta_int, const int pred_mode)
{
  const int width = 8;
  const __m128i v16s = _mm_set1_epi16(16);
  const __m128i vshuf = _mm_setr_epi8(
    0x00, 0x01, 0x01, 0x02, 0x02, 0x03, 0x03, 0x04,
    0x04, 0x05, 0x05, 0x06, 0x06, 0x07, 0x07, 0x08
  );
  const int wide_angle = pred_mode > 66 || pred_mode < 2;
  const int mode_idx = wide_angle ? (pred_mode < 2 ? 12 + pred_mode : 80 - pred_mode) : (pred_mode <= 34 ? (pred_mode - 2) : (66 - pred_mode));
  const int8_t* coeff_table = wide_angle ? intra_chroma_linear_interpolation_weights_w8_ver_wide_angle : intra_chroma_linear_interpolation_weights_w8_ver;
  const int coeff_table_offset = mode_idx * 64;

  // Height has to be at least 2, handle 2 lines at once
  for (int y = 0; y < height; y += 2) {
    const int16_t* coeff_tmp0 = (const int16_t*) &coeff_table[coeff_table_offset + (y << 1) + 0];
    const int16_t* coeff_tmp1 = (const int16_t*) &coeff_table[coeff_table_offset + (y << 1) + 2];
    
    __m128i vsrc0 = _mm_loadu_si128((const __m128i*) & ref[delta_int[y + 0] + 1]);
    __m128i vsrc1 = _mm_loadu_si128((const __m128i*) & ref[delta_int[y + 1] + 1]);

    vsrc0 = _mm_shuffle_epi8(vsrc0, vshuf);
    vsrc1 = _mm_shuffle_epi8(vsrc1, vshuf);

    const __m128i vcoeff0 = _mm_set1_epi16(*coeff_tmp0);
    const __m128i vcoeff1 = _mm_set1_epi16(*coeff_tmp1);

    __m128i res0 = _mm_maddubs_epi16(vsrc0, vcoeff0);
    __m128i res1 = _mm_maddubs_epi16(vsrc1, vcoeff1);
    res0 = _mm_add_epi16(res0, v16s);
    res1 = _mm_add_epi16(res1, v16s);
    res0 = _mm_srai_epi16(res0, 5);
    res1 = _mm_srai_epi16(res1, 5);

    _mm_store_si128((__m128i*)dst, _mm_packus_epi16(res0, res1));
    dst += 16;
  }
}


static void angular_pred_linear_filter_w16_ver_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int16_t* delta_int, const int pred_mode)
{
  const __m128i v16s = _mm_set1_epi16(16);
  const __m128i vshuf = _mm_setr_epi8(
    0x00, 0x01, 0x01, 0x02, 0x02, 0x03, 0x03, 0x04,
    0x04, 0x05, 0x05, 0x06, 0x06, 0x07, 0x07, 0x08
  );

  const int wide_angle = pred_mode > 66 || pred_mode < 2;
  const int mode_idx = wide_angle ? (pred_mode < 2 ? 12 + pred_mode : 80 - pred_mode) : (pred_mode <= 34 ? (pred_mode - 2) : (66 - pred_mode));
  const int8_t* coeff_table = wide_angle ? intra_chroma_linear_interpolation_weights_w8_ver_wide_angle : intra_chroma_linear_interpolation_weights_w8_ver;
  const int coeff_table_offset = mode_idx * 64;

  // Handle 1 line at a time
  for (int y = 0; y < height; ++y) {
    const int16_t* coeff_tmp = (const int16_t*)&coeff_table[coeff_table_offset + (y << 1)];
    __m128i vcoeff = _mm_set1_epi16(*coeff_tmp);

    __m128i vsrc0 = _mm_loadu_si128((const __m128i*)&ref[delta_int[y] + 0 + 1]);
    __m128i vsrc1 = _mm_loadu_si128((const __m128i*)&ref[delta_int[y] + 8 + 1]);
    
    vsrc0 = _mm_shuffle_epi8(vsrc0, vshuf);
    vsrc1 = _mm_shuffle_epi8(vsrc1, vshuf);

    __m128i res0 = _mm_maddubs_epi16(vsrc0, vcoeff);
    __m128i res1 = _mm_maddubs_epi16(vsrc1, vcoeff);
    res0 = _mm_add_epi16(res0, v16s);
    res1 = _mm_add_epi16(res1, v16s);
    res0 = _mm_srai_epi16(res0, 5);
    res1 = _mm_srai_epi16(res1, 5);

    _mm_store_si128((__m128i*)dst, _mm_packus_epi16(res0, res1));
    dst += 16;
  }
}


static void angular_pred_linear_filter_w32_ver_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int16_t* delta_int, const int pred_mode)
{
  const __m256i v16s = _mm256_set1_epi16(16);
  const __m256i vshuf = _mm256_setr_epi8(
    0x00, 0x01, 0x01, 0x02, 0x02, 0x03, 0x03, 0x04,
    0x04, 0x05, 0x05, 0x06, 0x06, 0x07, 0x07, 0x08,
    0x00, 0x01, 0x01, 0x02, 0x02, 0x03, 0x03, 0x04,
    0x04, 0x05, 0x05, 0x06, 0x06, 0x07, 0x07, 0x08
  );

  const int wide_angle = pred_mode > 66 || pred_mode < 2;
  const int mode_idx = wide_angle ? (pred_mode < 2 ? 12 + pred_mode : 80 - pred_mode) : (pred_mode <= 34 ? (pred_mode - 2) : (66 - pred_mode));
  const int8_t* coeff_table = wide_angle ? intra_chroma_linear_interpolation_weights_w8_ver_wide_angle : intra_chroma_linear_interpolation_weights_w8_ver;
  const int coeff_table_offset = mode_idx * 64;

  // Handle 1 line at a time
  for (int y = 0; y < height; ++y) {
    const int16_t* coeff_tmp = (const int16_t*)&coeff_table[coeff_table_offset + (y << 1)];
    __m256i vcoeff = _mm256_set1_epi16(*coeff_tmp);

    ALIGNED(32) __m128i vsrc[4];
    vsrc[0] = _mm_loadu_si128((const __m128i*) & ref[delta_int[y] + 0 + 1]);
    vsrc[1] = _mm_loadu_si128((const __m128i*) & ref[delta_int[y] + 16 + 1]); // Flip these two middle sources. They will be later flipped back into place by packus
    vsrc[2] = _mm_loadu_si128((const __m128i*) & ref[delta_int[y] + 8 + 1]);
    vsrc[3] = _mm_loadu_si128((const __m128i*) & ref[delta_int[y] + 24 + 1]);

    __m256i* vsrc256 = (__m256i*)vsrc;
    vsrc256[0] = _mm256_shuffle_epi8(vsrc256[0], vshuf);
    vsrc256[1] = _mm256_shuffle_epi8(vsrc256[1], vshuf);

    __m256i res0 = _mm256_maddubs_epi16(vsrc256[0], vcoeff);
    __m256i res1 = _mm256_maddubs_epi16(vsrc256[1], vcoeff);
    res0 = _mm256_add_epi16(res0, v16s);
    res1 = _mm256_add_epi16(res1, v16s);
    res0 = _mm256_srai_epi16(res0, 5);
    res1 = _mm256_srai_epi16(res1, 5);

    _mm256_store_si256((__m256i*)dst, _mm256_packus_epi16(res0, res1));
    dst += 32;
  }
}


static void angular_pred_linear_filter_w4_hor_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int mode, const int16_t* delta_int)
{
  const int16_t* dint = delta_int;
  const __m128i v16s = _mm_set1_epi16(16);

  const int16_t weigth_offset = mode_to_weight_table_offset_w4_hor[mode];
  const int16_t shuf_offset = mode_to_shuffle_vector_table_offset_w4_hor[mode];

  __m128i vcoeff = _mm_load_si128((const __m128i*) &intra_chroma_linear_interpolation_weights_w4_hor[weigth_offset]);
  __m128i vshuf0 = _mm_load_si128((const __m128i*) &intra_chroma_linear_interpolation_shuffle_vectors_w4_hor[shuf_offset + 0]);
  __m128i vshuf1 = _mm_load_si128((const __m128i*) &intra_chroma_linear_interpolation_shuffle_vectors_w4_hor[shuf_offset + 16]);

  // Load refs from smallest index onwards, shuffle will handle the rest. The smallest index will be at one of these delta int table indices
  const int16_t min_offset = 1 + MIN(dint[0], dint[3]);

  // Height has to be at least 4, handle 4 lines at once
  for (int y = 0; y < height; y += 4) {
    // Prepare sources
    __m128i vidx = _mm_set_epi64x((long long int)(min_offset + y + 2), (long long int)(min_offset + y + 0));
    __m128i vsrc_tmp = _mm_i64gather_epi64((const long long*)ref, vidx, 1);
    __m128i vsrc0 = _mm_shuffle_epi8(vsrc_tmp, vshuf0);
    __m128i vsrc1 = _mm_shuffle_epi8(vsrc_tmp, vshuf1);
  
    __m128i res0 = _mm_maddubs_epi16(vsrc0, vcoeff);
    __m128i res1 = _mm_maddubs_epi16(vsrc1, vcoeff);
    res0 = _mm_add_epi16(res0, v16s);
    res1 = _mm_add_epi16(res1, v16s);
    res0 = _mm_srai_epi16(res0, 5);
    res1 = _mm_srai_epi16(res1, 5);

    _mm_store_si128((__m128i*)dst, _mm_packus_epi16(res0, res1));
    dst += 16;
  }
}


static void angular_pred_linear_filter_w8_hor_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int mode, const int16_t* delta_int)
{
  const int16_t* dint = delta_int;
  const __m128i v16s = _mm_set1_epi16(16);
  const int16_t weigth_offset = (mode - 2) * 16;
  const int16_t shuf_offset = (mode - 2) * 32;

  __m128i vcoeff = _mm_load_si128((const __m128i*) & intra_chroma_linear_interpolation_weights_w8_hor[weigth_offset]);
  __m128i vshuf0 = _mm_load_si128((const __m128i*) & intra_chroma_linear_interpolation_shuffle_vectors_w8_hor[shuf_offset + 0]);
  __m128i vshuf1 = _mm_load_si128((const __m128i*) & intra_chroma_linear_interpolation_shuffle_vectors_w8_hor[shuf_offset + 16]);

  // Load refs from smallest index onwards, shuffle will handle the rest. The smallest index will be at one of these delta int table indices
  const int16_t min_offset = 1 + MIN(dint[0], dint[7]);

  // Height has to be at least 2, handle 2 lines at once
  for (int y = 0; y < height; y += 2) {
    // Prepare sources
    __m128i vsrc_tmp = _mm_loadu_si128((__m128i*)&ref[min_offset + y]);
    const __m128i vsrc0 = _mm_shuffle_epi8(vsrc_tmp, vshuf0);
    const __m128i vsrc1 = _mm_shuffle_epi8(vsrc_tmp, vshuf1);

    __m128i res0 = _mm_maddubs_epi16(vsrc0, vcoeff);
    __m128i res1 = _mm_maddubs_epi16(vsrc1, vcoeff);
    res0 = _mm_add_epi16(res0, v16s);
    res1 = _mm_add_epi16(res1, v16s);
    res0 = _mm_srai_epi16(res0, 5);
    res1 = _mm_srai_epi16(res1, 5);

    _mm_store_si128((__m128i*)dst, _mm_packus_epi16(res0, res1));
    dst += 16;
  }
}


static void angular_pred_linear_filter_w16_hor_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int mode, const int16_t* delta_int)
{
  const int16_t* dint = delta_int;
  const __m128i v16s = _mm_set1_epi16(16);
  const int16_t weigth_offset = (mode - 2) * 32;
  const int16_t shuf_offset = (mode - 2) * 64;

  __m128i vcoeff0 = _mm_load_si128((const __m128i*) & intra_chroma_linear_interpolation_weights_w16_hor[weigth_offset + 0]);
  __m128i vcoeff1 = _mm_load_si128((const __m128i*) & intra_chroma_linear_interpolation_weights_w16_hor[weigth_offset + 16]);
  __m128i vshuf0 = _mm_load_si128((const __m128i*) & intra_chroma_linear_interpolation_shuffle_vectors_w16_hor[shuf_offset + 0]);
  __m128i vshuf1 = _mm_load_si128((const __m128i*) & intra_chroma_linear_interpolation_shuffle_vectors_w16_hor[shuf_offset + 16]);
  __m128i vshuf2 = _mm_load_si128((const __m128i*) & intra_chroma_linear_interpolation_shuffle_vectors_w16_hor[shuf_offset + 32]);
  __m128i vshuf3 = _mm_load_si128((const __m128i*) & intra_chroma_linear_interpolation_shuffle_vectors_w16_hor[shuf_offset + 48]);

  // Load refs from smallest index onwards, shuffle will handle the rest. The smallest index will be at one of these delta int table indices
  const int16_t min_offset0 = 1 + MIN(dint[0], dint[7]);
  const int16_t min_offset1 = 1 + MIN(dint[8], dint[15]);

  // Height has to be at least 2, there is no 16x1 block for chroma.
  for (int y = 0; y < height; y += 2) {
    // Prepare sources
    __m128i vsrc_tmp0 = _mm_loadu_si128((__m128i*) &ref[min_offset0 + y]);
    __m128i vsrc_tmp1 = _mm_loadu_si128((__m128i*) &ref[min_offset1 + y]);
    const __m128i vsrc0 = _mm_shuffle_epi8(vsrc_tmp0, vshuf0);
    const __m128i vsrc1 = _mm_shuffle_epi8(vsrc_tmp1, vshuf1);
    const __m128i vsrc2 = _mm_shuffle_epi8(vsrc_tmp0, vshuf2);
    const __m128i vsrc3 = _mm_shuffle_epi8(vsrc_tmp1, vshuf3);

    __m128i res0 = _mm_maddubs_epi16(vsrc0, vcoeff0);
    __m128i res1 = _mm_maddubs_epi16(vsrc1, vcoeff1);
    __m128i res2 = _mm_maddubs_epi16(vsrc2, vcoeff0);
    __m128i res3 = _mm_maddubs_epi16(vsrc3, vcoeff1);
    res0 = _mm_add_epi16(res0, v16s);
    res1 = _mm_add_epi16(res1, v16s);
    res2 = _mm_add_epi16(res2, v16s);
    res3 = _mm_add_epi16(res3, v16s);
    res0 = _mm_srai_epi16(res0, 5);
    res1 = _mm_srai_epi16(res1, 5);
    res2 = _mm_srai_epi16(res2, 5);
    res3 = _mm_srai_epi16(res3, 5);

    _mm_store_si128((__m128i*)&dst[0], _mm_packus_epi16(res0, res1));
    _mm_store_si128((__m128i*)&dst[16], _mm_packus_epi16(res2, res3));
    dst += 32;
  }
}


static void angular_pred_linear_filter_w32_hor_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int mode, const int16_t* delta_int)
{
  const int16_t* dint = delta_int;
  const __m128i v16s = _mm_set1_epi16(16);
  const int16_t weigth_offset = (mode - 2) * 64;
  const int16_t shuf_offset = (mode - 2) * 64;

  __m128i vcoeff0 = _mm_load_si128((const __m128i*) &intra_chroma_linear_interpolation_weights_w32_hor[weigth_offset + 0]);
  __m128i vcoeff1 = _mm_load_si128((const __m128i*) &intra_chroma_linear_interpolation_weights_w32_hor[weigth_offset + 16]);
  __m128i vcoeff2 = _mm_load_si128((const __m128i*) &intra_chroma_linear_interpolation_weights_w32_hor[weigth_offset + 32]);
  __m128i vcoeff3 = _mm_load_si128((const __m128i*) &intra_chroma_linear_interpolation_weights_w32_hor[weigth_offset + 48]);
  __m128i vshuf0 = _mm_load_si128((const __m128i*) &intra_chroma_linear_interpolation_shuffle_vectors_w32_hor[shuf_offset + 0]);
  __m128i vshuf1 = _mm_load_si128((const __m128i*) &intra_chroma_linear_interpolation_shuffle_vectors_w32_hor[shuf_offset + 16]);
  __m128i vshuf2 = _mm_load_si128((const __m128i*) &intra_chroma_linear_interpolation_shuffle_vectors_w32_hor[shuf_offset + 32]);
  __m128i vshuf3 = _mm_load_si128((const __m128i*) &intra_chroma_linear_interpolation_shuffle_vectors_w32_hor[shuf_offset + 48]);

  // Load refs from smallest index onwards, shuffle will handle the rest. The smallest index will be at one of these delta int table indices
  // Due to width, two loads are needed, and therefore two offsets. Cannot use 256-bit loads due to alignment issues.
  const int16_t min_offset0 = 1 + MIN(dint[0], dint[15]);
  const int16_t min_offset1 = 1 + MIN(dint[16], dint[31]);

  // Height has to be at least 2. Due to width, handle 1 line at a time
  for (int y = 0; y < height; ++y) {
    // Prepare sources
    __m128i vsrc_tmp0 = _mm_loadu_si128((__m128i*) &ref[min_offset0 + y]);
    __m128i vsrc_tmp1 = _mm_loadu_si128((__m128i*) &ref[min_offset1 + y]);
    __m128i vsrc0 = _mm_shuffle_epi8(vsrc_tmp0, vshuf0);
    __m128i vsrc1 = _mm_shuffle_epi8(vsrc_tmp0, vshuf1);
    __m128i vsrc2 = _mm_shuffle_epi8(vsrc_tmp1, vshuf2);
    __m128i vsrc3 = _mm_shuffle_epi8(vsrc_tmp1, vshuf3);

    __m128i res0 = _mm_maddubs_epi16(vsrc0, vcoeff0);
    __m128i res1 = _mm_maddubs_epi16(vsrc1, vcoeff1);
    __m128i res2 = _mm_maddubs_epi16(vsrc2, vcoeff2);
    __m128i res3 = _mm_maddubs_epi16(vsrc3, vcoeff3);
    res0 = _mm_add_epi16(res0, v16s);
    res1 = _mm_add_epi16(res1, v16s);
    res2 = _mm_add_epi16(res2, v16s);
    res3 = _mm_add_epi16(res3, v16s);
    res0 = _mm_srai_epi16(res0, 5);
    res1 = _mm_srai_epi16(res1, 5);
    res2 = _mm_srai_epi16(res2, 5);
    res3 = _mm_srai_epi16(res3, 5);

    _mm_store_si128((__m128i*)&dst[0], _mm_packus_epi16(res0, res1));
    _mm_store_si128((__m128i*)&dst[16], _mm_packus_epi16(res2, res3));
    dst += 32;
  }
}


static void angular_pred_linear_filter_w8_ver_wide_angle_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int mode, const int16_t* delta_int, const int16_t* delta_fract)
{
  const int width = 8;
  const int16_t* dint = delta_int;
  const __m128i v16s = _mm_set1_epi16(16);
  // Height has to be at least 2, handle 2 lines at once
  for (int y = 0; y < height; y += 2) {
    uvg_pixel src[32];
    int16_t coeff_tmp[2];
    // TODO: get rid of this slow crap, this is just here to test the calculations
    for (int yy = 0; yy < 2; ++yy) {
      for (int x = 0, d = 0; x < width; ++x, d += 2) {
        src[yy * 16 + d + 0] = ref[dint[yy] + 1 + x + 0];
        src[yy * 16 + d + 1] = ref[dint[yy] + 1 + x + 1];
      }
      int8_t tmp[2] = { 32 - delta_fract[y + yy], delta_fract[y + yy] };
      coeff_tmp[yy] = *(int16_t*)tmp;
    }
    dint += 2;

    const __m128i vcoeff0 = _mm_set1_epi16(coeff_tmp[0]);
    const __m128i vcoeff1 = _mm_set1_epi16(coeff_tmp[1]);

    const __m128i* vsrc0 = (const __m128i*) & src[0];
    const __m128i* vsrc1 = (const __m128i*) & src[16];

    __m128i res0 = _mm_maddubs_epi16(*vsrc0, vcoeff0);
    __m128i res1 = _mm_maddubs_epi16(*vsrc1, vcoeff1);
    res0 = _mm_add_epi16(res0, v16s);
    res1 = _mm_add_epi16(res1, v16s);
    res0 = _mm_srai_epi16(res0, 5);
    res1 = _mm_srai_epi16(res1, 5);

    _mm_store_si128((__m128i*)dst, _mm_packus_epi16(res0, res1));
    dst += 16;
  }
}


static void angular_pred_linear_filter_w16_ver_wide_angle_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int mode, const int16_t* delta_int, const int16_t* delta_fract)
{
  const int width = 16;
  const int16_t* dint = delta_int;
  const __m128i v16s = _mm_set1_epi16(16);
  // Height has to be at least 2, handle 1 line at a time
  for (int y = 0; y < height; ++y) {
    uvg_pixel src[32];
    // TODO: get rid of this slow crap, this is just here to test the calculations
    for (int x = 0, d = 0; x < width; ++x, d += 2) {
      src[d + 0] = ref[*dint + 1 + x + 0];
      src[d + 1] = ref[*dint + 1 + x + 1];
    }
    dint++;

    int8_t tmp[2] = { 32 - delta_fract[y], delta_fract[y] };
    const int16_t coeff_tmp = *(int16_t*)tmp;
    const __m128i vcoeff = _mm_set1_epi16(coeff_tmp);
    
    const __m128i* vsrc0 = (const __m128i*) & src[0];
    const __m128i* vsrc1 = (const __m128i*) & src[16];

    __m128i res0 = _mm_maddubs_epi16(*vsrc0, vcoeff);
    __m128i res1 = _mm_maddubs_epi16(*vsrc1, vcoeff);
    res0 = _mm_add_epi16(res0, v16s);
    res1 = _mm_add_epi16(res1, v16s);
    res0 = _mm_srai_epi16(res0, 5);
    res1 = _mm_srai_epi16(res1, 5);

    _mm_store_si128((__m128i*)dst, _mm_packus_epi16(res0, res1));
    dst += 16;
  }
}


static void angular_pred_linear_filter_w32_ver_wide_angle_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int mode, const int16_t* delta_int, const int16_t* delta_fract)
{
  const int width = 32;
  const int16_t* dint = delta_int;
  const __m256i v16s = _mm256_set1_epi16(16);
  // Height has to be at least 2, handle 1 line at a time
  for (int y = 0; y < height; ++y) {
    uvg_pixel src[64];
    // TODO: get rid of this slow crap, this is just here to test the calculations
    for (int x = 0, d = 0; x < width; ++x, d += 2) {
      src[d + 0] = ref[*dint + 1 + x + 0];
      src[d + 1] = ref[*dint + 1 + x + 1];
    }
    dint++;

    int8_t tmp[2] = { 32 - delta_fract[y], delta_fract[y] };
    const int16_t coeff_tmp = *(int16_t*)tmp;
    const __m256i vcoeff = _mm256_set1_epi16(coeff_tmp);

    const __m256i* vsrc0 = (const __m256i*) & src[0];
    const __m256i* vsrc1 = (const __m256i*) & src[32];

    __m256i res0 = _mm256_maddubs_epi16(*vsrc0, vcoeff);
    __m256i res1 = _mm256_maddubs_epi16(*vsrc1, vcoeff);
    res0 = _mm256_add_epi16(res0, v16s);
    res1 = _mm256_add_epi16(res1, v16s);
    res0 = _mm256_srai_epi16(res0, 5);
    res1 = _mm256_srai_epi16(res1, 5);

    __m256i vfinal = _mm256_packus_epi16(res0, res1);
    vfinal = _mm256_permute4x64_epi64(vfinal, _MM_SHUFFLE(3, 1, 2, 0));

    _mm256_store_si256((__m256i*)dst, vfinal);
    dst += 32;
  }
}


static void angular_pred_linear_filter_w4_hor_wide_angle_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int mode, const int16_t* delta_int)
{
  const __m128i v16s = _mm_set1_epi16(16);

  const int mode_idx = mode < 2 ? mode + 12 : 80 - mode;
  const int table_offset = mode_idx * 128;

  const __m128i vcoeff0 = _mm_load_si128((const __m128i*) &intra_chroma_linear_interpolation_weights_w16_hor_wide_angle[table_offset + 0]);
  const __m128i vcoeff1 = _mm_load_si128((const __m128i*) &intra_chroma_linear_interpolation_weights_w16_hor_wide_angle[table_offset + 16]);

  const __m128i vshuf = _mm_setr_epi8(
    0x00, 0x01, 0x01, 0x02, 0x02, 0x03, 0x03, 0x04,
    0x08, 0x09, 0x09, 0x0a, 0x0a, 0x0b, 0x0b, 0x0c
  );

  const __m128i vtranspose = _mm_setr_epi8(
    0x00, 0x04, 0x08, 0x0c, 0x01, 0x05, 0x09, 0x0d,
    0x02, 0x06, 0x0a, 0x0e, 0x03, 0x07, 0x0b, 0x0f
  );

  const __m256i vidx = _mm256_setr_epi64x(delta_int[0], delta_int[1], delta_int[2], delta_int[3]);

  // Height has to be at least 4, handle 4 lines at once
  for (int y = 0; y < height; y += 4) {
    const __m256i vsrc_raw = _mm256_i64gather_epi64((const long long*)&ref[y + 1], vidx, 1);

    __m128i vsrc0 = _mm256_extracti128_si256(vsrc_raw, 0);
    __m128i vsrc1 = _mm256_extracti128_si256(vsrc_raw, 1);
    
    vsrc0 = _mm_shuffle_epi8(vsrc0, vshuf);
    vsrc1 = _mm_shuffle_epi8(vsrc1, vshuf);

    __m128i res0 = _mm_maddubs_epi16(vsrc0, vcoeff0);
    __m128i res1 = _mm_maddubs_epi16(vsrc1, vcoeff1);
    res0 = _mm_add_epi16(res0, v16s);
    res1 = _mm_add_epi16(res1, v16s);
    res0 = _mm_srai_epi16(res0, 5);
    res1 = _mm_srai_epi16(res1, 5);
    __m128i vfinal = _mm_packus_epi16(res0, res1);
    vfinal = _mm_shuffle_epi8(vfinal, vtranspose);

    _mm_store_si128((__m128i*)dst, vfinal);
    dst += 16;
  }
}


static void angular_pred_linear_filter_w8_hor_wide_angle_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int mode, const int16_t* delta_int, const int16_t* delta_fract)
{
  const int width = 8;
  const __m128i v16s = _mm_set1_epi16(16);

  const __m128i vshuf = _mm_setr_epi8(
    0x00, 0x01, 0x01, 0x02, 0x02, 0x03, 0x03, 0x04,
    0x08, 0x09, 0x09, 0x0a, 0x0a, 0x0b, 0x0b, 0x0c
  );

  const __m128i vtranspose = _mm_setr_epi8(
    0x00, 0x04, 0x08, 0x0c, 0x01, 0x05, 0x09, 0x0d,
    0x02, 0x06, 0x0a, 0x0e, 0x03, 0x07, 0x0b, 0x0f
  );

  const __m128i vidxshuf = _mm_setr_epi8(0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
                                         0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff); // Don't care
  __m128i vidx_raw = _mm_load_si128((__m128i*)delta_int);

  const __m256i vidx0 = _mm256_cvtepi16_epi64(vidx_raw);
  vidx_raw = _mm_shuffle_epi8(vidx_raw, vidxshuf);
  const __m256i vidx1 = _mm256_cvtepi16_epi64(vidx_raw);

  const int mode_idx = mode < 2 ? mode + 12 : 80 - mode;
  const int table_offset = mode_idx * 128;

  const __m128i vcoeff0 = _mm_load_si128((const __m128i*) & intra_chroma_linear_interpolation_weights_w16_hor_wide_angle[table_offset + 0]);
  const __m128i vcoeff1 = _mm_load_si128((const __m128i*) & intra_chroma_linear_interpolation_weights_w16_hor_wide_angle[table_offset + 16]);
  const __m128i vcoeff2 = _mm_load_si128((const __m128i*) & intra_chroma_linear_interpolation_weights_w16_hor_wide_angle[table_offset + 32]);
  const __m128i vcoeff3 = _mm_load_si128((const __m128i*) & intra_chroma_linear_interpolation_weights_w16_hor_wide_angle[table_offset + 48]);

  // Height has to be at least 2. Handle as 4x4 blocks. Special handling needed when height == 2.
  // TODO: make sure this function is not called when height is 2.
  for (int y = 0; y < height; y += 4) {
    const __m256i vsrc_raw0 = _mm256_i64gather_epi64((const long long*)&ref[y + 1], vidx0, 1);
    const __m256i vsrc_raw1 = _mm256_i64gather_epi64((const long long*)&ref[y + 1], vidx1, 1);

    __m128i vsrc0 = _mm256_extracti128_si256(vsrc_raw0, 0);
    __m128i vsrc1 = _mm256_extracti128_si256(vsrc_raw0, 1);
    __m128i vsrc2 = _mm256_extracti128_si256(vsrc_raw1, 0);
    __m128i vsrc3 = _mm256_extracti128_si256(vsrc_raw1, 1);

    vsrc0 = _mm_shuffle_epi8(vsrc0, vshuf);
    vsrc1 = _mm_shuffle_epi8(vsrc1, vshuf);
    vsrc2 = _mm_shuffle_epi8(vsrc2, vshuf);
    vsrc3 = _mm_shuffle_epi8(vsrc3, vshuf);

    __m128i res0 = _mm_maddubs_epi16(vsrc0, vcoeff0);
    __m128i res1 = _mm_maddubs_epi16(vsrc1, vcoeff1);
    __m128i res2 = _mm_maddubs_epi16(vsrc2, vcoeff2);
    __m128i res3 = _mm_maddubs_epi16(vsrc3, vcoeff3);
    res0 = _mm_add_epi16(res0, v16s);
    res1 = _mm_add_epi16(res1, v16s);
    res2 = _mm_add_epi16(res2, v16s);
    res3 = _mm_add_epi16(res3, v16s);
    res0 = _mm_srai_epi16(res0, 5);
    res1 = _mm_srai_epi16(res1, 5);
    res2 = _mm_srai_epi16(res2, 5);
    res3 = _mm_srai_epi16(res3, 5);

    __m128i vtmp0 = _mm_packus_epi16(res0, res1);
    __m128i vtmp1 = _mm_packus_epi16(res2, res3);
    vtmp0 = _mm_shuffle_epi8(vtmp0, vtranspose);
    vtmp1 = _mm_shuffle_epi8(vtmp1, vtranspose);

    __m128i vfinal0 = _mm_unpacklo_epi32(vtmp0, vtmp1);
    __m128i vfinal1 = _mm_unpackhi_epi32(vtmp0, vtmp1);
    

    _mm_store_si128((__m128i*)&dst[0],  vfinal0);
    _mm_store_si128((__m128i*)&dst[16], vfinal1);
    dst += 32;
  }
}


static void angular_pred_linear_filter_w16_hor_wide_angle_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int mode, const int16_t* delta_int, const int16_t* delta_fract)
{
  const int width = 16;
  const __m128i v16s = _mm_set1_epi16(16);

  const __m128i vshuf = _mm_setr_epi8(
    0x00, 0x01, 0x01, 0x02, 0x02, 0x03, 0x03, 0x04,
    0x08, 0x09, 0x09, 0x0a, 0x0a, 0x0b, 0x0b, 0x0c
  );

  const __m128i vtranspose = _mm_setr_epi8(
    0x00, 0x04, 0x08, 0x0c, 0x01, 0x05, 0x09, 0x0d,
    0x02, 0x06, 0x0a, 0x0e, 0x03, 0x07, 0x0b, 0x0f
  );

  const __m128i vidxshuf = _mm_setr_epi8(0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
                                         0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff); // Don't care

  __m128i vidx_raw0 = _mm_load_si128((__m128i*) & delta_int[0]);
  __m128i vidx_raw1 = _mm_load_si128((__m128i*) & delta_int[8]);

  __m256i vidx[4];
  vidx[0] = _mm256_cvtepi16_epi64(vidx_raw0);
  vidx_raw0 = _mm_shuffle_epi8(vidx_raw0, vidxshuf);
  vidx[1] = _mm256_cvtepi16_epi64(vidx_raw0);

  vidx[2] = _mm256_cvtepi16_epi64(vidx_raw1);
  vidx_raw1 = _mm_shuffle_epi8(vidx_raw1, vidxshuf);
  vidx[3] = _mm256_cvtepi16_epi64(vidx_raw1);

  const int mode_idx = mode < 2 ? mode + 12 : 80 - mode;
  const int table_offset = mode_idx * 128;

  __m128i vcoeff[8];
  for (int i = 0, o = 0; i < 8; ++i, o += 16) {
    vcoeff[i] = _mm_load_si128((const __m128i*) & intra_chroma_linear_interpolation_weights_w16_hor_wide_angle[table_offset + o]);
  }

  // Height has to be at least 2. Handle as 4x4 blocks. Special handling needed when height < 4.
  // TODO: make sure this function is not called when height is less than 4.
  for (int y = 0; y < height; y += 4) {
    __m128i vtmp[4];
    for (int x = 0, v = 0, c = 0; x < width; x += 8, v += 2, c += 4) {
      const __m256i vsrc_raw0 = _mm256_i64gather_epi64((const long long*)&ref[y + 1], vidx[v + 0], 1);
      const __m256i vsrc_raw1 = _mm256_i64gather_epi64((const long long*)&ref[y + 1], vidx[v + 1], 1);

      __m128i vsrc0 = _mm256_extracti128_si256(vsrc_raw0, 0);
      __m128i vsrc1 = _mm256_extracti128_si256(vsrc_raw0, 1);
      __m128i vsrc2 = _mm256_extracti128_si256(vsrc_raw1, 0);
      __m128i vsrc3 = _mm256_extracti128_si256(vsrc_raw1, 1);

      vsrc0 = _mm_shuffle_epi8(vsrc0, vshuf);
      vsrc1 = _mm_shuffle_epi8(vsrc1, vshuf);
      vsrc2 = _mm_shuffle_epi8(vsrc2, vshuf);
      vsrc3 = _mm_shuffle_epi8(vsrc3, vshuf);

      __m128i res0 = _mm_maddubs_epi16(vsrc0, vcoeff[c + 0]);
      __m128i res1 = _mm_maddubs_epi16(vsrc1, vcoeff[c + 1]);
      __m128i res2 = _mm_maddubs_epi16(vsrc2, vcoeff[c + 2]);
      __m128i res3 = _mm_maddubs_epi16(vsrc3, vcoeff[c + 3]);
      res0 = _mm_add_epi16(res0, v16s);
      res1 = _mm_add_epi16(res1, v16s);
      res2 = _mm_add_epi16(res2, v16s);
      res3 = _mm_add_epi16(res3, v16s);
      res0 = _mm_srai_epi16(res0, 5);
      res1 = _mm_srai_epi16(res1, 5);
      res2 = _mm_srai_epi16(res2, 5);
      res3 = _mm_srai_epi16(res3, 5);

      vtmp[v + 0] = _mm_packus_epi16(res0, res1);
      vtmp[v + 1] = _mm_packus_epi16(res2, res3);
    }
    vtmp[0] = _mm_shuffle_epi8(vtmp[0], vtranspose);
    vtmp[1] = _mm_shuffle_epi8(vtmp[1], vtranspose);
    vtmp[2] = _mm_shuffle_epi8(vtmp[2], vtranspose);
    vtmp[3] = _mm_shuffle_epi8(vtmp[3], vtranspose);

    __m128i vupk32_lo0 = _mm_unpacklo_epi32(vtmp[0], vtmp[1]);
    __m128i vupk32_hi0 = _mm_unpackhi_epi32(vtmp[0], vtmp[1]);
    __m128i vupk32_lo1 = _mm_unpacklo_epi32(vtmp[2], vtmp[3]);
    __m128i vupk32_hi1 = _mm_unpackhi_epi32(vtmp[2], vtmp[3]);

    __m128i vfinal0 = _mm_unpacklo_epi64(vupk32_lo0, vupk32_lo1);
    __m128i vfinal1 = _mm_unpackhi_epi64(vupk32_lo0, vupk32_lo1);
    __m128i vfinal2 = _mm_unpacklo_epi64(vupk32_hi0, vupk32_hi1);
    __m128i vfinal3 = _mm_unpackhi_epi64(vupk32_hi0, vupk32_hi1);

    _mm_store_si128((__m128i*) & dst[0], vfinal0);
    _mm_store_si128((__m128i*) & dst[16], vfinal1);
    _mm_store_si128((__m128i*) & dst[32], vfinal2);
    _mm_store_si128((__m128i*) & dst[48], vfinal3);
    dst += 64;
  }
}


static void angular_pred_linear_filter_hor_avx2(uvg_pixel* dst, uvg_pixel* ref, const int width, const int height, const int16_t* delta_int, const int16_t* delta_fract)
{
  // 2-tap linear filter

  // Handle filtering in 4x4 blocks
  for (int y = 0; y < height; y += 4) {
    const __m256i vref = _mm256_loadu_si256((const __m256i*) & ref[y + 1]);
    for (int x = 0; x < width; x += 4) {

    }
  }
}


static void angular_pred_non_fractional_angle_pxl_copy_ver_avx2(uvg_pixel* dst, uvg_pixel* ref, const int width, const int height, const int16_t* delta_int)
{
  // Note: this probably won't work for wide angle modes.
  for (int y = 0; y < height; ++y) {
    uvg_pixel* dst_row = dst + y * width;
    uvg_pixel* ref_row = ref + delta_int[y] + 1;
    switch (width) {
      case 4: memcpy(dst_row, ref_row, 4 * sizeof(uvg_pixel)); break;
      case 8: memcpy(dst_row, ref_row, 8 * sizeof(uvg_pixel)); break;
      case 16: memcpy(dst_row, ref_row, 16 * sizeof(uvg_pixel)); break;
      case 32: memcpy(dst_row, ref_row, 32 * sizeof(uvg_pixel)); break;
      case 64: memcpy(dst_row, ref_row, 64 * sizeof(uvg_pixel)); break;
    }
  }
}

static void angular_pred_non_fractional_angle_pxl_copy_hor_avx2(uvg_pixel* dst, uvg_pixel* ref, const int width, const int height, const int16_t* delta_int)
{
  // TODO: replace this generic solution after testing
  for (int y = 0; y < height; ++y) { 
    for (int x = 0; x < width; ++x) {
      dst[y * width + x] = ref[delta_int[x] + y + 1];
    }
  }
}


static void angular_pdpc_ver_w4_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int height, const int scale, const int mode_disp)
{
  const int width = 4;
  int16_t left[4][4];

  int limit = MIN(3 << scale, width);

  __m128i vseq = _mm_setr_epi32(0, 1, 2, 3);
  //__m128i vidx = _mm_slli_epi32(vseq, 2); // 2 is log2 width
  __m256i v32s = _mm256_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int offset = scale * 16;
  const __m256i vweight = _mm256_load_si256((const __m256i*)&intra_pdpc_w4_ver_weight[offset]);

  const int inv_angle_offset = mode_disp * 64;
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  // For a 4 width block, height must be at least 4. Handle 4 lines at once.
  for (int y = 0; y < height; y += 4) {
    for (int xx = 0; xx < width; ++xx) {
      for (int yy = 0; yy < 4; ++yy) {
        left[yy][xx] = ref_side[(y + yy) + shifted_inv_angle_sum[xx] + 1];
      }
    }

    __m128i vdst = _mm_i32gather_epi32((const int32_t*)(dst + y * width), vseq, 4);
    __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);
    __m256i vleft = _mm256_loadu_si256((__m256i*)left);

    __m256i accu = _mm256_sub_epi16(vleft, vdst16);
    accu = _mm256_mullo_epi16(vweight, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    _mm_store_si128((__m128i*)(dst + (y * width)), filtered);
  }
}


static void angular_pdpc_ver_w8_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int height, const int scale, const int mode_disp)
{
  const int width = 8;
  
  int limit = MIN(3 << scale, width);

  __m128i vseq = _mm_setr_epi32(0x00, 0x00, 0x01, 0x00);
  __m128i vidx = _mm_slli_epi64(vseq, 3); // 3 is log2 width
  __m256i v32s = _mm256_set1_epi16(32);

  const int offset = scale * 16;
  const __m256i vweight = _mm256_load_si256((const __m256i*)&intra_pdpc_w8_ver_weight[offset]);

  const int inv_angle_offset = mode_disp * 64;
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  // For width 8, height must be at least 2. Handle 2 lines at once.
  for (int y = 0; y < height; y += 2) {
    ALIGNED(32) int16_t left[16] = {0};
    for (int xx = 0; xx < limit; ++xx) {
      for (int yy = 0; yy < 2; ++yy) {
        left[yy * width +xx] = ref_side[(y + yy) + shifted_inv_angle_sum[xx] + 1];
      }
    }

    __m128i vdst = _mm_i64gather_epi64((const int64_t*)(dst + y * width), vseq, 8);
    __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);
    __m256i vleft = _mm256_loadu_si256((__m256i*)left);
    
    __m256i accu = _mm256_sub_epi16(vleft, vdst16);
    accu = _mm256_mullo_epi16(vweight, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    _mm_store_si128((__m128i*)(dst + (y * width)), filtered);
  }
}


static void angular_pdpc_ver_w16_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  __m256i v32s = _mm256_set1_epi16(32);
  const int scale = 2; // Other functions handle scales 0 and 1
  int limit = 12; // With scale 2, limit is always 12.

  const int offset = scale * 16;
  const __m256i vweight = _mm256_load_si256((const __m256i*)&intra_pdpc_w16_ver_weight[offset]);

  const int inv_angle_offset = mode_disp * 64;
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  for (int y = 0; y < height; ++y) {
    for (int  x = 0; x < limit; x += 16) {
      ALIGNED(32) int16_t left[16] = {0};
      for (int xx = 0; x + xx < limit; ++xx) {
        left[xx] = ref_side[y + shifted_inv_angle_sum[xx] + 1];
      }

      __m128i vdst = _mm_load_si128((const __m128i*)(dst + (y * width + x)));
      __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);
      __m256i vleft = _mm256_loadu_si256((__m256i*)left);

      __m256i accu = _mm256_sub_epi16(vleft, vdst16);
      accu = _mm256_mullo_epi16(vweight, accu);
      accu = _mm256_add_epi16(accu, v32s);
      accu = _mm256_srai_epi16(accu, 6);
      accu = _mm256_add_epi16(vdst16, accu);

      __m128i lo = _mm256_castsi256_si128(accu);
      __m128i hi = _mm256_extracti128_si256(accu, 1);
      __m128i filtered = _mm_packus_epi16(lo, hi);

      _mm_store_si128((__m128i*)(dst + (y * width + x)), filtered);
    }
  }
}

static void angular_pdpc_ver_w16_scale0_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  // NOTE: This function is just the w4 function, retrofitted to work with width 16 and up when scale is 0.
  // Since scale is 0, limit is 3 and therefore there is no meaningful work to be done when x > 3, so only the first column of 4x4 chunks is handled.
  // NOTE: This function also works with width 8 when scale is 0, the name w16 might be a bit misleading.
  const int scale = 0;
  int16_t left[4][4];
  const int log2_width = uvg_g_convert_to_log2[width];

  const int limit = 3;

  __m128i vseq = _mm_setr_epi32(0, 1, 2, 3);
  __m128i vidx = _mm_slli_epi32(vseq, log2_width);
  __m256i v32s = _mm256_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int offset = scale * 16;
  const __m256i vweight = _mm256_load_si256((const __m256i*) & intra_pdpc_w4_ver_weight[offset]);

  const int inv_angle_offset = mode_disp * 64;
  int16_t shifted_inv_angle_sum[64];
  memcpy(shifted_inv_angle_sum, &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset], height * sizeof(int16_t)); // TODO: would this be faster if the max amount (64) would be always loaded?

  // For a 4 width block, height must be at least 4. Handle 4 lines at once.
  for (int y = 0; y < height; y += 4) {
    for (int xx = 0; xx < 4; ++xx) {
      for (int yy = 0; yy < 4; ++yy) {
        left[yy][xx] = ref_side[(y + yy) + shifted_inv_angle_sum[xx] + 1];
      }
    }

    __m128i vdst = _mm_i32gather_epi32((const int32_t*)(dst + y * width), vidx, 1);
    __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);
    __m256i vleft = _mm256_loadu_si256((__m256i*)left);

    __m256i accu = _mm256_sub_epi16(vleft, vdst16);
    accu = _mm256_mullo_epi16(vweight, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    *(uint32_t*)(dst + (y + 0) * width) = _mm_extract_epi32(filtered, 0);
    *(uint32_t*)(dst + (y + 1) * width) = _mm_extract_epi32(filtered, 1);
    *(uint32_t*)(dst + (y + 2) * width) = _mm_extract_epi32(filtered, 2);
    *(uint32_t*)(dst + (y + 3) * width) = _mm_extract_epi32(filtered, 3);
  }
}

static void angular_pdpc_ver_w16_scale1_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  // NOTE: This function is just the w8 function, retrofitted to work with width 16 and up when scale is 1.
  // Since scale is 1, limit is 6 and therefore there is no meaningful work to be done when x > 6, so only the first column of 8x2 chunks is handled.
  const int scale = 1;
  const int log2_width = uvg_g_convert_to_log2[width];

  const int limit = 6;

  __m128i vseq = _mm_set_epi64x(1, 0);
  __m128i vidx = _mm_slli_epi32(vseq, log2_width);
  __m256i v32s = _mm256_set1_epi16(32);

  const int offset = scale * 16;
  const __m256i vweight = _mm256_load_si256((const __m256i*) &intra_pdpc_w8_ver_weight[offset]);

  const int inv_angle_offset = mode_disp * 64;
  const int16_t *shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  // For width 8, height must be at least 2. Handle 2 lines at once.
  for (int y = 0; y < height; y += 2) {
    ALIGNED(32) int16_t left[16] = { 0 };
    for (int yy = 0; yy < 2; ++yy) {
      for (int xx = 0; xx < limit; ++xx) {
        left[yy * 8 + xx] = ref_side[(y + yy) + shifted_inv_angle_sum[xx] + 1];
      }
    }

    __m128i vdst = _mm_i64gather_epi64((const int64_t*)(dst + y * width), vidx, 1);
    __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);
    __m256i vleft = _mm256_loadu_si256((__m256i*)left);

    __m256i accu = _mm256_sub_epi16(vleft, vdst16);
    accu = _mm256_mullo_epi16(vweight, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    *(uint64_t*)(dst + (y + 0) * width) = _mm_extract_epi64(filtered, 0);
    *(uint64_t*)(dst + (y + 1) * width) = _mm_extract_epi64(filtered, 1);
  }
}

static void angular_pdpc_ver_w16_scale2_high_angle_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  __m256i v32s = _mm256_set1_epi16(32);
  const int scale = 2; // Other functions handle scales 0 and 1
  int limit = 12; // With scale 2, limit is always 12.

  const int offset = scale * 16;
  const int inv_angle_offset = mode_disp * 64;
  const int shuf_offset = mode_disp * 16;
  
  const __m256i vweight = _mm256_load_si256((const __m256i*) &intra_pdpc_w16_ver_weight[offset]);
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];
  const __m128i vshuf = _mm_load_si128((const __m128i*) &intra_pdpc_shuffle_vectors_w16_scale2_ver[shuf_offset]);

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < limit; x += 16) {
      /*ALIGNED(32) int16_t left[16] = { 0 };
      for (int xx = 0; x + xx < limit; ++xx) {
        left[xx] = ref_side[y + shifted_inv_angle_sum[xx] + 1];
      }*/
      __m128i vleft = _mm_loadu_si128((__m128i*) & ref_side[y + shifted_inv_angle_sum[0] + 1]);
      vleft = _mm_shuffle_epi8(vleft, vshuf);

      __m128i vdst = _mm_load_si128((const __m128i*)(dst + (y * width + x)));
      __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);
      __m256i vleft16 = _mm256_cvtepu8_epi16(vleft);

      __m256i accu = _mm256_sub_epi16(vleft16, vdst16);
      accu = _mm256_mullo_epi16(vweight, accu);
      accu = _mm256_add_epi16(accu, v32s);
      accu = _mm256_srai_epi16(accu, 6);
      accu = _mm256_add_epi16(vdst16, accu);

      __m128i lo = _mm256_castsi256_si128(accu);
      __m128i hi = _mm256_extracti128_si256(accu, 1);
      __m128i filtered = _mm_packus_epi16(lo, hi);

      _mm_store_si128((__m128i*)(dst + (y * width + x)), filtered);
    }
  }
}


static void angular_pdpc_ver_w4_high_angle_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int height, const int scale, const int mode_disp)
{
  const int width = 4;
  int16_t left[4][4];

  int limit = MIN(3 << scale, width);

  __m128i vseq = _mm_setr_epi32(0, 1, 2, 3);
  //__m128i vidx = _mm_slli_epi32(vseq, 2); // 2 is log2 width
  __m256i v32s = _mm256_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int offset = scale * 16;
  const int inv_angle_offset = mode_disp * 64;
  const int shuf_offset = mode_disp * 16;
  
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];
  const __m256i vweight = _mm256_load_si256((const __m256i*) &intra_pdpc_w4_ver_weight[offset]);
  const __m128i vshuf = _mm_loadu_si128((__m128i*) &intra_pdpc_shuffle_vectors_w4_ver[shuf_offset]);

  // For a 4 width block, height must be at least 4. Handle 4 lines at once.
  for (int y = 0; y < height; y += 4) {
    /*for (int xx = 0; xx < width; ++xx) {
      for (int yy = 0; yy < 4; ++yy) {
        left[yy][xx] = ref_side[(y + yy) + shifted_inv_angle_sum[xx] + 1];
      }
    }*/
    __m128i vleft = _mm_loadu_si128((__m128i*)&ref_side[y + shifted_inv_angle_sum[0] + 1]);
    vleft = _mm_shuffle_epi8(vleft, vshuf);

    __m128i vdst = _mm_i32gather_epi32((const int32_t*)(dst + y * width), vseq, 4);
    __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);
    __m256i vleft16 = _mm256_cvtepu8_epi16(vleft);

    __m256i accu = _mm256_sub_epi16(vleft16, vdst16);
    accu = _mm256_mullo_epi16(vweight, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    _mm_store_si128((__m128i*)(dst + (y * width)), filtered);
  }
}

static void angular_pdpc_ver_4x4_scale0_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  // This function is just the w4 function, retrofitted to work with any width when scale is 0. If width is 4, use a specialized function instead.
  // Since scale is 0, limit is 3 and therefore there is no meaningful work to be done when x > 3, so only the first column of 4x4 chunks is handled.
  const int scale = 0;
  int16_t left[4][4];
  const int log2_width = uvg_g_convert_to_log2[width];

  const int limit = 3;

  __m128i vseq = _mm_setr_epi32(0, 1, 2, 3);
  __m128i vidx = _mm_slli_epi32(vseq, log2_width);
  __m256i v32s = _mm256_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int offset = scale * 16;
  const __m256i vweight = _mm256_load_si256((const __m256i*) &intra_pdpc_w4_ver_weight[offset]);

  const int inv_angle_offset = mode_disp * 64;
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  // For a 4 width block, height must be at least 4. Handle 4 lines at once.
  for (int y = 0; y < height; y += 4) {
    for (int xx = 0; xx < 4; ++xx) {
      for (int yy = 0; yy < 4; ++yy) {
        left[yy][xx] = ref_side[(y + yy) + shifted_inv_angle_sum[xx] + 1];
      }
    }

    __m128i vdst = _mm_i32gather_epi32((const int32_t*)(dst + y * width), vidx, 1);
    __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);
    __m256i vleft = _mm256_loadu_si256((__m256i*)left);

    __m256i accu = _mm256_sub_epi16(vleft, vdst16);
    accu = _mm256_mullo_epi16(vweight, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    *(uint32_t*)(dst + (y + 0) * width) = _mm_extract_epi32(filtered, 0);
    *(uint32_t*)(dst + (y + 1) * width) = _mm_extract_epi32(filtered, 1);
    *(uint32_t*)(dst + (y + 2) * width) = _mm_extract_epi32(filtered, 2);
    *(uint32_t*)(dst + (y + 3) * width) = _mm_extract_epi32(filtered, 3);
  }
}

static void angular_pdpc_ver_4x4_scale0_high_angle_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  // This function is just the w4 function, retrofitted to work with any width when scale is 0. If width is 4, use a specialized function instead.
  // Since scale is 0, limit is 3 and therefore there is no meaningful work to be done when x > 3, so only the first column of 4x4 chunks is handled.
  // This function handles cases where prediction angle is high. For PDPC, this means the needed reference samples are close together, enabling more effective loading.
  const int scale = 0;
  int16_t left[4][4];
  const int log2_width = uvg_g_convert_to_log2[width];

  const int limit = 3;

  __m128i vseq = _mm_setr_epi32(0, 1, 2, 3);
  __m128i vidx = _mm_slli_epi32(vseq, log2_width);
  __m256i v32s = _mm256_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int offset = scale * 16;
  const int inv_angle_offset = mode_disp * 64;
  const int shuf_offset = mode_disp * 16;
  
  const __m256i vweight = _mm256_load_si256((const __m256i*) &intra_pdpc_w4_ver_weight[offset]);
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];
  const __m128i vshuf = _mm_loadu_si128((__m128i*) &intra_pdpc_shuffle_vectors_w4_ver[shuf_offset]);

  // For a 4 width block, height must be at least 4. Handle 4 lines at once.
  for (int y = 0; y < height; y += 4) {
    __m128i vleft = _mm_loadu_si128((__m128i*) & ref_side[y + shifted_inv_angle_sum[0] + 1]);
    vleft = _mm_shuffle_epi8(vleft, vshuf);

    __m128i vdst = _mm_i32gather_epi32((const int32_t*)(dst + y * width), vidx, 1);
    __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);
    __m256i vleft16 = _mm256_cvtepu8_epi16(vleft);

    __m256i accu = _mm256_sub_epi16(vleft16, vdst16);
    accu = _mm256_mullo_epi16(vweight, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    *(uint32_t*)(dst + (y + 0) * width) = _mm_extract_epi32(filtered, 0);
    *(uint32_t*)(dst + (y + 1) * width) = _mm_extract_epi32(filtered, 1);
    *(uint32_t*)(dst + (y + 2) * width) = _mm_extract_epi32(filtered, 2);
    *(uint32_t*)(dst + (y + 3) * width) = _mm_extract_epi32(filtered, 3);
  }
}


static void angular_pdpc_ver_8x2_scale1_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  // NOTE: This function is just the w8 function, retrofitted to work with width 16 and up when scale is 1.
  // Since scale is 1, limit is 6 and therefore there is no meaningful work to be done when x > 6, so only the first column of 8x2 chunks is handled.
  const int scale = 1;
  const int log2_width = uvg_g_convert_to_log2[width];

  const int limit = 6;

  __m128i vseq = _mm_set_epi64x(1, 0);
  __m128i vidx = _mm_slli_epi32(vseq, log2_width);
  __m256i v32s = _mm256_set1_epi16(32);

  const int offset = scale * 16;
  const int inv_angle_offset = mode_disp * 64;
  const int shuf_offset = mode_disp * 16;

  const __m256i vweight = _mm256_load_si256((const __m256i*) & intra_pdpc_w8_ver_weight[offset]);
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];
  const __m128i vshuf = _mm_loadu_si128((__m128i*) & intra_pdpc_shuffle_vectors_8x2_scale1_ver[shuf_offset]);

  // For width 8, height must be at least 2. Handle 2 lines at once.
  for (int y = 0; y < height; y += 2) {
    /*ALIGNED(32) int16_t left[16] = { 0 };
    for (int yy = 0; yy < 2; ++yy) {
      for (int xx = 0; xx < limit; ++xx) {
        left[yy * 8 + xx] = ref_side[(y + yy) + shifted_inv_angle_sum[xx] + 1];
      }
    }*/
    __m128i vleft = _mm_loadu_si128((__m128i*) & ref_side[y + shifted_inv_angle_sum[0] + 1]);
    vleft = _mm_shuffle_epi8(vleft, vshuf);

    __m128i vdst = _mm_i64gather_epi64((const int64_t*)(dst + y * width), vidx, 1);
    __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);
    __m256i vleft16 = _mm256_cvtepu8_epi16(vleft);
    //__m256i vleft = _mm256_loadu_si256((__m256i*)left);

    __m256i accu = _mm256_sub_epi16(vleft16, vdst16);
    accu = _mm256_mullo_epi16(vweight, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    // TODO: if this if branch is deemed to cause slow down, make another version of this, where this check is not needed.
    // If this does not slow down significantly, make this same check in other functions to reduce the function call switch case complexity
    if (width == 8) { 
      _mm_store_si128((__m128i*)(dst + (y * width)), filtered);
    }
    else {
      *(uint64_t*)(dst + (y + 0) * width) = _mm_extract_epi64(filtered, 0);
      *(uint64_t*)(dst + (y + 1) * width) = _mm_extract_epi64(filtered, 1);
    }
  }
}

static void angular_pdpc_ver_8x2_scale2_high_angle_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  // NOTE: This function is just the w8 function, retrofitted to work with width 16 and up when scale is 1.
  // Since scale is 1, limit is 6 and therefore there is no meaningful work to be done when x > 6, so only the first column of 8x2 chunks is handled.
  // This function handles cases where prediction angle is high. For PDPC, this means the needed reference samples are close together, enabling more effective loading.
  const int scale = 1;
  const int log2_width = uvg_g_convert_to_log2[width];

  const int limit = 6;

  __m128i vseq = _mm_set_epi64x(1, 0);
  __m128i vidx = _mm_slli_epi32(vseq, log2_width);
  __m256i v32s = _mm256_set1_epi16(32);

  const int offset = scale * 16;
  const int inv_angle_offset = mode_disp * 64;
  const int shuf_offset = mode_disp * 16;

  const __m256i vweight = _mm256_load_si256((const __m256i*) & intra_pdpc_w8_ver_weight[offset]);
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];
  const __m128i vshuf = _mm_loadu_si128((__m128i*) & intra_pdpc_shuffle_vectors_8x2_scale2_ver[shuf_offset]);

  // For width 8, height must be at least 2. Handle 2 lines at once.
  for (int y = 0; y < height; y += 2) {
    /*ALIGNED(32) int16_t left[16] = { 0 };
    for (int yy = 0; yy < 2; ++yy) {
      for (int xx = 0; xx < limit; ++xx) {
        left[yy * 8 + xx] = ref_side[(y + yy) + shifted_inv_angle_sum[xx] + 1];
      }
    }*/
    __m128i vleft = _mm_loadu_si128((__m128i*) & ref_side[y + shifted_inv_angle_sum[0] + 1]);
    vleft = _mm_shuffle_epi8(vleft, vshuf);

    __m128i vdst = _mm_i64gather_epi64((const int64_t*)(dst + y * width), vidx, 1);
    __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);
    __m256i vleft16 = _mm256_cvtepu8_epi16(vleft);
    //__m256i vleft = _mm256_loadu_si256((__m256i*)left);

    __m256i accu = _mm256_sub_epi16(vleft16, vdst16);
    accu = _mm256_mullo_epi16(vweight, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    // TODO: if this if branch is deemed to cause slow down, make another version of this, where this check is not needed.
    // If this does not slow down significantly, make this same check in other functions to reduce the function call switch case complexity
    if (width == 8) {
      _mm_store_si128((__m128i*)(dst + (y * width)), filtered);
    }
    else {
      *(uint64_t*)(dst + (y + 0) * width) = _mm_extract_epi64(filtered, 0);
      *(uint64_t*)(dst + (y + 1) * width) = _mm_extract_epi64(filtered, 1);
    }
  }
}

static void angular_pdpc_ver_8x2_scale1_high_angle_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  // NOTE: This function is just the w8 function, retrofitted to work with width 16 and up when scale is 1.
  // Since scale is 1, limit is 6 and therefore there is no meaningful work to be done when x > 6, so only the first column of 8x2 chunks is handled.
  // This function handles cases where prediction angle is high. For PDPC, this means the needed reference samples are close together, enabling more effective loading.
  const int scale = 1;
  const int log2_width = uvg_g_convert_to_log2[width];

  const int limit = 6;

  __m128i vseq = _mm_set_epi64x(1, 0);
  __m128i vidx = _mm_slli_epi32(vseq, log2_width);
  __m256i v32s = _mm256_set1_epi16(32);

  const int offset = scale * 16;
  const __m256i vweight = _mm256_load_si256((const __m256i*) & intra_pdpc_w8_ver_weight[offset]);

  const int inv_angle_offset = mode_disp * 64;
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  // For width 8, height must be at least 2. Handle 2 lines at once.
  for (int y = 0; y < height; y += 2) {
    ALIGNED(32) int16_t left[16] = { 0 };
    for (int yy = 0; yy < 2; ++yy) {
      for (int xx = 0; xx < limit; ++xx) {
        left[yy * 8 + xx] = ref_side[(y + yy) + shifted_inv_angle_sum[xx] + 1];
      }
    }

    __m128i vdst = _mm_i64gather_epi64((const int64_t*)(dst + y * width), vidx, 1);
    __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);
    __m256i vleft = _mm256_loadu_si256((__m256i*)left);

    __m256i accu = _mm256_sub_epi16(vleft, vdst16);
    accu = _mm256_mullo_epi16(vweight, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    // TODO: if this if branch is deemed to cause slow down, make another version of this, where this check is not needed.
    // If this does not slow down significantly, make this same check in other functions to reduce the function call switch case complexity
    if (width == 8) {
      _mm_store_si128((__m128i*)(dst + (y * width)), filtered);
    }
    else {
      *(uint64_t*)(dst + (y + 0) * width) = _mm_extract_epi64(filtered, 0);
      *(uint64_t*)(dst + (y + 1) * width) = _mm_extract_epi64(filtered, 1);
    }
  }
}


// Height versions of vertical PDPC, these are unused but left here for archiving purposes. Maybe this method can be refined to be effective.

static void angular_pdpc_ver_h4_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int scale, const int16_t inv_sample_disp)
{
  const int height = 4;

  int limit = MIN(3 << scale, width);
  const int log2_width = uvg_g_convert_to_log2[width];

  const __m256i v32s = _mm256_set1_epi16(32);
  const __m256i wL_shuffle = _mm256_setr_epi8(
    0x00, 0x02, 0x00, 0x02, 0x00, 0x02, 0x00, 0x02,
    0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a,
    0x00, 0x02, 0x00, 0x02, 0x00, 0x02, 0x00, 0x02,
    0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a, 0x08, 0x0a
  );

  for (int x = 0; x < limit; x += 4) {
    int shifted_inv_angle_sum[4] = {0};
    int16_t wL[4] = {0};
    ALIGNED(32) uvg_pixel tmp[16];
    for (int xx = 0; xx < 4; ++xx) {
      shifted_inv_angle_sum[xx] = (256 + (x + xx + 1) * inv_sample_disp) >> 9;
      wL[xx] = (x + xx) < limit ? 32 >> ((2 * (x + xx)) >> scale) : 0;

      tmp[xx * 4 + 0]  = ref_side[0 + shifted_inv_angle_sum[xx] + 1];
      tmp[xx * 4 + 1]  = ref_side[1 + shifted_inv_angle_sum[xx] + 1];
      tmp[xx * 4 + 2]  = ref_side[2 + shifted_inv_angle_sum[xx] + 1];
      tmp[xx * 4 + 3]  = ref_side[3 + shifted_inv_angle_sum[xx] + 1];

    }

    int16_t tmp_dst[16];
    for (int yy = 0; yy < height; ++yy) {
      tmp_dst[0 + yy]  = dst[yy * width + x + 0];
      tmp_dst[4 + yy]  = dst[yy * width + x + 1];
      tmp_dst[8 + yy]  = dst[yy * width + x + 2];
      tmp_dst[12 + yy] = dst[yy * width + x + 3];
    }
    
    __m256i* vdst16 = (__m256i*)tmp_dst;
    __m128i vleft = _mm_load_si128((__m128i*)tmp);
    __m256i vleft16 = _mm256_cvtepu8_epi16(vleft);
    __m256i accu = _mm256_sub_epi16(vleft16, *vdst16);
    __m256i vwL = _mm256_setr_epi64x(wL[0], wL[1], wL[2], wL[3]);
    vwL = _mm256_shuffle_epi8(vwL, wL_shuffle);
    accu = _mm256_mullo_epi16(vwL, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(*vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    const uvg_pixel* result = (uvg_pixel*)&filtered;

    for (int yy = 0; yy < height; ++yy) {
      dst[yy * width + x + 0] = result[0 + yy];
      dst[yy * width + x + 1] = result[4 + yy];
      dst[yy * width + x + 2] = result[8 + yy];
      dst[yy * width + x + 3] = result[12 + yy];
    }
  }
}

static void angular_pdpc_ver_h8_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int scale, const int16_t inv_sample_disp)
{
  const int height = 8;

  int limit = MIN(3 << scale, width);
  __m256i v32s = _mm256_set1_epi16(32);

  for (int x = 0; x < limit; x += 2) {
    int shifted_inv_angle_sum0 = (256 + (x + 0 + 1) * inv_sample_disp) >> 9;
    int shifted_inv_angle_sum1 = (256 + (x + 1 + 1) * inv_sample_disp) >> 9;
    __m128i vwL[2];
    const int16_t wL0 = 32 >> ((2 * (x + 0)) >> scale);
    const int16_t wL1 = (x + 1) < limit ? 32 >> ((2 * (x + 1)) >> scale) : 0;
    vwL[0] = _mm_set1_epi16(wL0);
    vwL[1] = _mm_set1_epi16(wL1);

    ALIGNED(32) int16_t tmp_dst[16];
    for (int yy = 0; yy < height; ++yy) {
      tmp_dst[0 + yy] = dst[(yy) * width + x + 0];
      tmp_dst[8 + yy] = dst[(yy) * width + x + 1];
    }

    ALIGNED(32) uvg_pixel left[16];
    memcpy(&left[0], &ref_side[shifted_inv_angle_sum0 + 1], 8 * sizeof(uvg_pixel));
    memcpy(&left[8], &ref_side[shifted_inv_angle_sum1 + 1], 8 * sizeof(uvg_pixel));

    __m256i vdst16 = _mm256_load_si256((__m256i*)tmp_dst);
    __m128i vleft = _mm_load_si128((__m128i*)left);
    __m256i vleft16 = _mm256_cvtepu8_epi16(vleft);
    __m256i* vwL256 = (__m256i*)vwL;

    __m256i accu = _mm256_sub_epi16(vleft16, vdst16);
    accu = _mm256_mullo_epi16(*vwL256, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    const uvg_pixel* result = (uvg_pixel*)&filtered;
    for (int yy = 0; yy < height; ++yy) {
      dst[(yy) * width + x + 0] = result[0 + yy];
      dst[(yy) * width + x + 1] = result[8 + yy];
    }
    
  }
}

static void angular_pdpc_ver_h16_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int scale, const int16_t inv_sample_disp)
{
  int limit = MIN(3 << scale, width);
  __m256i v32s = _mm256_set1_epi16(32);

  for (int x = 0; x < limit; ++x) {
    int shifted_inv_angle_sum = (256 + (x + 1) * inv_sample_disp) >> 9;
    const int16_t wL = 32 >> ((2 * x) >> scale);
    const __m256i vwL = _mm256_set1_epi16(wL);

    for (int y = 0; y < height; y += 16) {
      ALIGNED(32) int16_t tmp_dst[16];
      for (int yy = 0; yy < 16; ++yy) {
        tmp_dst[yy] = dst[(y + yy) * width + x];
      }
      __m256i vdst16 = _mm256_load_si256((__m256i*)tmp_dst);
      __m128i vleft = _mm_loadu_si128((__m128i*)&ref_side[y + shifted_inv_angle_sum + 1]);
      __m256i vleft16 = _mm256_cvtepu8_epi16(vleft);

      __m256i accu = _mm256_sub_epi16(vleft16, vdst16);
      accu = _mm256_mullo_epi16(vwL, accu);
      accu = _mm256_add_epi16(accu, v32s);
      accu = _mm256_srai_epi16(accu, 6);
      accu = _mm256_add_epi16(vdst16, accu);

      __m128i lo = _mm256_castsi256_si128(accu);
      __m128i hi = _mm256_extracti128_si256(accu, 1);
      __m128i filtered = _mm_packus_epi16(lo, hi);

      const uvg_pixel* result = (uvg_pixel*)&filtered;
      for (int yy = 0; yy < 16; ++yy) {
        dst[(y + yy) * width + x] = result[yy];
      }
    }
  }
}

static void angular_pdpc_hor_w4_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int height, const int scale, const int mode_disp)
{
  const int width = 4;

  int16_t wT[4];
  int8_t ref_top[4][4];

  int limit = MIN(3 << scale, height);

  __m128i vseq = _mm_setr_epi32(0, 1, 2, 3);
  __m128i vidx = _mm_slli_epi32(vseq, 2); // 2 is log2_width
  __m256i v32s = _mm256_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int table_offset = scale * 64;
  const int inv_angle_offset = mode_disp * 64;
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  for (int y = 0, o = 0; y < limit; y += 4, o += 16) {
    for (int yy = 0; yy < 4; ++yy) {
      memcpy(ref_top[yy], &ref_side[shifted_inv_angle_sum[y + yy] + 1], 4 * sizeof(int8_t));
    }
    const int offset = table_offset + o;

    __m128i vpred = _mm_i32gather_epi32((const int32_t*)(dst + y * width), vidx, 1);
    __m256i vpred16 = _mm256_cvtepu8_epi16(vpred);
    __m128i vtop = _mm_loadu_si128((__m128i*)ref_top);
    __m256i vtop16 = _mm256_cvtepu8_epi16(vtop);
    __m256i vwT = _mm256_load_si256((const __m256i*)&intra_pdpc_w4_hor_weight[offset]);

    __m256i accu = _mm256_sub_epi16(vtop16, vpred16);
    accu = _mm256_mullo_epi16(vwT, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vpred16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    _mm_storeu_si128((__m128i*)(dst + (y * width)), filtered);
  }
}

static void angular_pdpc_hor_w8_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int height, const int scale, const int mode_disp)
{
  const int width = 8;

  int limit = MIN(3 << scale, height);

  __m128i vseq = _mm_setr_epi32(0x00, 0x00, 0x01, 0x00);
  __m128i vidx = _mm_slli_epi64(vseq, 3); // 3 is log2 width
  __m256i v32s = _mm256_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int table_offset = scale * 128;
  const int inv_angle_offset = mode_disp * 64;
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  for (int y = 0, o = table_offset; y < limit; y += 2, o += 16) {
    const __m256i vwT = _mm256_load_si256((const __m256i*)&intra_pdpc_w8_hor_weight[o]);
    
    ALIGNED(32) uvg_pixel tmp[16];
    memcpy(&tmp[0], &ref_side[shifted_inv_angle_sum[y + 0] + 1], 8 * sizeof(uvg_pixel));
    memcpy(&tmp[8], &ref_side[shifted_inv_angle_sum[y + 1] + 1], 8 * sizeof(uvg_pixel));
    
    __m128i vpred = _mm_i64gather_epi64((const int64_t*)(dst + y * width), vidx, 1);
    __m256i vpred16 = _mm256_cvtepu8_epi16(vpred);
    __m128i vtop = _mm_load_si128((__m128i*)tmp);
    __m256i vtop16 = _mm256_cvtepu8_epi16(vtop);
    
    __m256i accu = _mm256_sub_epi16(vtop16, vpred16);
    accu = _mm256_mullo_epi16(vwT, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vpred16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    _mm_storeu_si128((__m128i*)(dst + (y * width)), filtered);
  }
}

static void angular_pdpc_hor_w16_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int scale, const int mode_disp)
{
  int limit = MIN(3 << scale, height);
  __m256i v32s = _mm256_set1_epi16(32);

  const int inv_angle_offset = mode_disp * 64;
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  // Handle one line at a time. Skip line if vertical limit reached.
  for (int y = 0; y < limit; ++y) {
    const int16_t wT = 32 >> (2 * (y + 0) >> scale);
    __m256i vwT = _mm256_set1_epi16(wT);

    for (int x = 0; x < width; x += 16) {
      __m128i vpred = _mm_load_si128((__m128i*)(dst + (y * width + x)));
      __m256i vpred16 = _mm256_cvtepu8_epi16(vpred);
      __m128i vtop = _mm_loadu_si128((__m128i*)&ref_side[x + shifted_inv_angle_sum[y] + 1]);
      __m256i vtop16 = _mm256_cvtepu8_epi16(vtop);

      __m256i accu = _mm256_sub_epi16(vtop16, vpred16);
      accu = _mm256_mullo_epi16(vwT, accu);
      accu = _mm256_add_epi16(accu, v32s);
      accu = _mm256_srai_epi16(accu, 6);
      accu = _mm256_add_epi16(vpred16, accu);

      __m128i lo = _mm256_castsi256_si128(accu);
      __m128i hi = _mm256_extracti128_si256(accu, 1);
      __m128i filtered = _mm_packus_epi16(lo, hi);

      _mm_storeu_si128((__m128i*)(dst + (y * width + x)), filtered);
    }
  }
}


static void angular_pdpc_hor_w4_high_angle_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int height, const int scale, const int mode_disp)
{
  const int width = 4;

  int16_t wT[4];
  int8_t ref_top[4][4];

  int limit = MIN(3 << scale, height);

  __m128i vseq = _mm_setr_epi32(0, 1, 2, 3);
  __m128i vidx = _mm_slli_epi32(vseq, 2); // 2 is log2_width
  __m256i v32s = _mm256_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int table_offset = scale * 64;
  const int shuf_offset = mode_disp * 256;
  const int inv_angle_offset = mode_disp * 64;
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  for (int y = 0, o = 0; y < limit; y += 4, o += 16) {
    const __m128i vshuf = _mm_loadu_si128((__m128i*)&intra_pdpc_shuffle_vectors_w4_hor[shuf_offset + o]);
    /*for (int yy = 0; yy < 4; ++yy) {
      memcpy(ref_top[yy], &ref_side[shifted_inv_angle_sum[y + yy] + 1], 4 * sizeof(int8_t));
    }*/

    __m128i vtop = _mm_loadu_si128((__m128i*)&ref_side[shifted_inv_angle_sum[y] + 1]);
    vtop = _mm_shuffle_epi8(vtop, vshuf);

    const int offset = table_offset + o;

    __m128i vpred = _mm_i32gather_epi32((const int32_t*)(dst + y * width), vidx, 1);
    __m256i vpred16 = _mm256_cvtepu8_epi16(vpred);
    //__m128i vtop = _mm_loadu_si128((__m128i*)ref_top);
    __m256i vtop16 = _mm256_cvtepu8_epi16(vtop);
    __m256i vwT = _mm256_load_si256((const __m256i*) & intra_pdpc_w4_hor_weight[offset]);

    __m256i accu = _mm256_sub_epi16(vtop16, vpred16);
    accu = _mm256_mullo_epi16(vwT, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vpred16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    _mm_storeu_si128((__m128i*)(dst + (y * width)), filtered);
  }
}


// This is the non-vectorized version of pdpc mode 18. It is left here for archiving purposes.
static void angular_pdpc_mode18_avx2(uvg_pixel* dst, const uvg_pixel top_left, const uvg_pixel* ref_side, const int width, const int height, const int scale)
{
  const int limit = MIN(3 << scale, height);
  for (int_fast32_t x = 0; x < width; ++x) {
    const uvg_pixel ref_top = ref_side[1 + x];
    for (int yy = 0; yy < limit; ++yy) {
      const int wT = 32 >> ((yy * 2) >> scale);
      const uvg_pixel val = dst[yy * width + x];
      dst[yy * width + x] = CLIP_TO_PIXEL(val + (((ref_top - top_left) * wT + 32) >> 6));
    }
  }
}

static void angular_pdpc_mode18_w4_avx2(uvg_pixel* dst, const uvg_pixel top_left, const uvg_pixel* ref_side, const int height, const int scale)
{
  const int width = 4;
  const int limit = MIN(3 << scale, height);

  __m128i vseq = _mm_setr_epi32(0, 1, 2, 3);
  __m128i vidx = _mm_slli_epi32(vseq, 2); // 2 is log2_width
  __m256i v32s = _mm256_set1_epi16(32);

  const uint32_t ref4 = *(uint32_t*)&ref_side[1];
  
  __m128i vref = _mm_set1_epi32(ref4);
  __m256i vref16 = _mm256_cvtepu8_epi16(vref);

  __m256i vtopleft = _mm256_set1_epi16((uint16_t)top_left);

  // Weight table offset
  const int table_offset = scale * 64;

  for (int y = 0, o = 0; y < limit; y += 4, o += 16) {
    const int offset = table_offset + o;

    __m128i vpred = _mm_i32gather_epi32((const int32_t*)(dst + y * width), vidx, 1);
    __m256i vpred16 = _mm256_cvtepu8_epi16(vpred);
    __m256i vwT = _mm256_load_si256((const __m256i*) &intra_pdpc_w4_hor_weight[offset]);

    __m256i accu = _mm256_sub_epi16(vref16, vtopleft);
    accu = _mm256_mullo_epi16(vwT, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vpred16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    _mm_storeu_si128((__m128i*)(dst + (y * width)), filtered);
  }
}

static void angular_pdpc_mode18_w8_avx2(uvg_pixel* dst, const uvg_pixel top_left, const uvg_pixel* ref_side, const int height, const int scale)
{
  const int width = 8;

  int limit = MIN(3 << scale, height);

  __m128i vseq = _mm_setr_epi32(0x00, 0x00, 0x01, 0x00);
  __m128i vidx = _mm_slli_epi64(vseq, 3); // 3 is log2 width
  __m256i v32s = _mm256_set1_epi16(32);

  const uint64_t ref8 = *(uint64_t*)&ref_side[1];

  __m128i vref = _mm_set1_epi64x(ref8);
  __m256i vref16 = _mm256_cvtepu8_epi16(vref);

  __m256i vtopleft = _mm256_set1_epi16((uint16_t)top_left);

  // Weight table offset
  const int table_offset = scale * 128;

  for (int y = 0, o = table_offset; y < limit; y += 2, o += 16) {
    const __m256i vwT = _mm256_load_si256((const __m256i*) &intra_pdpc_w8_hor_weight[o]);

    __m128i vpred = _mm_i64gather_epi64((const int64_t*)(dst + y * width), vidx, 1);
    __m256i vpred16 = _mm256_cvtepu8_epi16(vpred);
    
    __m256i accu = _mm256_sub_epi16(vref16, vtopleft);
    accu = _mm256_mullo_epi16(vwT, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vpred16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    _mm_storeu_si128((__m128i*)(dst + (y * width)), filtered);
  }
}

static void angular_pdpc_mode18_w16_avx2(uvg_pixel* dst, const uvg_pixel top_left, const uvg_pixel* ref_side, const int height, const int scale)
{
  const int width = 16;
  int limit = MIN(3 << scale, height);
  __m256i v32s = _mm256_set1_epi16(32);

  __m128i vref = _mm_loadu_si128((const __m128i*)&ref_side[1]);
  __m256i vref16 = _mm256_cvtepu8_epi16(vref);

  __m256i vtopleft = _mm256_set1_epi16((uint16_t)top_left);

  // Handle one line at a time. Skip line if vertical limit reached.
  for (int y = 0; y < limit; ++y) {
    const int16_t wT = 32 >> (2 * (y + 0) >> scale);
    __m256i vwT = _mm256_set1_epi16(wT);

    for (int x = 0; x < width; x += 16) {
      __m128i vpred = _mm_load_si128((__m128i*)(dst + (y * width + x)));
      __m256i vpred16 = _mm256_cvtepu8_epi16(vpred);
      
      __m256i accu = _mm256_sub_epi16(vref16, vtopleft);
      accu = _mm256_mullo_epi16(vwT, accu);
      accu = _mm256_add_epi16(accu, v32s);
      accu = _mm256_srai_epi16(accu, 6);
      accu = _mm256_add_epi16(vpred16, accu);

      __m128i lo = _mm256_castsi256_si128(accu);
      __m128i hi = _mm256_extracti128_si256(accu, 1);
      __m128i filtered = _mm_packus_epi16(lo, hi);

      _mm_storeu_si128((__m128i*)(dst + (y * width + x)), filtered);
    }
  }
}

static void angular_pdpc_mode18_w32_avx2(uvg_pixel* dst, const uvg_pixel top_left, const uvg_pixel* ref_side, const int height, const int scale)
{
  const int width = 32;
  int limit = MIN(3 << scale, height);
  __m256i v32s = _mm256_set1_epi16(32);

  __m128i vrefa = _mm_loadu_si128((const __m128i*) &ref_side[1]);
  __m256i vref16a = _mm256_cvtepu8_epi16(vrefa);

  __m128i vrefb = _mm_loadu_si128((const __m128i*) &ref_side[17]);
  __m256i vref16b = _mm256_cvtepu8_epi16(vrefb);

  __m256i vtopleft = _mm256_set1_epi16((uint16_t)top_left);

  // Handle one line at a time. Skip line if vertical limit reached.
  for (int y = 0; y < limit; ++y) {
    const int16_t wT = 32 >> (2 * (y + 0) >> scale);
    __m256i vwT = _mm256_set1_epi16(wT);

    // Calculate first half
    __m128i vpred = _mm_load_si128((__m128i*)(dst + (y * width + 0)));
    __m256i vpred16 = _mm256_cvtepu8_epi16(vpred);

    __m256i accu = _mm256_sub_epi16(vref16a, vtopleft);
    accu = _mm256_mullo_epi16(vwT, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vpred16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    _mm_storeu_si128((__m128i*)(dst + (y * width + 0)), filtered);

    // Calculate second half
    vpred = _mm_load_si128((__m128i*)(dst + (y * width + 16)));
    vpred16 = _mm256_cvtepu8_epi16(vpred);

    accu = _mm256_sub_epi16(vref16b, vtopleft);
    accu = _mm256_mullo_epi16(vwT, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vpred16, accu);

    lo = _mm256_castsi256_si128(accu);
    hi = _mm256_extracti128_si256(accu, 1);
    filtered = _mm_packus_epi16(lo, hi);

    _mm_storeu_si128((__m128i*)(dst + (y * width + 16)), filtered);
  }
}

static void angular_pdpc_mode18_w64_avx2(uvg_pixel* dst, const uvg_pixel top_left, const uvg_pixel* ref_side, const int height, const int scale)
{
  const int width = 64;
  int limit = MIN(3 << scale, height);
  __m256i v32s = _mm256_set1_epi16(32);

  __m128i vrefa = _mm_loadu_si128((const __m128i*) &ref_side[0 + 1]);
  __m256i vref16a = _mm256_cvtepu8_epi16(vrefa);

  __m128i vrefb = _mm_loadu_si128((const __m128i*) &ref_side[16 + 1]);
  __m256i vref16b = _mm256_cvtepu8_epi16(vrefb);

  __m128i vrefc = _mm_loadu_si128((const __m128i*) &ref_side[32 + 1]);
  __m256i vref16c = _mm256_cvtepu8_epi16(vrefc);

  __m128i vrefd = _mm_loadu_si128((const __m128i*) &ref_side[48 + 1]);
  __m256i vref16d = _mm256_cvtepu8_epi16(vrefd);

  __m256i vtopleft = _mm256_set1_epi16((uint16_t)top_left);

  // Handle one line at a time. Skip line if vertical limit reached.
  for (int y = 0; y < limit; ++y) {
    const int16_t wT = 32 >> (2 * (y + 0) >> scale);
    __m256i vwT = _mm256_set1_epi16(wT);

    // Calculate first quarter
    __m128i vpred = _mm_load_si128((__m128i*)(dst + (y * width + 0)));
    __m256i vpred16 = _mm256_cvtepu8_epi16(vpred);

    __m256i accu = _mm256_sub_epi16(vref16a, vtopleft);
    accu = _mm256_mullo_epi16(vwT, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vpred16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    _mm_storeu_si128((__m128i*)(dst + (y * width + 0)), filtered);

    // Calculate second quarter
    vpred = _mm_load_si128((__m128i*)(dst + (y * width + 16)));
    vpred16 = _mm256_cvtepu8_epi16(vpred);

    accu = _mm256_sub_epi16(vref16b, vtopleft);
    accu = _mm256_mullo_epi16(vwT, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vpred16, accu);

    lo = _mm256_castsi256_si128(accu);
    hi = _mm256_extracti128_si256(accu, 1);
    filtered = _mm_packus_epi16(lo, hi);

    _mm_storeu_si128((__m128i*)(dst + (y * width + 16)), filtered);

    // Calculate third quarter
    vpred = _mm_load_si128((__m128i*)(dst + (y * width + 32)));
    vpred16 = _mm256_cvtepu8_epi16(vpred);

    accu = _mm256_sub_epi16(vref16c, vtopleft);
    accu = _mm256_mullo_epi16(vwT, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vpred16, accu);

    lo = _mm256_castsi256_si128(accu);
    hi = _mm256_extracti128_si256(accu, 1);
    filtered = _mm_packus_epi16(lo, hi);

    _mm_storeu_si128((__m128i*)(dst + (y * width + 32)), filtered);

    // Calculate fourth quarter
    vpred = _mm_load_si128((__m128i*)(dst + (y * width + 48)));
    vpred16 = _mm256_cvtepu8_epi16(vpred);

    accu = _mm256_sub_epi16(vref16d, vtopleft);
    accu = _mm256_mullo_epi16(vwT, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vpred16, accu);

    lo = _mm256_castsi256_si128(accu);
    hi = _mm256_extracti128_si256(accu, 1);
    filtered = _mm_packus_epi16(lo, hi);

    _mm_storeu_si128((__m128i*)(dst + (y * width + 48)), filtered);
  }
}


// This is the non-vectorized version of pdpc mode 50. It is left here for archiving purposes.
static void angular_pdpc_mode50_avx2(uvg_pixel* dst, const uvg_pixel top_left, const uvg_pixel* ref_side, const int width, const int height, const int scale)
{
  const int limit = MIN(3 << scale, width);
  for (int y = 0; y < height; ++y) {
    const uvg_pixel left = ref_side[1 + y];
    for (int x = 0; x < limit; x++) {
      const int wL = 32 >> (2 * x >> scale);
      const uvg_pixel val = dst[y * width + x];
      dst[y * width + x] = CLIP_TO_PIXEL(val + ((wL * (left - top_left) + 32) >> 6));
    }
  }
}

static void angular_pdpc_mode50_w4_avx2(uvg_pixel* dst, const uvg_pixel top_left, const uvg_pixel* ref_side, const int height, const int scale)
{
  const int width = 4;
  int limit = MIN(3 << scale, width); // Not used

  //__m128i vseq = _mm_setr_epi32(0, 1, 2, 3);
  //__m128i vidx = _mm_slli_epi32(vseq, 2); // 2 is log2 width
  __m256i v32s = _mm256_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int offset = scale * 16;
  const __m256i vweight = _mm256_load_si256((const __m256i*) &intra_pdpc_w4_ver_weight[offset]);
  const __m256i vtopleft = _mm256_set1_epi16((uint16_t)top_left);

  const __m128i vshuf = _mm_setr_epi8(
    0x00, 0x00, 0x00, 0x00, 0x01, 0x01, 0x01, 0x01,
    0x02, 0x02, 0x02, 0x02, 0x03, 0x03, 0x03, 0x03
  );

  // For a 4 width block, height must be at least 4. Handle 4 lines at once.
  for (int y = 0; y < height; y += 4) {
    const uint32_t ref4 = *(uint32_t*)&ref_side[1 + y];
    __m128i vref = _mm_set1_epi32(ref4);
    vref = _mm_shuffle_epi8(vref, vshuf);
    __m256i vref16 = _mm256_cvtepu8_epi16(vref);

    //__m128i vdst = _mm_i32gather_epi32((const int32_t*)(dst + y * width), vseq, 4);
    __m128i vdst = _mm_load_si128((const __m128i*)(dst + y * width));
    __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);

    __m256i accu = _mm256_sub_epi16(vref16, vtopleft);
    accu = _mm256_mullo_epi16(vweight, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    _mm_store_si128((__m128i*)(dst + (y * width)), filtered);
  }
}

static void angular_pdpc_mode50_w8_avx2(uvg_pixel* dst, const uvg_pixel top_left, const uvg_pixel* ref_side, const int height, const int scale)
{
  const int width = 8;
  int limit = MIN(3 << scale, width); // Not used.

  __m256i v32s = _mm256_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int offset = scale * 16;
  const __m256i vweight = _mm256_load_si256((const __m256i*) &intra_pdpc_w8_ver_weight[offset]);
  const __m256i vtopleft = _mm256_set1_epi16((uint16_t)top_left);

  const __m128i vshuf = _mm_setr_epi8(
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01
  );

  // For width 8, height must be at least 2. Handle 2 lines at once.
  for (int y = 0; y < height; y += 2) {
    const uint16_t ref2 = *(uint16_t*)&ref_side[1 + y];
    __m128i vref = _mm_set1_epi16(ref2);
    vref = _mm_shuffle_epi8(vref, vshuf);
    __m256i vref16 = _mm256_cvtepu8_epi16(vref);

    __m128i vdst = _mm_load_si128((const __m128i*)(dst + y * width));
    __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);

    __m256i accu = _mm256_sub_epi16(vref16, vtopleft);
    accu = _mm256_mullo_epi16(vweight, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    _mm_store_si128((__m128i*)(dst + (y * width)), filtered);
  }
}

static void angular_pdpc_mode50_w16_avx2(uvg_pixel* dst, const uvg_pixel top_left, const uvg_pixel* ref_side, const int width, const int height, const int scale)
{
  int limit = MIN(3 << scale, width); // Not used.

  __m256i v32s = _mm256_set1_epi16(32);

  const int offset = scale * 16;
  const __m256i vweight = _mm256_load_si256((const __m256i*) &intra_pdpc_w16_ver_weight[offset]);
  const __m256i vtopleft = _mm256_set1_epi16((uint16_t)top_left);

  for (int y = 0; y < height; ++y) {
    __m256i vref = _mm256_set1_epi16((int16_t)ref_side[1 + y]);
    
    __m128i vdst = _mm_load_si128((const __m128i*)(dst + y * width));
    __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);

    __m256i accu = _mm256_sub_epi16(vref, vtopleft);
    accu = _mm256_mullo_epi16(vweight, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    _mm_store_si128((__m128i*)(dst + y * width), filtered);
  }
}

static void angular_pdpc_mode50_scale1_avx2(uvg_pixel* dst, const uvg_pixel top_left, const uvg_pixel* ref_side, const int width, const int height)
{
  //const int scale = 1;
  //int limit = MIN(3 << scale, width); // Not used.

  __m256i v32s = _mm256_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int offset = 16; // scale * 16
  const __m256i vweight = _mm256_load_si256((const __m256i*) & intra_pdpc_w8_ver_weight[offset]);
  const __m256i vtopleft = _mm256_set1_epi16((uint16_t)top_left);

  const __m128i vshuf = _mm_setr_epi8(
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01
  );

  const int log2w = uvg_g_convert_to_log2[width];
  __m128i vseq = _mm_setr_epi32(0x00, 0x00, 0x01, 0x00);
  __m128i vidx = _mm_slli_epi64(vseq, log2w);

  // For width 8, height must be at least 2. Handle 2 lines at once.
  for (int y = 0; y < height; y += 2) {
    const uint16_t ref2 = *(uint16_t*)&ref_side[1 + y];
    __m128i vref = _mm_set1_epi16(ref2);
    vref = _mm_shuffle_epi8(vref, vshuf);
    __m256i vref16 = _mm256_cvtepu8_epi16(vref);

    //__m128i vdst = _mm_load_si128((const __m128i*)(dst + y * width));
    __m128i vdst = _mm_i64gather_epi64((const int64_t*)(dst + y * width), vidx, 1);
    __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);

    __m256i accu = _mm256_sub_epi16(vref16, vtopleft);
    accu = _mm256_mullo_epi16(vweight, accu);
    accu = _mm256_add_epi16(accu, v32s);
    accu = _mm256_srai_epi16(accu, 6);
    accu = _mm256_add_epi16(vdst16, accu);

    __m128i lo = _mm256_castsi256_si128(accu);
    __m128i hi = _mm256_extracti128_si256(accu, 1);
    __m128i filtered = _mm_packus_epi16(lo, hi);

    //_mm_store_si128((__m128i*)(dst + (y * width)), filtered);
    *(uint64_t*)(dst + ((y + 0) * width)) = _mm_extract_epi64(filtered, 0);
    *(uint64_t*)(dst + ((y + 1) * width)) = _mm_extract_epi64(filtered, 1);
  }
}

static void uvg_angular_pred_avx2(
  const cu_loc_t* const cu_loc,
  const int_fast8_t intra_mode,
  const int_fast8_t channel_type,
  const uvg_pixel* const in_ref_above,
  const uvg_pixel* const in_ref_left,
  uvg_pixel* const dst,
  const uint8_t multi_ref_idx,
  const uint8_t isp_mode,
  const int cu_dim)

{
  // ISP_TODO: non-square block implementation, height is passed but not used
  int width = channel_type == COLOR_Y ? cu_loc->width : cu_loc->chroma_width;
  int height = channel_type == COLOR_Y ? cu_loc->height : cu_loc->chroma_height;
  const int log2_width = uvg_g_convert_to_log2[width];
  const int log2_height = uvg_g_convert_to_log2[height];

  assert((log2_width >= 2 && log2_width <= 6) && (log2_height >= 0 && log2_height <= 6));

  // For chroma blocks, height has to be at least 2
  if (channel_type != COLOR_Y) {
    assert(log2_height >= 1);
  }

  // Modes [-1, -14] and [67, 81] are wide angle modes
  assert(intra_mode >= -14 && intra_mode <= 81);

  uint8_t multi_ref_index = channel_type == COLOR_Y ? multi_ref_idx : 0;
  uint8_t isp = isp_mode;

  static const int16_t modedisp2sampledisp[32] = { 0,    1,    2,    3,    4,    6,     8,   10,   12,   14,   16,   18,   20,   23,   26,   29,   32,   35,   39,  45,  51,  57,  64,  73,  86, 102, 128, 171, 256, 341, 512, 1024 };
  static const int16_t modedisp2invsampledisp[32] = { 0, 16384, 8192, 5461, 4096, 2731, 2048, 1638, 1365, 1170, 1024, 910, 819, 712, 630, 565, 512, 468, 420, 364, 321, 287, 256, 224, 191, 161, 128, 96, 64, 48, 32, 16 }; // (512 * 32) / sampledisp
  static const int32_t pre_scale[] = { 8, 7, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, -1, -1, -2, -3 };

  // Temporary buffer for modes 11-25.
  // It only needs to be big enough to hold indices from -width to width-1.
  uvg_pixel temp_main[2 * 128 + 3 + 33 * MAX_REF_LINE_IDX] = { 0 };
  uvg_pixel temp_side[2 * 128 + 3 + 33 * MAX_REF_LINE_IDX] = { 0 };

  int32_t pred_mode = intra_mode; // ToDo: handle WAIP

  // Whether to swap references to always project on the left reference row.
  const bool vertical_mode = intra_mode >= 34;
  
  // Modes distance to horizontal or vertical mode. Possible values: [-16, 16]
  // For pure vertical or horizontal modes, this is 0. For pure diagonal modes, this is either -16 or 16.
  const int_fast8_t mode_disp = vertical_mode ? pred_mode - 50 : -(pred_mode - 18);
  const bool wide_angle_mode = mode_disp > 16;

  // Sample displacement per column in fractions of 32.
  const int_fast16_t sample_disp = (mode_disp < 0 ? -1 : 1) * modedisp2sampledisp[abs(mode_disp)];

  const int side_size = vertical_mode ? log2_height : log2_width;
  int scale = MIN(2, side_size - pre_scale[abs(mode_disp)]);

  // Pointer for the reference we are interpolating from.
  uvg_pixel* ref_main;
  // Pointer for the other reference.
  const uvg_pixel* ref_side;

  const int top_ref_length = isp_mode == ISP_MODE_VER ? width + cu_dim : width << 1;
  const int left_ref_length = isp_mode == ISP_MODE_HOR ? height + cu_dim : height << 1;

  // Set ref_main and ref_side such that, when indexed with 0, they point to
  // index 0 in block coordinates.
  if (sample_disp < 0) {
    memcpy(&temp_main[height], &in_ref_above[0], (width + 2 + multi_ref_index) * sizeof(uvg_pixel));
    memcpy(&temp_side[width], &in_ref_left[0], (height + 2 + multi_ref_index) * sizeof(uvg_pixel));

    ref_main = vertical_mode ? temp_main + height : temp_side + width;
    ref_side = vertical_mode ? temp_side + width  : temp_main + height;

    int size_side = vertical_mode ? height : width;
    for (int i = -size_side; i <= -1; i++) {
      ref_main[i] = ref_side[MIN((-i * modedisp2invsampledisp[abs(mode_disp)] + 256) >> 9, size_side)];
    }
  }
  else {
    memcpy(&temp_main[0], &in_ref_above[0], (top_ref_length + 1 + multi_ref_index) * sizeof(uvg_pixel));
    memcpy(&temp_side[0], &in_ref_left[0], (left_ref_length + 1 + multi_ref_index) * sizeof(uvg_pixel));

    ref_main = vertical_mode ? temp_main : temp_side;
    ref_side = vertical_mode ? temp_side : temp_main;

    const int log2_ratio = log2_width - log2_height;
    const int s = MAX(0, vertical_mode ? log2_ratio : -log2_ratio);
    const int max_index = (multi_ref_index << s) + 2;
    int ref_length;
    if (isp_mode) {
      ref_length = vertical_mode ? top_ref_length : left_ref_length;
    }
    else {
      ref_length = vertical_mode ? width << 1 : height << 1;
    }
    const uvg_pixel val = ref_main[ref_length + multi_ref_index];
    for (int j = 1; j <= max_index; j++) {
      ref_main[ref_length + multi_ref_index + j] = val;
    }
  }

  // compensate for line offset in reference line buffers
  ref_main += multi_ref_index;
  ref_side += multi_ref_index;
  //if (!vertical_mode) { SWAP(width, height, int) }

  static const int uvg_intra_hor_ver_dist_thres[8] = { 24, 24, 24, 14, 2, 0, 0, 0 };
  int filter_threshold = uvg_intra_hor_ver_dist_thres[(log2_width + log2_height) >> 1];
  int dist_from_vert_or_hor = MIN(abs((int32_t)pred_mode - 50), abs((int32_t)pred_mode - 18));

  bool use_cubic = true; // Default to cubic filter
  if (dist_from_vert_or_hor > filter_threshold) {
    if ((abs(sample_disp) & 0x1F) != 0)
    {
      use_cubic = false;
    }
  }
  // Cubic must be used if ref line != 0 or if isp mode != 0
  if (multi_ref_index || isp) {
    use_cubic = true;
  }



  if (sample_disp != 0) {
    // The mode is not horizontal or vertical, we have to do interpolation.

    // Set delta table pointers
    const int table_offset = wide_angle_mode ? (pred_mode < 2 ? (pred_mode + 13) * 64 : (81 - pred_mode) * 64) : (pred_mode <= 34 ? (pred_mode - 2) * 64 : (66 - pred_mode) * 64);
    const int16_t* delta_int   = wide_angle_mode ? &delta_int_wide_angle_table[table_offset] : &delta_int_table[table_offset];
    delta_int += multi_ref_index; // TODO: This are not necessarily large enough for 64 dimension blocks
    const int16_t* delta_fract = wide_angle_mode ? &delta_fract_wide_angle_table[table_offset] : &delta_fract_table[table_offset];
    delta_fract += multi_ref_index;

    // Check if the angle is fractional. If yes, interpolation is needed
    if ((abs(sample_disp) & 0x1F) != 0) {

      // Luma Channel
      if (channel_type == 0) {
        if (vertical_mode) {
          switch (width) {
            case  4: angular_pred_w4_ver_avx2(dst, ref_main, delta_int, delta_fract, height, use_cubic); break;
            case  8: angular_pred_w8_ver_avx2(dst, ref_main, delta_int, delta_fract, height, use_cubic); break;
            case 16: angular_pred_w16_ver_avx2(dst, ref_main, delta_int, delta_fract, width, height, use_cubic); break;
            case 32: angular_pred_w16_ver_avx2(dst, ref_main, delta_int, delta_fract, width, height, use_cubic); break;
            case 64: angular_pred_w16_ver_avx2(dst, ref_main, delta_int, delta_fract, width, height, use_cubic); break;
            default:
              assert(false && "Intra angular predicion: illegal width.\n");
              break;
          }
        }
        else {
          switch (width) {
            case  4: angular_pred_w4_hor_avx2(dst, ref_main, delta_int, delta_fract, height, use_cubic); break;
            case  8: angular_pred_w8_hor_avx2(dst, ref_main, delta_int, delta_fract, height, use_cubic); break;
            case 16: angular_pred_w16_hor_avx2(dst, ref_main, delta_int, delta_fract, width, height, use_cubic); break;
            case 32: angular_pred_w16_hor_avx2(dst, ref_main, delta_int, delta_fract, width, height, use_cubic); break;
            case 64: angular_pred_w16_hor_avx2(dst, ref_main, delta_int, delta_fract, width, height, use_cubic); break;
            default:
              assert(false && "Intra angular predicion: illegal width.\n");
              break;
          }
        }
      }
      // Chroma channels
      else {
        // Do 2-tap linear filtering for chroma channels
        
        if (vertical_mode) {
          switch (width) {
          // No wide angle handling for w4 is needed.
          case  4: angular_pred_linear_filter_w4_ver_avx2(dst, ref_main, height, delta_int, pred_mode); break;
          case  8: angular_pred_linear_filter_w8_ver_avx2(dst, ref_main, height, delta_int, pred_mode); break;
          case 16: angular_pred_linear_filter_w16_ver_avx2(dst, ref_main, height, delta_int, pred_mode); break;
          case 32: angular_pred_linear_filter_w32_ver_avx2(dst, ref_main, height, delta_int, pred_mode); break;
          default:
            assert(false && "Intra angular predicion: illegal chroma width.\n");
            break;
          }
        }
        else {
          if (wide_angle_mode) {
            switch (width) {
            case  4: angular_pred_linear_filter_w4_hor_wide_angle_avx2(dst, ref_main, height, pred_mode, delta_int); break;
            case  8: angular_pred_linear_filter_w8_hor_wide_angle_avx2(dst, ref_main, height, pred_mode, delta_int, delta_fract); break;
            case 16: angular_pred_linear_filter_w16_hor_wide_angle_avx2(dst, ref_main, height, pred_mode, delta_int, delta_fract); break;
            case 32: assert(false && "This code branch only works with UVG_FORMAT_P420."); break; // This branch is never executed with UVG_FORMAT_P420, due to chroma being only 32 width or height.
            default:
              assert(false && "Intra angular predicion: illegal chroma width.\n");
              break;
            }
          }
          else {
            switch (width) {
            case  4: angular_pred_linear_filter_w4_hor_avx2(dst, ref_main, height, pred_mode, delta_int); break;
            case  8: angular_pred_linear_filter_w8_hor_avx2(dst, ref_main, height, pred_mode, delta_int); break;
            case 16: angular_pred_linear_filter_w16_hor_avx2(dst, ref_main, height, pred_mode, delta_int); break;
            case 32: angular_pred_linear_filter_w32_hor_avx2(dst, ref_main, height, pred_mode, delta_int); break;
            default:
              assert(false && "Intra angular predicion: illegal chroma width.\n");
              break;
            }
          }
        }
      }
    }
    else {
      // No interpolation or filtering needed, just copy the integer samples
      if (vertical_mode) {
        angular_pred_non_fractional_angle_pxl_copy_ver_avx2(dst, ref_main, width, height, delta_int);
      }
      else {
        angular_pred_non_fractional_angle_pxl_copy_hor_avx2(dst, ref_main, width, height, delta_int);
      }
    }
  }
  else {
    // Mode is horizontal or vertical, just copy the pixels.
    if (vertical_mode) {
      for (int_fast32_t y = 0; y < height; ++y) {
        switch (width) {
          case 4:  memcpy(&dst[y * 4],  &ref_main[1],  4 * sizeof(uvg_pixel)); break;
          case 8:  memcpy(&dst[y * 8],  &ref_main[1],  8 * sizeof(uvg_pixel)); break;
          case 16: memcpy(&dst[y * 16], &ref_main[1], 16 * sizeof(uvg_pixel)); break;
          case 32: memcpy(&dst[y * 32], &ref_main[1], 32 * sizeof(uvg_pixel)); break;
          case 64: memcpy(&dst[y * 64], &ref_main[1], 64 * sizeof(uvg_pixel)); break;
        }
      }
    }
    else {
    #define UNROLL(w, h) \
      if ((h) == height && (w) == width) { \
        for (int y = 0; y < (h); ++y) { \
          const __m128i vdst = _mm_set1_epi8(ref_main[y + 1]); \
          switch ((w)) {\
            case 4:  _mm_storeu_si32((__m128i*) &dst[y * 4], vdst); break;\
            case 8:  _mm_storeu_si64((__m128i*) &dst[y * 8], vdst); break;\
            case 16: _mm_store_si128((__m128i*) &dst[y * 16], vdst); break;\
            case 32:\
              _mm_store_si128((__m128i*) &dst[y * 32 +  0], vdst);\
              _mm_store_si128((__m128i*) &dst[y * 32 + 16], vdst);\
              break;\
            case 64: \
              _mm_store_si128((__m128i*) &dst[y * 64 +  0], vdst);\
              _mm_store_si128((__m128i*) &dst[y * 64 + 16], vdst);\
              _mm_store_si128((__m128i*) &dst[y * 64 + 32], vdst);\
              _mm_store_si128((__m128i*) &dst[y * 64 + 48], vdst);\
              break;  \
            default:\
              assert(false && "Intra angular predicion: illegal width.\n");\
              break;\
          }\
        } \
      }
      UNROLL(4, 4);
      UNROLL(4, 8);
      UNROLL(4, 16);
      UNROLL(4, 32);
      UNROLL(4, 64);
      UNROLL(8, 2);
      UNROLL(8, 4);
      UNROLL(8, 8);
      UNROLL(8, 16);
      UNROLL(8, 32);
      UNROLL(8, 64);
      UNROLL(16, 1);
      UNROLL(16, 2);
      UNROLL(16, 4);
      UNROLL(16, 8);
      UNROLL(16, 16);
      UNROLL(16, 32);
      UNROLL(16, 64);
      UNROLL(32, 1);
      UNROLL(32, 2);
      UNROLL(32, 4);
      UNROLL(32, 8);
      UNROLL(32, 16);
      UNROLL(32, 32);
      UNROLL(32, 64);
      UNROLL(64, 1);
      UNROLL(64, 2);
      UNROLL(64, 4);
      UNROLL(64, 8);
      UNROLL(64, 16);
      UNROLL(64, 32);
      UNROLL(64, 64);
      #undef UNROLL
    }
  }

  
  bool PDPC_filter = (width >= TR_MIN_WIDTH && height >= TR_MIN_WIDTH);
  if (pred_mode > 1 && pred_mode < 67) {
    // Disable PDPC filter if both references are used or if MRL is used
    if (mode_disp < 0 || multi_ref_index) {
      PDPC_filter = false;
    }
    else if (mode_disp > 0) {
      // If scale is negative, PDPC filtering has no effect, therefore disable it.
      PDPC_filter &= (scale >= 0);
    }
  }
  if (PDPC_filter) {
    // Handle pure horizontal and vertical with separate PDPC solution
    if (pred_mode == 18) {
      scale = (log2_width + log2_height - 2) >> 2;
      const uvg_pixel top_left = ref_main[0];
      
      switch (width) {
        case 4:  angular_pdpc_mode18_w4_avx2(dst, top_left, ref_side, height, scale); break;
        case 8:  angular_pdpc_mode18_w8_avx2(dst, top_left, ref_side, height, scale); break;
        case 16: angular_pdpc_mode18_w16_avx2(dst, top_left, ref_side, height, scale); break;
        case 32: angular_pdpc_mode18_w32_avx2(dst, top_left, ref_side, height, scale); break;
        case 64: angular_pdpc_mode18_w64_avx2(dst, top_left, ref_side, height, scale); break;
        default:
          assert(false && "Intra PDPC, invalid width.\n");
          break;
      }
    }
    else if (pred_mode == 50) {
      scale = (log2_width + log2_height - 2) >> 2;
      const uvg_pixel top_left = ref_main[0];
      switch (width) {
        case 4:  angular_pdpc_mode50_w4_avx2(dst, top_left, ref_side, height, scale); break;
        case 8:  angular_pdpc_mode50_w8_avx2(dst, top_left, ref_side, height, scale); break;
        case 16: // 16 and higher handled by same functions.
        case 32: 
        case 64: 
          if (scale == 1) {
            angular_pdpc_mode50_scale1_avx2(dst, top_left, ref_side, width, height);
          }
          else {
            angular_pdpc_mode50_w16_avx2(dst, top_left, ref_side, width, height, scale);
          }
          break;
        default:
          assert(false && "Intra PDPC, invalid width.\n");
          break;
      }
    }
    else {
      if (vertical_mode) {
        // Note: no need to check for negative mode_disp, since it is already checked before.
        switch (width) {
        case 4:
          // Low mode disp -> low angle. For pdpc, this causes the needed references to be extremely sparse making loads without using gathers impossible.
          // Handle high angles with more tight reference spacing with separate functions with more optimized loads.
          if (mode_disp < 6)
            angular_pdpc_ver_w4_avx2(dst, ref_side, height, scale, mode_disp);
          else
            angular_pdpc_ver_w4_high_angle_avx2(dst, ref_side, height, scale, mode_disp);
          break;
        case 8:
          if (scale == 0) {
            if (mode_disp < 6)
              angular_pdpc_ver_4x4_scale0_avx2(dst, ref_side, width, height, mode_disp);
            else
              angular_pdpc_ver_4x4_scale0_high_angle_avx2(dst, ref_side, width, height, mode_disp);
          }
          else if (scale == 1) {
            if (mode_disp < 8)
              angular_pdpc_ver_8x2_scale1_avx2(dst, ref_side, width, height, mode_disp);
            else
              angular_pdpc_ver_8x2_scale1_high_angle_avx2(dst, ref_side, width, height, mode_disp);
          }
          else {
            if (mode_disp < 10)
              angular_pdpc_ver_w8_avx2(dst, ref_side, height, scale, mode_disp);
            else
              angular_pdpc_ver_8x2_scale2_high_angle_avx2(dst, ref_side, width, height, mode_disp);
          }
          break;
        case 16: // 16 width and higher done with the same functions
        case 32:
        case 64:
          switch (scale) {
          case 0:
            if (mode_disp < 6)
              angular_pdpc_ver_4x4_scale0_avx2(dst, ref_side, width, height, mode_disp);
            else
              angular_pdpc_ver_4x4_scale0_high_angle_avx2(dst, ref_side, width, height, mode_disp);
            break;
          case 1:
            if (mode_disp < 8)
              angular_pdpc_ver_8x2_scale1_avx2(dst, ref_side, width, height, mode_disp);
            else
              angular_pdpc_ver_8x2_scale1_high_angle_avx2(dst, ref_side, width, height, mode_disp);
            break;
          case 2:
            if (mode_disp < 14)
              angular_pdpc_ver_w16_avx2(dst, ref_side, width, height, mode_disp);
            else
              angular_pdpc_ver_w16_scale2_high_angle_avx2(dst, ref_side, width, height, mode_disp);
            break;
          default:
            assert(false && "Intra PDPC: Invalid scale.\n");
          }
          break;
        default:
          assert(false && "Intra PDPC: Invalid width.\n");
        }
      }
      else {
        switch (width) {
        case 4:
          // Low mode disp -> low angle. For pdpc, this causes the needed references to be extremely sparse making loads without using gathers impossible.
          // Handle high angles with more tight reference spacing with separate functions with more optimized loads.
          if (mode_disp < 6)
            angular_pdpc_hor_w4_avx2(dst, ref_side, height, scale, mode_disp);
          else
            angular_pdpc_hor_w4_high_angle_avx2(dst, ref_side, height, scale, mode_disp);
          break;
        case 8:  angular_pdpc_hor_w8_avx2(dst, ref_side, height, scale, mode_disp); break;
        case 16: // 16 width and higher done with the same function
        case 32:
        case 64: angular_pdpc_hor_w16_avx2(dst, ref_side, width, height, scale, mode_disp); break;
        default:
          assert(false && "Intra PDPC: Invalid width.\n");
        }
      }
    }
  }
}


typedef void (intra_planar_half_func)(const uvg_pixel* ref_main, const uvg_pixel* ref_side, const int line, const int shift, __m256i* dst);

// w1 and w2 for planar horizontal do not exist, since intra prediction must be at least of width 4
// Also worth noting is that minimum amount of samples must be 16, 
// therefore the smallest possible predictions are 4x4, 8x2 and 16x1
static void intra_pred_planar_hor_w4(const uvg_pixel* ref, const uvg_pixel* ref_side, const int line, const int shift, __m256i* dst)
{
  const __m256i v_last_ref = _mm256_set1_epi16(ref_side[4 + 1]);

  const __m256i v_ref_coeff = _mm256_setr_epi16(3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0);
  const __m256i v_last_ref_coeff = _mm256_setr_epi16(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4);

  const __m256i v_last_ref_mul = _mm256_mullo_epi16(v_last_ref, v_last_ref_coeff);
  const __m256i shuffle_mask = _mm256_setr_epi8(0, -1, 0, -1, 0, -1, 0, -1, 8, -1, 8, -1, 8, -1, 8, -1, 0, -1, 0, -1, 0, -1, 0, -1, 8, -1, 8, -1, 8, -1, 8, -1);

  // Handle 4 lines at a time
  #define UNROLL_LOOP(num) \
  for (int i = 0, d = 0; i < (num); i += 4, ++d) { \
    /* | ref1 | ref2 | ref3 | ref4 | Don't care*/ \
    __m128i v_ref_0 = _mm_loadu_si128((__m128i const*)& ref[i + 1]); \
    /* | ref1 | 0 * 7 | ref2 | 0 * 7 | ref3 | 0 * 7 | ref4 | 0* 7 | */ \
    __m256i v_ref = _mm256_cvtepu8_epi64(v_ref_0); \
    /* | ref1_l | ref1_h | ref1_l | ref1_h | ... */ \
    v_ref = _mm256_shuffle_epi8(v_ref, shuffle_mask); \
    \
    __m256i v_tmp = _mm256_mullo_epi16(v_ref, v_ref_coeff); \
    \
    dst[d] = _mm256_add_epi16(v_last_ref_mul, v_tmp); \
  }

  switch (line) {
  case 1: UNROLL_LOOP(1); break;
  case 2: UNROLL_LOOP(2); break;
  case 4: UNROLL_LOOP(4); break;
  case 8: UNROLL_LOOP(8); break;
  case 16: UNROLL_LOOP(16); break;
  case 32: UNROLL_LOOP(32); break;
  case 64: UNROLL_LOOP(64); break;
  default:
    assert(false && "Invalid dimension.");
    break;
  }
  #undef UNROLL_LOOP
}
static void intra_pred_planar_hor_w8(const uvg_pixel* ref, const uvg_pixel* ref_side, const int line, const int shift, __m256i* dst)
{
  const __m256i v_last_ref = _mm256_set1_epi16(ref_side[8 + 1]);

  const __m256i v_ref_coeff = _mm256_setr_epi16(7, 6, 5, 4, 3, 2, 1, 0, 7, 6, 5, 4, 3, 2, 1, 0);
  const __m256i v_last_ref_coeff = _mm256_setr_epi16(1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8);

  const __m256i v_last_ref_mul = _mm256_mullo_epi16(v_last_ref, v_last_ref_coeff);

  // Handle 2 lines at a time
  #define UNROLL_LOOP(num) \
  for (int i = 0, d = 0; i < (num); i += 2, ++d) { \
    __m128i v_ref0 = _mm_set1_epi16(ref[i + 1]); \
    __m128i v_ref1 = _mm_set1_epi16(ref[i + 2]); \
    \
    __m256i v_ref = _mm256_castsi128_si256(v_ref0); \
    v_ref = _mm256_inserti128_si256(v_ref, v_ref1, 1); \
    \
    __m256i v_tmp = _mm256_mullo_epi16(v_ref, v_ref_coeff); \
    \
    dst[d] = _mm256_add_epi16(v_last_ref_mul, v_tmp); \
  }

  switch (line) {
  case 1: UNROLL_LOOP(1); break;
  case 2: UNROLL_LOOP(2); break;
  case 4: UNROLL_LOOP(4); break;
  case 8: UNROLL_LOOP(8); break;
  case 16: UNROLL_LOOP(16); break;
  case 32: UNROLL_LOOP(32); break;
  case 64: UNROLL_LOOP(64); break;
  default:
    assert(false && "Invalid dimension.");
    break;
  }
  #undef UNROLL_LOOP
}
static void intra_pred_planar_hor_w16(const uvg_pixel* ref, const uvg_pixel* ref_side, const int line, const int shift, __m256i* dst)
{
  const __m256i v_last_ref = _mm256_set1_epi16(ref_side[16 + 1]);

  const __m256i v_ref_coeff = _mm256_setr_epi16(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
  const __m256i v_last_ref_coeff = _mm256_setr_epi16(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);

  const __m256i v_last_ref_mul = _mm256_mullo_epi16(v_last_ref, v_last_ref_coeff);

  #define UNROLL_LOOP(num) \
  for (int i = 0, d = 0; i < (num); ++i, ++d) { \
    __m256i v_ref = _mm256_set1_epi16(ref[i + 1]); \
    __m256i v_tmp = _mm256_mullo_epi16(v_ref, v_ref_coeff); \
    dst[d] = _mm256_add_epi16(v_last_ref_mul, v_tmp); \
  }

  switch (line) {
  case 1: UNROLL_LOOP(1); break;
  case 2: UNROLL_LOOP(2); break;
  case 4: UNROLL_LOOP(4); break;
  case 8: UNROLL_LOOP(8); break;
  case 16: UNROLL_LOOP(16); break;
  case 32: UNROLL_LOOP(32); break;
  case 64: UNROLL_LOOP(64); break;
  default:
    assert(false && "Invalid dimension.");
    break;
  }
  #undef UNROLL_LOOP
}
static void intra_pred_planar_hor_w32(const uvg_pixel* ref, const uvg_pixel* ref_side, const int line, const int shift, __m256i* dst)
{
  const __m256i v_last_ref = _mm256_set1_epi16(ref_side[32 + 1]);

  const __m256i v_ref_coeff0 = _mm256_setr_epi16(31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16);
  const __m256i v_ref_coeff1 = _mm256_setr_epi16(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);

  const __m256i v_last_ref_coeff0 = _mm256_setr_epi16(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
  const __m256i v_last_ref_coeff1 = _mm256_setr_epi16(17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32);

  const __m256i v_last_ref_mul0 = _mm256_mullo_epi16(v_last_ref, v_last_ref_coeff0);
  const __m256i v_last_ref_mul1 = _mm256_mullo_epi16(v_last_ref, v_last_ref_coeff1);

  #define UNROLL_LOOP(num) \
  for (int i = 0, d = 0; i < (num); ++i, d += 2) { \
    __m256i v_ref = _mm256_set1_epi16(ref[i + 1]); \
    __m256i v_tmp0 = _mm256_mullo_epi16(v_ref, v_ref_coeff0); \
    __m256i v_tmp1 = _mm256_mullo_epi16(v_ref, v_ref_coeff1); \
    dst[d + 0] = _mm256_add_epi16(v_last_ref_mul0, v_tmp0); \
    dst[d + 1] = _mm256_add_epi16(v_last_ref_mul1, v_tmp1); \
  }

  switch (line) {
  case 1: UNROLL_LOOP(1); break;
  case 2: UNROLL_LOOP(2); break;
  case 4: UNROLL_LOOP(4); break;
  case 8: UNROLL_LOOP(8); break;
  case 16: UNROLL_LOOP(16); break;
  case 32: UNROLL_LOOP(32); break;
  case 64: UNROLL_LOOP(64); break;
  default:
    assert(false && "Invalid dimension.");
    break;
  }
  #undef UNROLL_LOOP
}
static void intra_pred_planar_hor_w64(const uvg_pixel* ref, const uvg_pixel* ref_side, const int line, const int shift, __m256i* dst)
{
  const __m256i v_last_ref = _mm256_set1_epi16(ref_side[64 + 1]);

  const __m256i v_ref_coeff0 = _mm256_setr_epi16(63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48);
  const __m256i v_ref_coeff1 = _mm256_setr_epi16(47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32);
  const __m256i v_ref_coeff2 = _mm256_setr_epi16(31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16);
  const __m256i v_ref_coeff3 = _mm256_setr_epi16(15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0);

  const __m256i v_last_ref_coeff0 = _mm256_setr_epi16( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16);
  const __m256i v_last_ref_coeff1 = _mm256_setr_epi16(17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32);
  const __m256i v_last_ref_coeff2 = _mm256_setr_epi16(33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48);
  const __m256i v_last_ref_coeff3 = _mm256_setr_epi16(49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64);

  const __m256i v_last_ref_mul0 = _mm256_mullo_epi16(v_last_ref, v_last_ref_coeff0);
  const __m256i v_last_ref_mul1 = _mm256_mullo_epi16(v_last_ref, v_last_ref_coeff1);
  const __m256i v_last_ref_mul2 = _mm256_mullo_epi16(v_last_ref, v_last_ref_coeff2);
  const __m256i v_last_ref_mul3 = _mm256_mullo_epi16(v_last_ref, v_last_ref_coeff3);

  for (int i = 0, d = 0; i < line; ++i, d += 4) {
    __m256i v_ref = _mm256_set1_epi16(ref[i + 1]);
    __m256i v_tmp0 = _mm256_mullo_epi16(v_ref, v_ref_coeff0);
    __m256i v_tmp1 = _mm256_mullo_epi16(v_ref, v_ref_coeff1);
    __m256i v_tmp2 = _mm256_mullo_epi16(v_ref, v_ref_coeff2);
    __m256i v_tmp3 = _mm256_mullo_epi16(v_ref, v_ref_coeff3);
    dst[d + 0] = _mm256_add_epi16(v_last_ref_mul0, v_tmp0);
    dst[d + 1] = _mm256_add_epi16(v_last_ref_mul1, v_tmp1);
    dst[d + 2] = _mm256_add_epi16(v_last_ref_mul2, v_tmp2);
    dst[d + 3] = _mm256_add_epi16(v_last_ref_mul3, v_tmp3);
  }
}

static void intra_pred_planar_ver_w4(const uvg_pixel* ref, const uvg_pixel* ref_side, const int line, const int shift, __m256i* dst)
{
  const __m256i v_last_ref = _mm256_set1_epi8(ref_side[line + 1]);

  // Overflow possible for this width if line > 32
  const bool overflow = line > 32;

  // Got four 8-bit references, or 32 bits of data. Duplicate to fill a whole 256-bit vector.
  const uint32_t* tmp = (const uint32_t*)&ref[1]; // Cast to 32 bit int to load 4 refs at the same time
  const __m256i v_ref = _mm256_set1_epi32(*tmp);

  const __m256i* v_ys = (const __m256i*)planar_avx2_ver_w4ys;

  // Table offset
  int offset;
  switch (line) {
    case 64: offset = 0;  break;
    case 32: offset = 16; break;
    case 16: offset = 24; break;
    case  8: offset = 28; break;
    case  4: offset = 30; break;
    default:
      assert(false && "Invalid height for width 4.");
      break;
  }

  // Handle 4 lines at a time
  #define UNROLL_LOOP(num) \
  for (int y = 0, s = offset, d = 0; y < (num); y += 4, ++s, ++d) { \
    __m256i v_lo = _mm256_unpacklo_epi8(v_ref, v_last_ref); \
    dst[d] = _mm256_maddubs_epi16(v_lo, v_ys[s]); \
  }

  switch (line) {
    case 1: UNROLL_LOOP(1); break;
    case 2: UNROLL_LOOP(2); break;
    case 4: UNROLL_LOOP(4); break;
    case 8: UNROLL_LOOP(8); break;
    case 16: UNROLL_LOOP(16); break;
    case 32: UNROLL_LOOP(32); break;
    case 64: UNROLL_LOOP(64); break;
    default:
      assert(false && "Invalid dimension.");
      break;
  }
  #undef UNROLL_LOOP
}
static void intra_pred_planar_ver_w8(const uvg_pixel* ref, const uvg_pixel* ref_side, const int line, const int shift, __m256i* dst)
{
  const __m256i v_last_ref = _mm256_set1_epi8(ref_side[line + 1]);
  
  // Got eight 8-bit samples, or 64 bits of data. Duplicate to fill a whole 256-bit vector.
  const __m128i v_ref_raw = _mm_loadu_si128((const __m128i*)&ref[1]);
  __m256i v_ref = _mm256_castsi128_si256(v_ref_raw);
  v_ref = _mm256_inserti128_si256(v_ref, v_ref_raw, 1);
  v_ref = _mm256_shuffle_epi32(v_ref, _MM_SHUFFLE(1, 1, 0, 0));

  const __m256i* v_ys = (const __m256i*)planar_avx2_ver_w4ys;

  // Table offset
  int offset;
  switch (line) {
  case 64: offset = 0;  break;
  case 32: offset = 16; break;
  case 16: offset = 24; break;
  case  8: offset = 28; break;
  case  4: offset = 30; break;
  case  2: offset = 31; break;
  default:
    assert(false && "Invalid height for width 8.");
    break;
  }

  // Handle 4 lines at a time
  #define UNROLL_LOOP(num) \
  for (int y = 0, s = offset, d = 0; y < (num); y += 4, ++s, d += 2) { \
    __m256i v_lo = _mm256_unpacklo_epi8(v_ref, v_last_ref); \
    __m256i v_hi = _mm256_unpackhi_epi8(v_ref, v_last_ref); \
    \
    __m256i v_madd_lo = _mm256_maddubs_epi16(v_lo, v_ys[s]); \
    __m256i v_madd_hi = _mm256_maddubs_epi16(v_hi, v_ys[s]); \
    __m256i v_tmp0 = _mm256_permute2x128_si256(v_madd_lo, v_madd_hi, 0x20); \
    __m256i v_tmp1 = _mm256_permute2x128_si256(v_madd_lo, v_madd_hi, 0x31); \
    \
    dst[d + 0] = _mm256_permute4x64_epi64(v_tmp0, _MM_SHUFFLE(3, 1, 2, 0)); \
    dst[d + 1] = _mm256_permute4x64_epi64(v_tmp1, _MM_SHUFFLE(3, 1, 2, 0)); \
  }

  switch (line) {
  case 1: UNROLL_LOOP(1); break;
  case 2: UNROLL_LOOP(2); break;
  case 4: UNROLL_LOOP(4); break;
  case 8: UNROLL_LOOP(8); break;
  case 16: UNROLL_LOOP(16); break;
  case 32: UNROLL_LOOP(32); break;
  case 64: UNROLL_LOOP(64); break;
  default:
    assert(false && "Invalid dimension.");
    break;
  }
  #undef UNROLL_LOOP
}
static void intra_pred_planar_ver_w16(const uvg_pixel* ref, const uvg_pixel* ref_side, const int line, const int shift, __m256i* dst)
{
  const __m256i v_last_ref = _mm256_set1_epi8(ref_side[line + 1]);

  // Got 16 8-bit samples, or 128 bits of data. Duplicate to fill a whole 256-bit vector.
  const __m128i v_ref_raw = _mm_loadu_si128((const __m128i*) &ref[1]);
  __m256i v_ref = _mm256_castsi128_si256(v_ref_raw);
  v_ref = _mm256_inserti128_si256(v_ref, v_ref_raw, 1);

  const __m256i* v_ys = (const __m256i*)planar_avx2_ver_w8ys;

  // Table offset
  int offset;
  switch (line) {
  case 64: offset = 0;  break;
  case 32: offset = 32; break;
  case 16: offset = 48; break;
  case  8: offset = 56; break;
  case  4: offset = 60; break;
  case  2: offset = 62; break;
  case  1: offset = 64; break;
  default:
    assert(false && "Invalid height for width 16.");
    break;
  }

  // Calculations for cases where line > 2
  // These stay constant through the loop
  const __m256i v_lo = _mm256_unpacklo_epi8(v_ref, v_last_ref);
  const __m256i v_hi = _mm256_unpackhi_epi8(v_ref, v_last_ref);

  // Handle 2 lines at a time
  #define UNROLL_LOOP(num) \
  for (int y = 0, s = offset; y < line; y += 2, ++s) { \
    __m256i v_madd_lo = _mm256_maddubs_epi16(v_lo, v_ys[s]); \
    __m256i v_madd_hi = _mm256_maddubs_epi16(v_hi, v_ys[s]); \
    dst[y + 0] = _mm256_permute2x128_si256(v_madd_lo, v_madd_hi, 0x20); \
    dst[y + 1] = _mm256_permute2x128_si256(v_madd_lo, v_madd_hi, 0x31); \
  }

  __m256i v_tmp;
  switch (line) {
  case 1:
    // Specialized calculation for line == 1
    v_tmp = _mm256_permute2x128_si256(v_lo, v_hi, 0x20);
    dst[0] = _mm256_maddubs_epi16(v_tmp, v_ys[offset + 0]);
    break;
  case 2: 
    // Specialized calculation for line == 2
    v_tmp = _mm256_permute2x128_si256(v_lo, v_hi, 0x20);
    dst[0] = _mm256_maddubs_epi16(v_tmp, v_ys[offset + 0]);
    dst[1] = _mm256_maddubs_epi16(v_tmp, v_ys[offset + 1]);
    break;
  case 4: UNROLL_LOOP(4); break;
  case 8: UNROLL_LOOP(8); break;
  case 16: UNROLL_LOOP(16); break;
  case 32: UNROLL_LOOP(32); break;
  case 64: UNROLL_LOOP(64); break;
  default:
    assert(false && "Invalid dimension.");
    break;
  }
#undef UNROLL_LOOP
}
static void intra_pred_planar_ver_w32(const uvg_pixel* ref, const uvg_pixel* ref_side, const int line, const int shift, __m256i* dst)
{
  const __m256i v_last_ref = _mm256_set1_epi8(ref_side[line + 1]);

  // Got 32 8-bit samples, or 256 bits of data. Load into a single vector
  const __m256i v_ref = _mm256_loadu_si256((const __m256i*) &ref[1]);

  // These stay constant through the loop
  const __m256i v_lo = _mm256_unpacklo_epi8(v_ref, v_last_ref);
  const __m256i v_hi = _mm256_unpackhi_epi8(v_ref, v_last_ref);

  #define UNROLL_LOOP(num) \
  for (uint8_t y = 0, a = (num) - 1, b = 1, d = 0; y < (num); ++y, --a, ++b, d += 2) { \
    uint8_t tmp[2] = {a, b}; \
    uint16_t* tmp2 = (uint16_t*)tmp; \
    const __m256i v_ys = _mm256_set1_epi16(*tmp2); \
    \
    __m256i v_madd_lo = _mm256_maddubs_epi16(v_lo, v_ys); \
    __m256i v_madd_hi = _mm256_maddubs_epi16(v_hi, v_ys); \
    dst[d + 0] = _mm256_permute2x128_si256(v_madd_lo, v_madd_hi, 0x20); \
    dst[d + 1] = _mm256_permute2x128_si256(v_madd_lo, v_madd_hi, 0x31); \
  }

  switch (line) {
  case 1: UNROLL_LOOP(1); break;
  case 2: UNROLL_LOOP(2); break;
  case 4: UNROLL_LOOP(4); break;
  case 8: UNROLL_LOOP(8); break;
  case 16: UNROLL_LOOP(16); break;
  case 32: UNROLL_LOOP(32); break;
  case 64: UNROLL_LOOP(64); break;
  default:
    assert(false && "Invalid dimension.");
    break;
  }
  #undef UNROLL_LOOP
}
static void intra_pred_planar_ver_w64(const uvg_pixel* ref, const uvg_pixel* ref_side, const int line, const int shift, __m256i* dst)
{
  const __m256i v_last_ref = _mm256_set1_epi8(ref_side[line + 1]);

  // Got 64 8-bit samples, or 512 bits of data. Load into two vectors
  const __m256i v_ref0 = _mm256_loadu_si256((const __m256i*) &ref[1]);
  const __m256i v_ref1 = _mm256_loadu_si256((const __m256i*) &ref[33]);

  // These stay constant through the loop
  const __m256i v_lo0 = _mm256_unpacklo_epi8(v_ref0, v_last_ref);
  const __m256i v_lo1 = _mm256_unpacklo_epi8(v_ref1, v_last_ref);
  const __m256i v_hi0 = _mm256_unpackhi_epi8(v_ref0, v_last_ref);
  const __m256i v_hi1 = _mm256_unpackhi_epi8(v_ref1, v_last_ref);

  for (uint8_t y = 0, a = line - 1, b = 1, d = 0; y < line; ++y, --a, ++b, d += 4) {
    uint8_t tmp[2] = {a, b};
    uint16_t* tmp2 = (uint16_t*)tmp;
    const __m256i v_ys = _mm256_set1_epi16(*tmp2);
    
    __m256i v_madd_lo0 = _mm256_maddubs_epi16(v_lo0, v_ys);
    __m256i v_madd_lo1 = _mm256_maddubs_epi16(v_lo1, v_ys);
    __m256i v_madd_hi0 = _mm256_maddubs_epi16(v_hi0, v_ys);
    __m256i v_madd_hi1 = _mm256_maddubs_epi16(v_hi1, v_ys);

    dst[d + 0] = _mm256_permute2x128_si256(v_madd_lo0, v_madd_hi0, 0x20);
    dst[d + 1] = _mm256_permute2x128_si256(v_madd_lo0, v_madd_hi0, 0x31);
    dst[d + 2] = _mm256_permute2x128_si256(v_madd_lo1, v_madd_hi1, 0x20);
    dst[d + 3] = _mm256_permute2x128_si256(v_madd_lo1, v_madd_hi1, 0x31);
  }
}


static intra_planar_half_func* planar_func_table[2][7] = {
  {                    NULL,                      NULL,  intra_pred_planar_hor_w4,  intra_pred_planar_hor_w8, intra_pred_planar_hor_w16, intra_pred_planar_hor_w32, intra_pred_planar_hor_w64},
  {                    NULL,                      NULL,  intra_pred_planar_ver_w4,  intra_pred_planar_ver_w8, intra_pred_planar_ver_w16, intra_pred_planar_ver_w32, intra_pred_planar_ver_w64}
};


static void uvg_intra_pred_planar_avx2(const cu_loc_t* const cu_loc,
  color_t color,
  const uint8_t* const ref_top,
  const uint8_t* const ref_left,
  uvg_pixel* dst)
{
  const int16_t width = color == COLOR_Y ? cu_loc->width : cu_loc->chroma_width;
  const int16_t height = color == COLOR_Y ? cu_loc->height : cu_loc->chroma_height;
  const int samples = width * height;
  const __m256i v_samples = _mm256_set1_epi32(samples);

  const int log2_width = uvg_g_convert_to_log2[width];
  const int log2_height = uvg_g_convert_to_log2[height];
  const int shift_r = log2_width + log2_height + 1;
  
  __m256i v_pred_hor[256];
  __m256i v_pred_ver[256];

  intra_planar_half_func* planar_hor = planar_func_table[0][log2_width];
  intra_planar_half_func* planar_ver = planar_func_table[1][log2_width];

  planar_hor(ref_left, ref_top, height, log2_height, v_pred_hor);
  planar_ver(ref_top, ref_left, height, log2_width, v_pred_ver);

  // debug
  int16_t* hor_res = (int16_t*)v_pred_hor;
  int16_t* ver_res = (int16_t*)v_pred_ver;

  // Cast two 16-bit values to 32-bit and fill a 256-bit vector
  int16_t tmp[2] = {height, width};
  int32_t* tmp2 = (int32_t*)tmp;
  const __m256i v_madd_shift = _mm256_set1_epi32(*tmp2);

  __m256i v_res[256];
  // Old loop
  /*for (int i = 0, d = 0; i < samples; i += 16, ++d) {
    v_res[d] = _mm256_add_epi16(v_pred_ver[d], v_pred_hor[d]);
    v_res[d] = _mm256_add_epi16(v_res[d], v_samples);
    v_res[d] = _mm256_srli_epi16(v_res[d], shift_r);
  }*/

  // New loop
  __m128i shift_r_v = _mm_setzero_si128();
  shift_r_v = _mm_insert_epi32(shift_r_v, shift_r, 0);
  for (int i = 0, d = 0; i < samples; i += 16, ++d) {
    __m256i v_lo = _mm256_unpacklo_epi16(v_pred_hor[d], v_pred_ver[d]);
    __m256i v_hi = _mm256_unpackhi_epi16(v_pred_hor[d], v_pred_ver[d]);

    // madd will extend the intermediate results to 32-bit to avoid overflows
    __m256i v_madd_lo = _mm256_madd_epi16(v_lo, v_madd_shift);
    __m256i v_madd_hi = _mm256_madd_epi16(v_hi, v_madd_shift);

    v_madd_lo = _mm256_add_epi32(v_madd_lo, v_samples);
    v_madd_hi = _mm256_add_epi32(v_madd_hi, v_samples);

    v_madd_lo = _mm256_srl_epi32(v_madd_lo, shift_r_v);
    v_madd_hi = _mm256_srl_epi32(v_madd_hi, shift_r_v);

    v_res[d] = _mm256_packs_epi32(v_madd_lo, v_madd_hi);
  }

  // debug
  int16_t* res = (int16_t*)v_res;

  if (samples == 16) {
    __m256i v_tmp = _mm256_packus_epi16(v_res[0], v_res[0]);
    v_tmp = _mm256_permute4x64_epi64(v_tmp, _MM_SHUFFLE(3, 1, 2, 0));
    __m128i v_tmp2 = _mm256_castsi256_si128(v_tmp);
    _mm_store_si128((__m128i*)dst, v_tmp2);
  }
  else {
    for (int i = 0, s = 0; i < samples; i += 32, s += 2) {
      __m256i v_tmp = _mm256_packus_epi16(v_res[s + 0], v_res[s + 1]);
      v_tmp = _mm256_permute4x64_epi64(v_tmp, _MM_SHUFFLE(3, 1, 2, 0));

      _mm256_storeu_si256((__m256i*)&dst[i], v_tmp);
    }
  }
}


// Calculate the DC value for a 4x4 block. The algorithm uses slightly
// different addends, multipliers etc for different pixels in the block,
// but for a fixed-size implementation one vector wide, all the weights,
// addends etc can be preinitialized for each position.
static void pred_filtered_dc_4x4(const uint8_t *ref_top,
                                 const uint8_t *ref_left,
                                       uint8_t *out_block,
                                 const uint8_t multi_ref_idx)
{
  const uint32_t rt_u32 = *(const uint32_t *)(ref_top  + 1);
  const uint32_t rl_u32 = *(const uint32_t *)(ref_left + 1);

  const __m128i zero    = _mm_setzero_si128();
  const __m128i twos    = _mm_set1_epi8(2);

  // Hack. Move 4 u8's to bit positions 0, 64, 128 and 192 in two regs, to
  // expand them to 16 bits sort of "for free". Set highest bits on all the
  // other bytes in vectors to zero those bits in the result vector.
  const __m128i rl_shuf_lo = _mm_setr_epi32(0x80808000, 0x80808080,
                                            0x80808001, 0x80808080);
  const __m128i rl_shuf_hi = _mm_add_epi8  (rl_shuf_lo, twos);

  // Every second multiplier is 1, because we want maddubs to calculate
  // a + bc = 1 * a + bc (actually 2 + bc). We need to fill a vector with
  // ((u8)2)'s for other stuff anyway, so that can also be used here.
  const __m128i mult_lo = _mm_setr_epi32(0x01030102, 0x01030103,
                                         0x01040103, 0x01040104);
  const __m128i mult_hi = _mm_setr_epi32(0x01040103, 0x01040104,
                                         0x01040103, 0x01040104);
  __m128i four         = _mm_cvtsi32_si128  (4);
  __m128i rt           = _mm_cvtsi32_si128  (rt_u32);
  __m128i rl           = _mm_cvtsi32_si128  (rl_u32);
  __m128i rtrl         = _mm_unpacklo_epi32 (rt, rl);

  __m128i sad0         = _mm_sad_epu8       (rtrl, zero);
  __m128i sad1         = _mm_shuffle_epi32  (sad0, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad2         = _mm_add_epi64      (sad0, sad1);
  __m128i sad3         = _mm_add_epi64      (sad2, four);

  __m128i dc_64        = _mm_srli_epi64     (sad3, 3);
  __m128i dc_8         = _mm_broadcastb_epi8(dc_64);

  __m128i rl_lo        = _mm_shuffle_epi8   (rl, rl_shuf_lo);
  __m128i rl_hi        = _mm_shuffle_epi8   (rl, rl_shuf_hi);

  __m128i rt_lo        = _mm_unpacklo_epi8  (rt, zero);
  __m128i rt_hi        = zero;

  __m128i dc_addend    = _mm_unpacklo_epi8(dc_8, twos);

  __m128i dc_multd_lo  = _mm_maddubs_epi16(dc_addend,    mult_lo);
  __m128i dc_multd_hi  = _mm_maddubs_epi16(dc_addend,    mult_hi);

  __m128i rl_rt_lo     = _mm_add_epi16    (rl_lo,        rt_lo);
  __m128i rl_rt_hi     = _mm_add_epi16    (rl_hi,        rt_hi);

  __m128i res_lo       = _mm_add_epi16    (dc_multd_lo,  rl_rt_lo);
  __m128i res_hi       = _mm_add_epi16    (dc_multd_hi,  rl_rt_hi);

          res_lo       = _mm_srli_epi16   (res_lo,       2);
          res_hi       = _mm_srli_epi16   (res_hi,       2);

  __m128i final        = _mm_packus_epi16 (res_lo,       res_hi);
  _mm_storeu_si128((__m128i *)out_block, final);
}

static void pred_filtered_dc_8x8(const uint8_t *ref_top,
                                 const uint8_t *ref_left,
                                       uint8_t *out_block,
                                 const uint8_t multi_ref_idx)
{
  const uint64_t rt_u64 = *(const uint64_t *)(ref_top  + 1);
  const uint64_t rl_u64 = *(const uint64_t *)(ref_left + 1);

  const __m128i zero128 = _mm_setzero_si128();
  const __m256i twos    = _mm256_set1_epi8(2);

  // DC multiplier is 2 at (0, 0), 3 at (*, 0) and (0, *), and 4 at (*, *).
  // There is a constant addend of 2 on each pixel, use values from the twos
  // register and multipliers of 1 for that, to use maddubs for an (a*b)+c
  // operation.
  const __m256i mult_up_lo = _mm256_setr_epi32(0x01030102, 0x01030103,
                                               0x01030103, 0x01030103,
                                               0x01040103, 0x01040104,
                                               0x01040104, 0x01040104);

  // The 6 lowest rows have same multipliers, also the DC values and addends
  // are the same so this works for all of those
  const __m256i mult_rest  = _mm256_permute4x64_epi64(mult_up_lo, _MM_SHUFFLE(3, 2, 3, 2));

  // Every 8-pixel row starts with the next pixel of ref_left. Along with
  // doing the shuffling, also expand u8->u16, ie. move bytes 0 and 1 from
  // ref_left to bit positions 0 and 128 in rl_up_lo, 2 and 3 to rl_up_hi,
  // etc. The places to be zeroed out are 0x80 instead of the usual 0xff,
  // because this allows us to form new masks on the fly by adding 0x02-bytes
  // to this mask and still retain the highest bits as 1 where things should
  // be zeroed out.
  const __m256i rl_shuf_up_lo = _mm256_setr_epi32(0x80808000, 0x80808080,
                                                  0x80808080, 0x80808080,
                                                  0x80808001, 0x80808080,
                                                  0x80808080, 0x80808080);
  // And don't waste memory or architectural regs, hope these instructions
  // will be placed in between the shuffles by the compiler to only use one
  // register for the shufmasks, and executed way ahead of time because their
  // regs can be renamed.
  const __m256i rl_shuf_up_hi = _mm256_add_epi8 (rl_shuf_up_lo, twos);
  const __m256i rl_shuf_dn_lo = _mm256_add_epi8 (rl_shuf_up_hi, twos);
  const __m256i rl_shuf_dn_hi = _mm256_add_epi8 (rl_shuf_dn_lo, twos);

  __m128i eight         = _mm_cvtsi32_si128     (8);
  __m128i rt            = _mm_cvtsi64_si128     (rt_u64);
  __m128i rl            = _mm_cvtsi64_si128     (rl_u64);
  __m128i rtrl          = _mm_unpacklo_epi64    (rt, rl);

  __m128i sad0          = _mm_sad_epu8          (rtrl, zero128);
  __m128i sad1          = _mm_shuffle_epi32     (sad0, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad2          = _mm_add_epi64         (sad0, sad1);
  __m128i sad3          = _mm_add_epi64         (sad2, eight);

  __m128i dc_64         = _mm_srli_epi64        (sad3, 4);
  __m256i dc_8          = _mm256_broadcastb_epi8(dc_64);

  __m256i dc_addend     = _mm256_unpacklo_epi8  (dc_8, twos);

  __m256i dc_up_lo      = _mm256_maddubs_epi16  (dc_addend, mult_up_lo);
  __m256i dc_rest       = _mm256_maddubs_epi16  (dc_addend, mult_rest);

  // rt_dn is all zeros, as is rt_up_hi. This'll get us the rl and rt parts
  // in A|B, C|D order instead of A|C, B|D that could be packed into abcd
  // order, so these need to be permuted before adding to the weighed DC
  // values.
  __m256i rt_up_lo      = _mm256_cvtepu8_epi16   (rt);

  __m256i rlrlrlrl      = _mm256_broadcastq_epi64(rl);
  __m256i rl_up_lo      = _mm256_shuffle_epi8    (rlrlrlrl, rl_shuf_up_lo);

  // Everything ref_top is zero except on the very first row
  __m256i rt_rl_up_hi   = _mm256_shuffle_epi8    (rlrlrlrl, rl_shuf_up_hi);
  __m256i rt_rl_dn_lo   = _mm256_shuffle_epi8    (rlrlrlrl, rl_shuf_dn_lo);
  __m256i rt_rl_dn_hi   = _mm256_shuffle_epi8    (rlrlrlrl, rl_shuf_dn_hi);

  __m256i rt_rl_up_lo   = _mm256_add_epi16       (rt_up_lo, rl_up_lo);

  __m256i rt_rl_up_lo_2 = _mm256_permute2x128_si256(rt_rl_up_lo, rt_rl_up_hi, 0x20);
  __m256i rt_rl_up_hi_2 = _mm256_permute2x128_si256(rt_rl_up_lo, rt_rl_up_hi, 0x31);
  __m256i rt_rl_dn_lo_2 = _mm256_permute2x128_si256(rt_rl_dn_lo, rt_rl_dn_hi, 0x20);
  __m256i rt_rl_dn_hi_2 = _mm256_permute2x128_si256(rt_rl_dn_lo, rt_rl_dn_hi, 0x31);

  __m256i up_lo = _mm256_add_epi16(rt_rl_up_lo_2, dc_up_lo);
  __m256i up_hi = _mm256_add_epi16(rt_rl_up_hi_2, dc_rest);
  __m256i dn_lo = _mm256_add_epi16(rt_rl_dn_lo_2, dc_rest);
  __m256i dn_hi = _mm256_add_epi16(rt_rl_dn_hi_2, dc_rest);

          up_lo = _mm256_srli_epi16(up_lo, 2);
          up_hi = _mm256_srli_epi16(up_hi, 2);
          dn_lo = _mm256_srli_epi16(dn_lo, 2);
          dn_hi = _mm256_srli_epi16(dn_hi, 2);

  __m256i res_up = _mm256_packus_epi16(up_lo, up_hi);
  __m256i res_dn = _mm256_packus_epi16(dn_lo, dn_hi);

  _mm256_storeu_si256(((__m256i *)out_block) + 0, res_up);
  _mm256_storeu_si256(((__m256i *)out_block) + 1, res_dn);
}

static INLINE __m256i cvt_u32_si256(const uint32_t u)
{
  const __m256i zero = _mm256_setzero_si256();
  return _mm256_insert_epi32(zero, u, 0);
}

static void pred_filtered_dc_16x16(const uint8_t *ref_top,
                                   const uint8_t *ref_left,
                                         uint8_t *out_block,
                                   const uint8_t multi_ref_idx)
{
  const __m128i rt_128 = _mm_loadu_si128((const __m128i *)(ref_top  + 1));
  const __m128i rl_128 = _mm_loadu_si128((const __m128i *)(ref_left + 1));

  const __m128i zero_128 = _mm_setzero_si128();
  const __m256i zero     = _mm256_setzero_si256();
  const __m256i twos     = _mm256_set1_epi8(2);

  const __m256i mult_r0  = _mm256_setr_epi32(0x01030102, 0x01030103,
                                             0x01030103, 0x01030103,
                                             0x01030103, 0x01030103,
                                             0x01030103, 0x01030103);

  const __m256i mult_left = _mm256_set1_epi16(0x0103);

  // Leftmost bytes' blend mask, to move bytes (pixels) from the leftmost
  // column vector to the result row
  const __m256i lm8_bmask = _mm256_setr_epi32(0xff, 0, 0, 0, 0xff, 0, 0, 0);

  __m128i sixteen = _mm_cvtsi32_si128(16);
  __m128i sad0_t  = _mm_sad_epu8 (rt_128, zero_128);
  __m128i sad0_l  = _mm_sad_epu8 (rl_128, zero_128);
  __m128i sad0    = _mm_add_epi64(sad0_t, sad0_l);

  __m128i sad1    = _mm_shuffle_epi32      (sad0, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad2    = _mm_add_epi64          (sad0, sad1);
  __m128i sad3    = _mm_add_epi64          (sad2, sixteen);

  __m128i dc_64   = _mm_srli_epi64         (sad3, 5);
  __m256i dc_8    = _mm256_broadcastb_epi8 (dc_64);

  __m256i rt      = _mm256_cvtepu8_epi16   (rt_128);
  __m256i rl      = _mm256_cvtepu8_epi16   (rl_128);

  uint8_t rl0       = *(uint8_t *)(ref_left + 1);
  __m256i rl_r0     = cvt_u32_si256((uint32_t)rl0);

  __m256i rlrt_r0   = _mm256_add_epi16(rl_r0, rt);

  __m256i dc_addend = _mm256_unpacklo_epi8(dc_8, twos);
  __m256i r0        = _mm256_maddubs_epi16(dc_addend, mult_r0);
  __m256i left_dcs  = _mm256_maddubs_epi16(dc_addend, mult_left);

          r0        = _mm256_add_epi16    (r0,       rlrt_r0);
          r0        = _mm256_srli_epi16   (r0, 2);
  __m256i r0r0      = _mm256_packus_epi16 (r0, r0);
          r0r0      = _mm256_permute4x64_epi64(r0r0, _MM_SHUFFLE(3, 1, 2, 0));

  __m256i leftmosts = _mm256_add_epi16    (left_dcs,  rl);
          leftmosts = _mm256_srli_epi16   (leftmosts, 2);

  // Contain the leftmost column's bytes in both lanes of lm_8
  __m256i lm_8      = _mm256_packus_epi16 (leftmosts, zero);
          lm_8      = _mm256_permute4x64_epi64(lm_8,  _MM_SHUFFLE(2, 0, 2, 0));

  __m256i lm8_r1    = _mm256_srli_epi32       (lm_8, 8);
  __m256i r1r1      = _mm256_blendv_epi8      (dc_8, lm8_r1, lm8_bmask);
  __m256i r0r1      = _mm256_blend_epi32      (r0r0, r1r1, 0xf0);

  _mm256_storeu_si256((__m256i *)out_block, r0r1);

  // Starts from 2 because row 0 (and row 1) is handled separately
  __m256i lm8_l     = _mm256_bsrli_epi128     (lm_8, 2);
  __m256i lm8_h     = _mm256_bsrli_epi128     (lm_8, 3);
          lm_8      = _mm256_blend_epi32      (lm8_l, lm8_h, 0xf0);

  for (uint32_t y = 2; y < 16; y += 2) {
    __m256i curr_row = _mm256_blendv_epi8 (dc_8, lm_8, lm8_bmask);
    _mm256_storeu_si256((__m256i *)(out_block + (y << 4)), curr_row);
    lm_8 = _mm256_bsrli_epi128(lm_8, 2);
  }
}

static void pred_filtered_dc_32x32(const uint8_t *ref_top,
                                   const uint8_t *ref_left,
                                         uint8_t *out_block,
                                   const uint8_t multi_ref_idx)
{
  const __m256i rt = _mm256_loadu_si256((const __m256i *)(ref_top  + 1));
  const __m256i rl = _mm256_loadu_si256((const __m256i *)(ref_left + 1));

  const __m256i zero = _mm256_setzero_si256();
  const __m256i twos = _mm256_set1_epi8(2);

  const __m256i mult_r0lo = _mm256_setr_epi32(0x01030102, 0x01030103,
                                              0x01030103, 0x01030103,
                                              0x01030103, 0x01030103,
                                              0x01030103, 0x01030103);

  const __m256i mult_left = _mm256_set1_epi16(0x0103);
  const __m256i lm8_bmask = cvt_u32_si256    (0xff);

  const __m256i bshif_msk = _mm256_setr_epi32(0x04030201, 0x08070605,
                                              0x0c0b0a09, 0x800f0e0d,
                                              0x03020100, 0x07060504,
                                              0x0b0a0908, 0x0f0e0d0c);
  __m256i debias = cvt_u32_si256(32);
  __m256i sad0_t = _mm256_sad_epu8         (rt,     zero);
  __m256i sad0_l = _mm256_sad_epu8         (rl,     zero);
  __m256i sad0   = _mm256_add_epi64        (sad0_t, sad0_l);

  __m256i sad1   = _mm256_permute4x64_epi64(sad0,   _MM_SHUFFLE(1, 0, 3, 2));
  __m256i sad2   = _mm256_add_epi64        (sad0,   sad1);
  __m256i sad3   = _mm256_shuffle_epi32    (sad2,   _MM_SHUFFLE(1, 0, 3, 2));
  __m256i sad4   = _mm256_add_epi64        (sad2,   sad3);
  __m256i sad5   = _mm256_add_epi64        (sad4,   debias);
  __m256i dc_64  = _mm256_srli_epi64       (sad5,   6);

  __m128i dc_64_ = _mm256_castsi256_si128  (dc_64);
  __m256i dc_8   = _mm256_broadcastb_epi8  (dc_64_);

  __m256i rtlo   = _mm256_unpacklo_epi8    (rt, zero);
  __m256i rllo   = _mm256_unpacklo_epi8    (rl, zero);
  __m256i rthi   = _mm256_unpackhi_epi8    (rt, zero);
  __m256i rlhi   = _mm256_unpackhi_epi8    (rl, zero);

  __m256i dc_addend = _mm256_unpacklo_epi8 (dc_8, twos);
  __m256i r0lo   = _mm256_maddubs_epi16    (dc_addend, mult_r0lo);
  __m256i r0hi   = _mm256_maddubs_epi16    (dc_addend, mult_left);
  __m256i c0dc   = r0hi;

          r0lo   = _mm256_add_epi16        (r0lo, rtlo);
          r0hi   = _mm256_add_epi16        (r0hi, rthi);

  __m256i rlr0   = _mm256_blendv_epi8      (zero, rl, lm8_bmask);
          r0lo   = _mm256_add_epi16        (r0lo, rlr0);

          r0lo   = _mm256_srli_epi16       (r0lo, 2);
          r0hi   = _mm256_srli_epi16       (r0hi, 2);
  __m256i r0     = _mm256_packus_epi16     (r0lo, r0hi);

  _mm256_storeu_si256((__m256i *)out_block, r0);

  __m256i c0lo   = _mm256_add_epi16        (c0dc, rllo);
  __m256i c0hi   = _mm256_add_epi16        (c0dc, rlhi);
          c0lo   = _mm256_srli_epi16       (c0lo, 2);
          c0hi   = _mm256_srli_epi16       (c0hi, 2);

  __m256i c0     = _mm256_packus_epi16     (c0lo, c0hi);

  // r0 already handled!
  for (uint32_t y = 1; y < 32; y++) {
    if (y == 16) {
      c0 = _mm256_permute4x64_epi64(c0, _MM_SHUFFLE(1, 0, 3, 2));
    } else {
      c0 = _mm256_shuffle_epi8     (c0, bshif_msk);
    }
    __m256i curr_row = _mm256_blendv_epi8 (dc_8, c0, lm8_bmask);
    _mm256_storeu_si256(((__m256i *)out_block) + y, curr_row);
  }
}

/**
* \brief Generage intra DC prediction with post filtering applied.
* \param log2_width    Log2 of width, range 2..5.
* \param in_ref_above  Pointer to -1 index of above reference, length=width*2+1.
* \param in_ref_left   Pointer to -1 index of left reference, length=width*2+1.
* \param dst           Buffer of size width*width.
* \param multi_ref_idx Reference line index. May be non-zero when MRL is used.
*/
static void uvg_intra_pred_filtered_dc_avx2(
  const int_fast8_t log2_width,
  const uint8_t *ref_top,
  const uint8_t *ref_left,
        uint8_t *out_block,
  const uint8_t multi_ref_idx)
{
  assert(log2_width >= 2 && log2_width <= 5);

  // TODO: implement multi reference index for all subfunctions
  if (log2_width == 2) {
    pred_filtered_dc_4x4(ref_top, ref_left, out_block, multi_ref_idx);
  } else if (log2_width == 3) {
    pred_filtered_dc_8x8(ref_top, ref_left, out_block, multi_ref_idx);
  } else if (log2_width == 4) {
    pred_filtered_dc_16x16(ref_top, ref_left, out_block, multi_ref_idx);
  } else if (log2_width == 5) {
    pred_filtered_dc_32x32(ref_top, ref_left, out_block, multi_ref_idx);
  }
}

// TODO: update all ranges (in comments, etc.) from HEVC to VVC

/**
* \brief Position Dependent Prediction Combination for Planar and DC modes.
* \param log2_width    Log2 of width, range 2..5.
* \param width         Block width matching log2_width.
* \param used_ref      Pointer used reference pixel struct.
* \param dst           Buffer of size width*width.
*/
// TODO: does not work with blocks with height 1 and 2
// TODO: also has width someplaces where height should be
static void uvg_pdpc_planar_dc_avx2(
  const int mode,
  const cu_loc_t* const cu_loc,
  const color_t color,
  const uvg_intra_ref *const used_ref,
  uvg_pixel *const dst)
{
  // ISP_TODO: non-square block implementation, height is passed but not used
  assert(mode == 0 || mode == 1);  // planar or DC
  const int width = color == COLOR_Y ? cu_loc->width : cu_loc->chroma_width;
  const int height = color == COLOR_Y ? cu_loc->height : cu_loc->chroma_height;
  const int log2_width =  uvg_g_convert_to_log2[width];
  const int log2_height = uvg_g_convert_to_log2[height];

  __m256i shuf_mask_byte = _mm256_setr_epi8(
    0, -1, 0, -1, 0, -1, 0, -1,
    1, -1, 1, -1, 1, -1, 1, -1,
    2, -1, 2, -1, 2, -1, 2, -1,
    3, -1, 3, -1, 3, -1, 3, -1
  );

  __m256i shuf_mask_word = _mm256_setr_epi8(
    0, 1, 0, 1, 0, 1, 0, 1,
    2, 3, 2, 3, 2, 3, 2, 3,
    4, 5, 4, 5, 4, 5, 4, 5,
    6, 7, 6, 7, 6, 7, 6, 7
  );

  // TODO: replace latter log2_width with log2_height
  const int scale = ((log2_width - 2 + log2_width - 2 + 2) >> 2);

  // Same weights regardless of axis, compute once
  int16_t w[LCU_WIDTH];
  for (int i = 0; i < width; i += 4) {
    __m128i base = _mm_set1_epi32(i);
    __m128i offs = _mm_setr_epi32(0, 1, 2, 3);
    __m128i idxs = _mm_add_epi32(base, offs);
    __m128i unclipped = _mm_slli_epi32(idxs, 1);
    unclipped = _mm_srli_epi32(unclipped, scale);
    __m128i clipped = _mm_min_epi32( _mm_set1_epi32(31), unclipped);
    __m128i weights = _mm_srlv_epi32(_mm_set1_epi32(32), clipped);
    weights = _mm_packus_epi32(weights, weights);
    _mm_storel_epi64((__m128i*)&w[i], weights);
  }

  // Process in 4x4 blocks
  for (int y = 0; y < height; y += 4) {
    for (int x = 0; x < width; x += 4) {

      uint32_t dw_left;
      uint32_t dw_top;
      memcpy(&dw_left, &used_ref->left[y + 1], sizeof(dw_left));
      memcpy(&dw_top , &used_ref->top [x + 1], sizeof(dw_top));
      __m256i vleft = _mm256_set1_epi32(dw_left);
      __m256i vtop  = _mm256_set1_epi32(dw_top);
      vleft  = _mm256_shuffle_epi8(vleft, shuf_mask_byte);
      vtop = _mm256_cvtepu8_epi16(_mm256_castsi256_si128(vtop));

      __m128i vseq = _mm_setr_epi32(0, 1, 2, 3);
      __m128i vidx = _mm_slli_epi32(vseq, log2_width);
      __m128i vdst = _mm_i32gather_epi32((const int32_t*)(dst + y * width + x), vidx, 1);
      __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);
      uint64_t quad_wL;
      uint64_t quad_wT;
      memcpy(&quad_wL, &w[x], sizeof(quad_wL));
      memcpy(&quad_wT, &w[y], sizeof(quad_wT));
      __m256i vwL = _mm256_set1_epi64x(quad_wL); 
      __m256i vwT = _mm256_set1_epi64x(quad_wT);
      vwT = _mm256_shuffle_epi8(vwT, shuf_mask_word);
      __m256i diff_left = _mm256_sub_epi16(vleft, vdst16);
      __m256i diff_top  = _mm256_sub_epi16(vtop , vdst16);
      __m256i prod_left = _mm256_mullo_epi16(vwL, diff_left);
      __m256i prod_top  = _mm256_mullo_epi16(vwT, diff_top);
      __m256i accu = _mm256_add_epi16(prod_left, prod_top);
      accu = _mm256_add_epi16(accu, _mm256_set1_epi16(32));
      accu = _mm256_srai_epi16(accu, 6);
      accu = _mm256_add_epi16(vdst16, accu);

      __m128i lo = _mm256_castsi256_si128(accu);
      __m128i hi = _mm256_extracti128_si256(accu, 1);
      vdst = _mm_packus_epi16(lo, hi);

      *(uint32_t*)(dst + (y + 0) * width + x) = _mm_extract_epi32(vdst, 0);
      *(uint32_t*)(dst + (y + 1) * width + x) = _mm_extract_epi32(vdst, 1);
      *(uint32_t*)(dst + (y + 2) * width + x) = _mm_extract_epi32(vdst, 2);
      *(uint32_t*)(dst + (y + 3) * width + x) = _mm_extract_epi32(vdst, 3);
    }
  }
}

static INLINE void mip_ref_downsampling_4x4_4to2_avx2(uvg_pixel* reduced_dst, const uvg_pixel* const ref_top, const uvg_pixel* const ref_left)
{
  const uint8_t down_smp_factor = 2; // width / red_bdry_size
  const int log2_factor = uvg_g_convert_to_log2[down_smp_factor];
  const int rounding_offset = (1 << (log2_factor - 1));

  const __m128i vrnd = _mm_set1_epi16(rounding_offset);

  ALIGNED(16) uint32_t ref[2];
  ref[0] = *(uint32_t*)ref_top;
  ref[1] = *(uint32_t*)ref_left;

  __m128i vref = _mm_load_si128((__m128i*)ref);
  vref = _mm_cvtepu8_epi16(vref);

  __m128i vres = _mm_hadd_epi16(vref, vref);

  vres = _mm_add_epi16(vres, vrnd);
  vres = _mm_srli_epi16(vres, log2_factor);
  __m128i vout = _mm_packus_epi16(vres, vres);

  *(uint32_t*)reduced_dst = _mm_extract_epi32(vout, 0);
}

static INLINE void mip_ref_downsampling_8x8_8to4_avx2(uvg_pixel* reduced_dst, const uvg_pixel* const ref_top, const uvg_pixel* const ref_left)
{
  const uint8_t down_smp_factor = 2; // width / red_bdry_size
  const int log2_factor = uvg_g_convert_to_log2[down_smp_factor];
  const int rounding_offset = (1 << (log2_factor - 1));

  const __m256i vrnd = _mm256_set1_epi16(rounding_offset);

  ALIGNED(16) uint64_t ref[2];
  ref[0] = *(uint64_t*)ref_top;
  ref[1] = *(uint64_t*)ref_left;

  __m128i vref = _mm_load_si128((__m128i*)ref);
  __m256i vref256 = _mm256_cvtepu8_epi16(vref);

  __m256i vres = _mm256_hadd_epi16(vref256, vref256);
  vres = _mm256_permute4x64_epi64(vres, _MM_SHUFFLE(3, 1, 2, 0));

  vres = _mm256_add_epi16(vres, vrnd);
  vres = _mm256_srli_epi16(vres, log2_factor);
  __m256i vout = _mm256_packus_epi16(vres, vres);

  *(uint64_t*)reduced_dst = _mm256_extract_epi64(vout, 0);
}

static INLINE void mip_ref_downsampling_1D_8to4_avx2(uvg_pixel* reduced_dst, const uvg_pixel* const ref_src)
{
  const uint8_t down_smp_factor = 2; // width / red_bdry_size
  const int log2_factor = uvg_g_convert_to_log2[down_smp_factor];
  const int rounding_offset = (1 << (log2_factor - 1));

  const __m128i vrnd = _mm_set1_epi16(rounding_offset);

  __m128i vref = _mm_loadu_si128((__m128i*)ref_src); // Half the data is garbage and will be ignored.
  vref = _mm_cvtepu8_epi16(vref);
  __m128i vres = _mm_hadd_epi16(vref, vref);
  vres = _mm_add_epi16(vres, vrnd);
  vres = _mm_srli_epi16(vres, log2_factor);
  __m128i vout = _mm_packus_epi16(vres, vres);

  *(int32_t*)reduced_dst = _mm_extract_epi32(vout, 0);
}

static INLINE void mip_ref_downsampling_1D_16to4_avx2(uvg_pixel* reduced_dst, const uvg_pixel* const ref_src)
{
  const uint8_t down_smp_factor = 4; // width / red_bdry_size
  const int log2_factor = uvg_g_convert_to_log2[down_smp_factor];
  const int rounding_offset = (1 << (log2_factor - 1));

  const __m256i vrnd = _mm256_set1_epi16(rounding_offset);

  // TODO: try _mm256_dpbuud. 
  // NOTE: ignore this TODO for now, using dpbuud causes error 0xC000001D: Illegal Instruction.
  // The instruction requires a newer CPU.
  __m128i vref = _mm_loadu_si128((__m128i*)ref_src);
  __m256i vref256 = _mm256_cvtepu8_epi16(vref);
  __m256i vres = _mm256_hadd_epi16(vref256, vref256);
  vres = _mm256_permute4x64_epi64(vres, _MM_SHUFFLE(3, 1, 2, 0));
  vres = _mm256_hadd_epi16(vres, vres);
  vres = _mm256_add_epi16(vres, vrnd);
  vres = _mm256_srli_epi16(vres, log2_factor);
  __m256i vout = _mm256_packus_epi16(vres, vres);

  *(int32_t*)(reduced_dst + 0) = _mm256_extract_epi32(vout, 0);
}

static INLINE void mip_ref_downsampling_1D_32to4_avx2(uvg_pixel* reduced_dst, const uvg_pixel* const ref_src)
{
  const uint8_t down_smp_factor = 8; // width / red_bdry_size
  const int log2_factor = uvg_g_convert_to_log2[down_smp_factor];
  const int rounding_offset = (1 << (log2_factor - 1));

  const __m256i vrnd = _mm256_set1_epi16(rounding_offset);

  __m128i vrefa = _mm_loadu_si128((__m128i*)(ref_src + 0));
  __m128i vrefb = _mm_loadu_si128((__m128i*)(ref_src + 16));

  __m256i vref256a = _mm256_cvtepu8_epi16(vrefa);
  __m256i vref256b = _mm256_cvtepu8_epi16(vrefb);

  __m256i vres = _mm256_hadd_epi16(vref256a, vref256b);
  vres = _mm256_permute4x64_epi64(vres, _MM_SHUFFLE(3, 1, 2, 0));
  vres = _mm256_hadd_epi16(vres, vres);
  vres = _mm256_permute4x64_epi64(vres, _MM_SHUFFLE(3, 1, 2, 0));
  vres = _mm256_hadd_epi16(vres, vres);

  vres = _mm256_add_epi16(vres, vrnd);
  vres = _mm256_srli_epi16(vres, log2_factor);
  __m256i vout = _mm256_packus_epi16(vres, vres);

  *(int32_t*)(reduced_dst + 0) = _mm256_extract_epi32(vout, 0);
}

static INLINE void mip_ref_downsampling_1D_64to4_avx2(uvg_pixel* reduced_dst, const uvg_pixel* const ref_src)
{
  const uint8_t down_smp_factor = 16; // width / red_bdry_size
  const int log2_factor = uvg_g_convert_to_log2[down_smp_factor];
  const int rounding_offset = (1 << (log2_factor - 1));

  const __m256i vrnd = _mm256_set1_epi16(rounding_offset);

  __m128i vrefa = _mm_loadu_si128((__m128i*)(ref_src + 0));
  __m128i vrefb = _mm_loadu_si128((__m128i*)(ref_src + 16));
  __m128i vrefc = _mm_loadu_si128((__m128i*)(ref_src + 32));
  __m128i vrefd = _mm_loadu_si128((__m128i*)(ref_src + 48));

  __m256i vref256a = _mm256_cvtepu8_epi16(vrefa);
  __m256i vref256b = _mm256_cvtepu8_epi16(vrefb);
  __m256i vref256c = _mm256_cvtepu8_epi16(vrefc);
  __m256i vref256d = _mm256_cvtepu8_epi16(vrefd);


  __m256i vres0 = _mm256_hadd_epi16(vref256a, vref256b);
  __m256i vres1 = _mm256_hadd_epi16(vref256c, vref256d);
  vres0 = _mm256_permute4x64_epi64(vres0, _MM_SHUFFLE(3, 1, 2, 0));
  vres1 = _mm256_permute4x64_epi64(vres1, _MM_SHUFFLE(3, 1, 2, 0));

  vres0 = _mm256_hadd_epi16(vres0, vres1);
  vres0 = _mm256_permute4x64_epi64(vres0, _MM_SHUFFLE(3, 1, 2, 0));

  vres0 = _mm256_hadd_epi16(vres0, vres0);
  vres0 = _mm256_permute4x64_epi64(vres0, _MM_SHUFFLE(3, 1, 2, 0));

  vres0 = _mm256_hadd_epi16(vres0, vres0);

  vres0 = _mm256_add_epi16(vres0, vrnd);
  vres0 = _mm256_srli_epi16(vres0, log2_factor);
  __m256i vout = _mm256_packus_epi16(vres0, vres0);

  *(int32_t*)(reduced_dst + 0) = _mm256_extract_epi32(vout, 0);
  //*(int32_t*)(reduced_dst + 2) = _mm_extract_epi16(vout, 8);
}


// Size ID 0
static INLINE void mip_reduced_pred_sid0_avx2(uvg_pixel* const output,
  const int16_t* const input,
  const uint16_t* matrix,
  const bool transpose,
  const int in_offset,
  const int in_offset_tr)
{
  const int input_size = 4;
  // const int pred_size = 4;
  // const int size_id = 0;

  int sum = 0;
  for (int i = 0; i < input_size; i++) {
    sum += input[i];
  }
  const int offset = (1 << (MIP_SHIFT_MATRIX - 1)) - MIP_OFFSET_MATRIX * sum;

  const __m128i vofs = _mm_set1_epi32(offset);

  const uint16_t* weight = matrix;
  const int input_offset = transpose ? in_offset_tr : in_offset;

  const __m128i vinofs = _mm_set1_epi32(input_offset);

  const __m128i vshuf = _mm_setr_epi8(
    0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
    0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07
  );

  const __m128i vtranspose = _mm_setr_epi8(
    0x00, 0x04, 0x08, 0x0c, 0x01, 0x05, 0x09, 0x0d,
    0x02, 0x06, 0x0a, 0x0e, 0x03, 0x07, 0x0b, 0x0f
  );
  
  const __m128i vinraw = _mm_loadu_si128((__m128i*)input);
  const __m128i vin = _mm_shuffle_epi8(vinraw, vshuf);

  // Calculate first half
  __m128i vweight0 = _mm_loadu_si128((__m128i*) &weight[0]);
  __m128i vweight1 = _mm_loadu_si128((__m128i*) &weight[8]);
  __m128i vweight2 = _mm_loadu_si128((__m128i*) &weight[16]);
  __m128i vweight3 = _mm_loadu_si128((__m128i*) &weight[24]);

  weight += 32;

  __m128i vmadd0 = _mm_madd_epi16(vin, vweight0);
  __m128i vmadd1 = _mm_madd_epi16(vin, vweight1);
  __m128i vmadd2 = _mm_madd_epi16(vin, vweight2);
  __m128i vmadd3 = _mm_madd_epi16(vin, vweight3);

  __m128i vresult0 = _mm_hadd_epi32(vmadd0, vmadd1);
  __m128i vresult1 = _mm_hadd_epi32(vmadd2, vmadd3);

  vresult0 = _mm_add_epi32(vresult0, vofs);
  vresult0 = _mm_srai_epi32(vresult0, MIP_SHIFT_MATRIX);
  vresult0 = _mm_add_epi32(vresult0, vinofs);

  vresult1 = _mm_add_epi32(vresult1, vofs);
  vresult1 = _mm_srai_epi32(vresult1, MIP_SHIFT_MATRIX);
  vresult1 = _mm_add_epi32(vresult1, vinofs);

  __m128i vres16_a = _mm_packus_epi32(vresult0, vresult1);

  // Calculate second half
  vweight0 = _mm_loadu_si128((__m128i*) & weight[0]);
  vweight1 = _mm_loadu_si128((__m128i*) & weight[8]);
  vweight2 = _mm_loadu_si128((__m128i*) & weight[16]);
  vweight3 = _mm_loadu_si128((__m128i*) & weight[24]);

  vmadd0 = _mm_madd_epi16(vin, vweight0);
  vmadd1 = _mm_madd_epi16(vin, vweight1);
  vmadd2 = _mm_madd_epi16(vin, vweight2);
  vmadd3 = _mm_madd_epi16(vin, vweight3);

  vresult0 = _mm_hadd_epi32(vmadd0, vmadd1);
  vresult1 = _mm_hadd_epi32(vmadd2, vmadd3);

  vresult0 = _mm_add_epi32(vresult0, vofs);
  vresult0 = _mm_srai_epi32(vresult0, MIP_SHIFT_MATRIX);
  vresult0 = _mm_add_epi32(vresult0, vinofs);

  vresult1 = _mm_add_epi32(vresult1, vofs);
  vresult1 = _mm_srai_epi32(vresult1, MIP_SHIFT_MATRIX);
  vresult1 = _mm_add_epi32(vresult1, vinofs);

  __m128i vres16_b = _mm_packus_epi32(vresult0, vresult1);
  __m128i vres8 = _mm_packus_epi16(vres16_a, vres16_b);

  if (transpose) {
    vres8 = _mm_shuffle_epi8(vres8, vtranspose);
    _mm_storeu_si128((__m128i*)output, vres8);
  }
  else {
    _mm_storeu_si128((__m128i*)output, vres8);
  }
}

// Size ID 1
static void INLINE mip_reduced_pred_sid1_avx2(uvg_pixel* const output,
  const int16_t* const input,
  const uint16_t* matrix,
  const bool transpose,
  const int in_offset,
  const int in_offset_tr)
{
  const int input_size = 8;
  const int pred_size = 4;
  const int size_id = 1;

  int sum = 0;
  for (int i = 0; i < input_size; i++) {
    sum += input[i];
  }
  const int offset = (1 << (MIP_SHIFT_MATRIX - 1)) - MIP_OFFSET_MATRIX * sum;

  const __m128i vofs = _mm_set1_epi32(offset);

  const uint16_t* weight = matrix;
  const int input_offset = transpose ? in_offset_tr : in_offset;

  const __m128i vinofs = _mm_set1_epi32(input_offset);

  const __m128i vshuf0 = _mm_setr_epi8(
    0x00, 0x01, 0x02, 0x03, 0x00, 0x01, 0x02, 0x03,
    0x00, 0x01, 0x02, 0x03, 0x00, 0x01, 0x02, 0x03);
  const __m128i vshuf1 = _mm_setr_epi8(
    0x04, 0x05, 0x06, 0x07, 0x04, 0x05, 0x06, 0x07,
    0x04, 0x05, 0x06, 0x07, 0x04, 0x05, 0x06, 0x07);
  const __m128i vshuf2 = _mm_setr_epi8(
    0x08, 0x09, 0x0a, 0x0b, 0x08, 0x09, 0x0a, 0x0b,
    0x08, 0x09, 0x0a, 0x0b, 0x08, 0x09, 0x0a, 0x0b);
  const __m128i vshuf3 = _mm_setr_epi8(
    0x0c, 0x0d, 0x0e, 0x0f, 0x0c, 0x0d, 0x0e, 0x0f,
    0x0c, 0x0d, 0x0e, 0x0f, 0x0c, 0x0d, 0x0e, 0x0f);
  const __m128i vtranspose = _mm_setr_epi8(
    0x00, 0x04, 0x08, 0x0c, 0x01, 0x05, 0x09, 0x0d,
    0x02, 0x06, 0x0a, 0x0e, 0x03, 0x07, 0x0b, 0x0f
  );

  const __m128i vinraw = _mm_loadu_si128((__m128i*)input);

  const __m128i vin0 = _mm_shuffle_epi8(vinraw, vshuf0);
  const __m128i vin1 = _mm_shuffle_epi8(vinraw, vshuf1);
  const __m128i vin2 = _mm_shuffle_epi8(vinraw, vshuf2);
  const __m128i vin3 = _mm_shuffle_epi8(vinraw, vshuf3);


  // Calculate row 1, first 4
  __m128i vweight0 = _mm_loadu_si128((__m128i*)&weight[0]);
  __m128i vweight1 = _mm_loadu_si128((__m128i*)&weight[8]);
  __m128i vweight2 = _mm_loadu_si128((__m128i*)&weight[16]);
  __m128i vweight3 = _mm_loadu_si128((__m128i*)&weight[24]);
  __m128i vmadd0 = _mm_madd_epi16(vin0, vweight0);
  __m128i vmadd1 = _mm_madd_epi16(vin1, vweight1);
  __m128i vmadd2 = _mm_madd_epi16(vin2, vweight2);
  __m128i vmadd3 = _mm_madd_epi16(vin3, vweight3);
  __m128i vadd0 = _mm_add_epi32(vmadd0, vmadd1);
  __m128i vadd1 = _mm_add_epi32(vmadd2, vmadd3);
  __m128i result0 = _mm_add_epi32(vadd0, vadd1);
  result0 = _mm_add_epi32(result0, vofs);
  result0 = _mm_srai_epi32(result0, MIP_SHIFT_MATRIX);
  result0 = _mm_add_epi32(result0, vinofs);

  weight += input_size * 4;

  // Calculate row 1, last 4
  vweight0 = _mm_loadu_si128((__m128i*)&weight[0]);
  vweight1 = _mm_loadu_si128((__m128i*)&weight[8]);
  vweight2 = _mm_loadu_si128((__m128i*)&weight[16]);
  vweight3 = _mm_loadu_si128((__m128i*)&weight[24]);
  vmadd0 = _mm_madd_epi16(vin0, vweight0);
  vmadd1 = _mm_madd_epi16(vin1, vweight1);
  vmadd2 = _mm_madd_epi16(vin2, vweight2);
  vmadd3 = _mm_madd_epi16(vin3, vweight3);
  vadd0 = _mm_add_epi32(vmadd0, vmadd1);
  vadd1 = _mm_add_epi32(vmadd2, vmadd3);
  __m128i result1 = _mm_add_epi32(vadd0, vadd1);
  result1 = _mm_add_epi32(result1, vofs);
  result1 = _mm_srai_epi32(result1, MIP_SHIFT_MATRIX);
  result1 = _mm_add_epi32(result1, vinofs);

  __m128i vres16_a = _mm_packus_epi32(result0, result1);


  weight += input_size * 4;
  // Calculate row 2, first 4
  vweight0 = _mm_loadu_si128((__m128i*)&weight[0]);
  vweight1 = _mm_loadu_si128((__m128i*)&weight[8]);
  vweight2 = _mm_loadu_si128((__m128i*)&weight[16]);
  vweight3 = _mm_loadu_si128((__m128i*)&weight[24]);

  vmadd0 = _mm_madd_epi16(vin0, vweight0);
  vmadd1 = _mm_madd_epi16(vin1, vweight1);
  vmadd2 = _mm_madd_epi16(vin2, vweight2);
  vmadd3 = _mm_madd_epi16(vin3, vweight3);

  vadd0 = _mm_add_epi32(vmadd0, vmadd1);
  vadd1 = _mm_add_epi32(vmadd2, vmadd3);

  result0 = _mm_add_epi32(vadd0, vadd1);

  result0 = _mm_add_epi32(result0, vofs);
  result0 = _mm_srai_epi32(result0, MIP_SHIFT_MATRIX);
  result0 = _mm_add_epi32(result0, vinofs);

  weight += input_size * 4;
  // Calculate row 2, last 4
  vweight0 = _mm_loadu_si128((__m128i*)&weight[0]);
  vweight1 = _mm_loadu_si128((__m128i*)&weight[8]);
  vweight2 = _mm_loadu_si128((__m128i*)&weight[16]);
  vweight3 = _mm_loadu_si128((__m128i*)&weight[24]);
  vmadd0 = _mm_madd_epi16(vin0, vweight0);
  vmadd1 = _mm_madd_epi16(vin1, vweight1);
  vmadd2 = _mm_madd_epi16(vin2, vweight2);
  vmadd3 = _mm_madd_epi16(vin3, vweight3);
  vadd0 = _mm_add_epi32(vmadd0, vmadd1);
  vadd1 = _mm_add_epi32(vmadd2, vmadd3);
  result1 = _mm_add_epi32(vadd0, vadd1);
  result1 = _mm_add_epi32(result1, vofs);
  result1 = _mm_srai_epi32(result1, MIP_SHIFT_MATRIX);
  result1 = _mm_add_epi32(result1, vinofs);
  __m128i vres16_b = _mm_packus_epi32(result0, result1);
  __m128i vres8 = _mm_packus_epi16(vres16_a, vres16_b);
  if (transpose) {
    vres8 = _mm_shuffle_epi8(vres8, vtranspose);
    _mm_storeu_si128((__m128i*)output, vres8);

  } else {
    _mm_storeu_si128((__m128i*)output, vres8);
  }
  /*if (transpose) {
    for (int y = 0; y < pred_size; y++) {
      for (int x = 0; x < pred_size; x++) {
        output[y * pred_size + x] = out_ptr[x * pred_size + y];
      }
    }
  }*/
}

// Size ID 2
static void INLINE mip_reduced_pred_sid2_avx2(uvg_pixel* const output,
  const int16_t* const input,
  const uint16_t* matrix,
  const bool transpose,
  const int in_offset,
  const int in_offset_tr)
{
  const int input_size = 8;
  const int pred_size = 8;
  const int size_id = 2;

  uvg_pixel * out_ptr = output;

  int sum = 0;
  for (int i = 0; i < input_size; i++) {
    sum += input[i];
  }
  const int offset = (1 << (MIP_SHIFT_MATRIX - 1)) - MIP_OFFSET_MATRIX * sum;

  const __m128i vofs = _mm_set1_epi32(offset);

  const uint16_t* weight = matrix;
  const int input_offset = transpose ? in_offset_tr : in_offset;

  const __m128i vinofs = _mm_set1_epi32(input_offset);

  const __m128i vshuf0 = _mm_setr_epi8(
    0x00, 0x01, 0x02, 0x03, 0x00, 0x01, 0x02, 0x03,
    0x00, 0x01, 0x02, 0x03, 0x00, 0x01, 0x02, 0x03);
  const __m128i vshuf1 = _mm_setr_epi8(
    0x04, 0x05, 0x06, 0x07, 0x04, 0x05, 0x06, 0x07,
    0x04, 0x05, 0x06, 0x07, 0x04, 0x05, 0x06, 0x07);
  const __m128i vshuf2 = _mm_setr_epi8(
    0x08, 0x09, 0x0a, 0x0b, 0x08, 0x09, 0x0a, 0x0b,
    0x08, 0x09, 0x0a, 0x0b, 0x08, 0x09, 0x0a, 0x0b);
  const __m128i vshuf3 = _mm_setr_epi8(
    0x0c, 0x0d, 0x0e, 0x0f, 0x0c, 0x0d, 0x0e, 0x0f,
    0x0c, 0x0d, 0x0e, 0x0f, 0x0c, 0x0d, 0x0e, 0x0f);

  const __m128i vinraw = _mm_loadu_si128((__m128i*)input);

  const __m128i vin0 = _mm_shuffle_epi8(vinraw, vshuf0);
  const __m128i vin1 = _mm_shuffle_epi8(vinraw, vshuf1);
  const __m128i vin2 = _mm_shuffle_epi8(vinraw, vshuf2);
  const __m128i vin3 = _mm_shuffle_epi8(vinraw, vshuf3);
  __m128i vtranspose = _mm_setr_epi8(
    0x00, 0x08, 0x01, 0x09, 0x02, 0x0a, 0x03, 0x0b,
    0x04, 0x0c, 0x05, 0x0d, 0x06, 0x0e, 0x07, 0x0f
  );
  
   __m128i vtmpres[4];
  
  for (int y = 0, tmp = 0; y < pred_size; y += 2, ++tmp) {
    // Calculate row 1, first 4
    __m128i vweight0 = _mm_loadu_si128((__m128i*) &weight[0]);
    __m128i vweight1 = _mm_loadu_si128((__m128i*) &weight[8]);
    __m128i vweight2 = _mm_loadu_si128((__m128i*) &weight[16]);
    __m128i vweight3 = _mm_loadu_si128((__m128i*) &weight[24]);

    __m128i vmadd0 = _mm_madd_epi16(vin0, vweight0);
    __m128i vmadd1 = _mm_madd_epi16(vin1, vweight1);
    __m128i vmadd2 = _mm_madd_epi16(vin2, vweight2);
    __m128i vmadd3 = _mm_madd_epi16(vin3, vweight3);

    __m128i vadd0 = _mm_add_epi32(vmadd0, vmadd1);
    __m128i vadd1 = _mm_add_epi32(vmadd2, vmadd3);

    __m128i result0 = _mm_add_epi32(vadd0, vadd1);

    result0 = _mm_add_epi32(result0, vofs);
    result0 = _mm_srai_epi32(result0, MIP_SHIFT_MATRIX);
    result0 = _mm_add_epi32(result0, vinofs);

    weight += input_size * 4;

    // Calculate row 1, last 4
    vweight0 = _mm_loadu_si128((__m128i*) &weight[0]);
    vweight1 = _mm_loadu_si128((__m128i*) &weight[8]);
    vweight2 = _mm_loadu_si128((__m128i*) &weight[16]);
    vweight3 = _mm_loadu_si128((__m128i*) &weight[24]);

    vmadd0 = _mm_madd_epi16(vin0, vweight0);
    vmadd1 = _mm_madd_epi16(vin1, vweight1);
    vmadd2 = _mm_madd_epi16(vin2, vweight2);
    vmadd3 = _mm_madd_epi16(vin3, vweight3);

    vadd0 = _mm_add_epi32(vmadd0, vmadd1);
    vadd1 = _mm_add_epi32(vmadd2, vmadd3);

    __m128i result1 = _mm_add_epi32(vadd0, vadd1);

    result1 = _mm_add_epi32(result1, vofs);
    result1 = _mm_srai_epi32(result1, MIP_SHIFT_MATRIX);
    result1 = _mm_add_epi32(result1, vinofs);

    __m128i vres16_a = _mm_packus_epi32(result0, result1);

    weight += input_size * 4;

    // Calculate row 2, first 4
    vweight0 = _mm_loadu_si128((__m128i*) &weight[0]);
    vweight1 = _mm_loadu_si128((__m128i*) &weight[8]);
    vweight2 = _mm_loadu_si128((__m128i*) &weight[16]);
    vweight3 = _mm_loadu_si128((__m128i*) &weight[24]);

    vmadd0 = _mm_madd_epi16(vin0, vweight0);
    vmadd1 = _mm_madd_epi16(vin1, vweight1);
    vmadd2 = _mm_madd_epi16(vin2, vweight2);
    vmadd3 = _mm_madd_epi16(vin3, vweight3);

    vadd0 = _mm_add_epi32(vmadd0, vmadd1);
    vadd1 = _mm_add_epi32(vmadd2, vmadd3);

    result0 = _mm_add_epi32(vadd0, vadd1);

    result0 = _mm_add_epi32(result0, vofs);
    result0 = _mm_srai_epi32(result0, MIP_SHIFT_MATRIX);
    result0 = _mm_add_epi32(result0, vinofs);

    weight += input_size * 4;

    // Calculate row 2, last 4
    vweight0 = _mm_loadu_si128((__m128i*) &weight[0]);
    vweight1 = _mm_loadu_si128((__m128i*) &weight[8]);
    vweight2 = _mm_loadu_si128((__m128i*) &weight[16]);
    vweight3 = _mm_loadu_si128((__m128i*) &weight[24]);

    vmadd0 = _mm_madd_epi16(vin0, vweight0);
    vmadd1 = _mm_madd_epi16(vin1, vweight1);
    vmadd2 = _mm_madd_epi16(vin2, vweight2);
    vmadd3 = _mm_madd_epi16(vin3, vweight3);

    vadd0 = _mm_add_epi32(vmadd0, vmadd1);
    vadd1 = _mm_add_epi32(vmadd2, vmadd3);

    result1 = _mm_add_epi32(vadd0, vadd1);

    result1 = _mm_add_epi32(result1, vofs);
    result1 = _mm_srai_epi32(result1, MIP_SHIFT_MATRIX);
    result1 = _mm_add_epi32(result1, vinofs);

    __m128i vres16_b = _mm_packus_epi32(result0, result1);
    __m128i vres8 = _mm_packus_epi16(vres16_a, vres16_b);

    if (transpose) {
      // Store into temporary storage, transpose later
      vtmpres[tmp] = vres8;
    }
    else {
     _mm_storeu_si128((__m128i*)out_ptr, vres8);
      out_ptr += 16;
    }
  }

  if (transpose) {
    vtmpres[0] = _mm_shuffle_epi8(vtmpres[0], vtranspose);
    vtmpres[1] = _mm_shuffle_epi8(vtmpres[1], vtranspose);
    vtmpres[2] = _mm_shuffle_epi8(vtmpres[2], vtranspose);
    vtmpres[3] = _mm_shuffle_epi8(vtmpres[3], vtranspose);

    __m128i v16lo0 = _mm_unpacklo_epi16(vtmpres[0], vtmpres[1]);
    __m128i v16lo1 = _mm_unpacklo_epi16(vtmpres[2], vtmpres[3]);
    __m128i v16hi0 = _mm_unpackhi_epi16(vtmpres[0], vtmpres[1]);
    __m128i v16hi1 = _mm_unpackhi_epi16(vtmpres[2], vtmpres[3]);

    __m128i v32lo0 = _mm_unpacklo_epi32(v16lo0, v16lo1);
    __m128i v32lo1 = _mm_unpacklo_epi32(v16hi0, v16hi1);
    __m128i v32hi0 = _mm_unpackhi_epi32(v16lo0, v16lo1);
    __m128i v32hi1 = _mm_unpackhi_epi32(v16hi0, v16hi1);
    
    /*__m128i vout0 = _mm_unpacklo_epi64(v32lo0, v32hi0);
    __m128i vout1 = _mm_unpacklo_epi64(v32lo1, v32hi1);
    __m128i vout2 = _mm_unpackhi_epi64(v32lo0, v32hi0);
    __m128i vout3 = _mm_unpackhi_epi64(v32lo1, v32hi1);*/

    _mm_store_si128((__m128i*)(output +  0), v32lo0);
    _mm_store_si128((__m128i*)(output + 16), v32hi0);
    _mm_store_si128((__m128i*)(output + 32), v32lo1);
    _mm_store_si128((__m128i*)(output + 48), v32hi1);
  }
}


// 8x8, size id 1 hor upscale params [4, 4, 1, 4, 1, 16, 2, 2]
static void mip_upsampling_w8_ups2_hor_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref, const uint16_t dst_step, const uint8_t ref_step)
{
  const uint8_t red_pred_size = 4;
  const uint8_t ups_factor = 2; // width / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];
  const int rounding_offset = 1 << (log2_factor - 1);
  
  // Shuffles for result lines 0 and 1
  __m128i vshuf0 = _mm_setr_epi8(
    0xff, 0x00, 0x00, 0x01, 0x01, 0x02, 0x02, 0x03,
    0xff, 0x04, 0x04, 0x05, 0x05, 0x06, 0x06, 0x07
  );

  __m128i vshuf1 = _mm_setr_epi8(
    0x00, 0x00, 0x01, 0x01, 0x02, 0x02, 0x03, 0x03,
    0x04, 0x04, 0x05, 0x05, 0x06, 0x06, 0x07, 0x07
  );

  // Shuffles for result lines 2 and 3
  __m128i vshuf2 = _mm_setr_epi8(
    0xff, 0x08, 0x08, 0x09, 0x09, 0x0a, 0x0a, 0x0b,
    0xff, 0x0c, 0x0c, 0x0d, 0x0d, 0x0e, 0x0e, 0x0f
  );

  __m128i vshuf3 = _mm_setr_epi8(
    0x08, 0x08, 0x09, 0x09, 0x0a, 0x0a, 0x0b, 0x0b,
    0x0c, 0x0c, 0x0d, 0x0d, 0x0e, 0x0e, 0x0f, 0x0f
  );

  __m256i vrnd = _mm256_set1_epi16(rounding_offset);

  uvg_pixel ref0 = *(ref + (ref_step * 1) - 1);
  uvg_pixel ref1 = *(ref + (ref_step * 2) - 1);
  uvg_pixel ref2 = *(ref + (ref_step * 3) - 1);
  uvg_pixel ref3 = *(ref + (ref_step * 4) - 1);

  __m128i vsrc = _mm_loadu_si128((__m128i*)src);

  __m128i vadd0 = _mm_shuffle_epi8(vsrc, vshuf0);
  __m128i vadd1 = _mm_shuffle_epi8(vsrc, vshuf1);
  __m128i vadd2 = _mm_shuffle_epi8(vsrc, vshuf2);
  __m128i vadd3 = _mm_shuffle_epi8(vsrc, vshuf3);

  vadd0 = _mm_insert_epi8(vadd0, ref0, 0x00);
  vadd0 = _mm_insert_epi8(vadd0, ref1, 0x08);
  vadd2 = _mm_insert_epi8(vadd2, ref2, 0x00);
  vadd2 = _mm_insert_epi8(vadd2, ref3, 0x08);

  // Extend to 16-bit
  __m256i vadd16_0 = _mm256_cvtepu8_epi16(vadd0);
  __m256i vadd16_1 = _mm256_cvtepu8_epi16(vadd1);
  __m256i vadd16_2 = _mm256_cvtepu8_epi16(vadd2);
  __m256i vadd16_3 = _mm256_cvtepu8_epi16(vadd3);

  __m256i vtmp0 = _mm256_add_epi16(vadd16_0, vadd16_1);
  __m256i vtmp1 = _mm256_add_epi16(vadd16_2, vadd16_3);

  vtmp0 = _mm256_add_epi16(vtmp0, vrnd);
  vtmp1 = _mm256_add_epi16(vtmp1, vrnd);

  vtmp0 = _mm256_srli_epi16(vtmp0, log2_factor);
  vtmp1 = _mm256_srli_epi16(vtmp1, log2_factor);

  __m256i vres = _mm256_packus_epi16(vtmp0, vtmp1);
  vres = _mm256_permute4x64_epi64(vres, _MM_SHUFFLE(3, 1, 2, 0));

  // Dst step is never 8, since this is only called for 8x8 blocks
  *(uint64_t*)&dst[dst_step * 0] = _mm256_extract_epi64(vres, 0);
  *(uint64_t*)&dst[dst_step * 1] = _mm256_extract_epi64(vres, 1);
  *(uint64_t*)&dst[dst_step * 2] = _mm256_extract_epi64(vres, 2);
  *(uint64_t*)&dst[dst_step * 3] = _mm256_extract_epi64(vres, 3);

  /*if (dst_step == 8) {
    _mm256_storeu_si256((__m256i*)dst, vres);
  }
  else {
    *(uint64_t*)&dst[dst_step * 0] = _mm256_extract_epi64(vres, 0);
    *(uint64_t*)&dst[dst_step * 1] = _mm256_extract_epi64(vres, 1);
    *(uint64_t*)&dst[dst_step * 2] = _mm256_extract_epi64(vres, 2);
    *(uint64_t*)&dst[dst_step * 3] = _mm256_extract_epi64(vres, 3);
  }*/
}

static void mip_upsampling_w8_ups2_hor_avx2_alt(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref, const uint16_t dst_step, const uint8_t ref_step)
{
  const uint8_t red_pred_size = 4;
  const uint8_t ups_factor = 2; // width / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];

  const uvg_pixel* ref_ptr = ref + ref_step - 1;
  const uvg_pixel* src_ptr = src;
  const uvg_pixel* dst_ptr = dst;

  ALIGNED(16) uint8_t before[17];
  memcpy(&before[1], src_ptr, 16);
  before[0] = ref_ptr[ref_step * 0];
  before[4] = ref_ptr[ref_step * 1];
  before[8] = ref_ptr[ref_step * 2];
  before[12] = ref_ptr[ref_step * 3];

  __m128i vbefore = _mm_load_si128((__m128i*)before);
  __m128i vbehind = _mm_load_si128((__m128i*)src_ptr);

  __m128i vavg = _mm_avg_epu8(vbefore, vbehind);

  __m128i vreslo = _mm_unpacklo_epi8(vavg, vbehind);
  __m128i vreshi = _mm_unpackhi_epi8(vavg, vbehind);

  // Dst step is never 8, since this is only called for 8x8 blocks
  *(uint64_t*)&dst[dst_step * 0] = _mm_extract_epi64(vreslo, 0);
  *(uint64_t*)&dst[dst_step * 1] = _mm_extract_epi64(vreslo, 1);
  *(uint64_t*)&dst[dst_step * 2] = _mm_extract_epi64(vreshi, 0);
  *(uint64_t*)&dst[dst_step * 3] = _mm_extract_epi64(vreshi, 1);

  /*if (dst_step == 8) {
    _mm256_storeu_si256((__m256i*)dst, vres);
  }
  else {
    *(uint64_t*)&dst[dst_step * 0] = _mm256_extract_epi64(vres, 0);
    *(uint64_t*)&dst[dst_step * 1] = _mm256_extract_epi64(vres, 1);
    *(uint64_t*)&dst[dst_step * 2] = _mm256_extract_epi64(vres, 2);
    *(uint64_t*)&dst[dst_step * 3] = _mm256_extract_epi64(vres, 3);
  }*/
}

static void mip_upsampling_w16_ups2_hor_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref, const uint16_t dst_step, const uint8_t ref_step)
{
  const uint8_t red_pred_size = 4;
  const uint8_t ups_factor = 2; // width / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];

  const uvg_pixel* ref_ptr = ref + ref_step - 1;
  const uvg_pixel* src_ptr = src;
  const uvg_pixel* dst_ptr = dst;

  for (int i = 0; i < 2; ++i) {
    ALIGNED(32) uint8_t before[33];
    memcpy(&before[1], src_ptr, 32);
    before[0] = ref_ptr[ref_step * 0];
    before[8] = ref_ptr[ref_step * 1];
    before[16] = ref_ptr[ref_step * 2];
    before[24] = ref_ptr[ref_step * 3];

    __m256i vbefore = _mm256_load_si256((__m256i*)before);
    __m256i vbehind = _mm256_loadu_si256((__m256i*)src_ptr);

    __m256i vavg = _mm256_avg_epu8(vbefore, vbehind);

    __m256i vreslo = _mm256_unpacklo_epi8(vavg, vbehind);
    __m256i vreshi = _mm256_unpackhi_epi8(vavg, vbehind);

    _mm_store_si128((__m128i*) & dst_ptr[dst_step * 0], _mm256_extracti128_si256(vreslo, 0));
    _mm_store_si128((__m128i*) & dst_ptr[dst_step * 1], _mm256_extracti128_si256(vreshi, 0));
    _mm_store_si128((__m128i*) & dst_ptr[dst_step * 2], _mm256_extracti128_si256(vreslo, 1));
    _mm_store_si128((__m128i*) & dst_ptr[dst_step * 3], _mm256_extracti128_si256(vreshi, 1));

    src_ptr += 32;
    dst_ptr += dst_step * 4;
    ref_ptr += ref_step * 4;
  }
}

static void mip_upsampling_w16_ups4_hor_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref, const uint16_t dst_step, const uint8_t ref_step)
{
  const uint8_t red_pred_size = 4;
  const uint8_t ups_factor = 4; // width / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];

  const uvg_pixel* ref_ptr = ref + ref_step - 1;
  const uvg_pixel* src_ptr = src;
  const uvg_pixel* dst_ptr = dst;

  int step = ref_step;
  __m128i ones   = _mm_set1_epi8(1);
  __m128i threes = _mm_set1_epi8(3);

  __m256i permute_mask = _mm256_setr_epi32(0, 2, 4, 6, 1, 3, 5, 7);

  // Assign references by hand after copying sources. This will avoid the use of inserts later.
  // Before buffer length is 33 since we need to copy reference value into the first index.
  // Copying 32 samples is faster than copying 31. First indices of each 8 wide row will be replaced
  // with a reference value.
  ALIGNED(16) uint8_t before[17];
  memcpy(&before[1], src_ptr, 16);
  before[0] = ref_ptr[ref_step * 0];
  before[4] = ref_ptr[ref_step * 1];
  before[8] = ref_ptr[ref_step * 2];
  before[12] = ref_ptr[ref_step * 3];

  __m128i vbefore = _mm_load_si128((__m128i*)before);
  __m128i vbehind = _mm_load_si128((__m128i*)src_ptr);

  // Permute the input values to get the result in correct order.
  //vbefore = _mm256_permutevar8x32_epi32(vbefore, permute_mask);
  //vbehind = _mm256_permutevar8x32_epi32(vbehind, permute_mask);

  // Calculate the 3 interpolated values between before and behind, middle, left and right.
  __m128i vmiddle = _mm_avg_epu8(vbefore, vbehind);
  __m128i vleft =   _mm_avg_epu8(vmiddle, vbefore);
  __m128i vright =  _mm_avg_epu8(vmiddle, vbehind);

  // Calculate the two last bits of difference between before and behind. These bits are used to determine if there will be rounding error.
  // Rounding error occurs in the left interpolated value if the two last bits of the difference between before and behind is 0b01.
  __m128i diff = _mm_sub_epi8(vbehind, vbefore);
  diff = _mm_and_si128(diff, threes);
  __m128i mask = _mm_cmpeq_epi8(diff, ones); // The rounding error mask will be generated based on the calculated last bits.
  __m128i sub_amount = _mm_blendv_epi8(_mm_set1_epi8(0), ones, mask);

  vleft = _mm_sub_epi8(vleft, sub_amount);

  // Same rounding error handling for right interpolated values. 
  // Error happens if the two last bits of the difference between before and behind is 0b11.
  mask = _mm_cmpeq_epi8(diff, threes);
  sub_amount = _mm_blendv_epi8(_mm_set1_epi8(0), ones, mask);

  vright = _mm_sub_epi8(vright, sub_amount);

  // Interleave results.
  __m128i left_temp0 = _mm_unpacklo_epi8(vleft, vmiddle);
  __m128i left_temp1 = _mm_unpackhi_epi8(vleft, vmiddle);
  __m128i right_temp0 = _mm_unpacklo_epi8(vright, vbehind);
  __m128i right_temp1 = _mm_unpackhi_epi8(vright, vbehind);

  __m128i vtmp0 = _mm_unpacklo_epi16(left_temp0, right_temp0);
  __m128i vtmp1 = _mm_unpackhi_epi16(left_temp0, right_temp0);
  __m128i vtmp2 = _mm_unpacklo_epi16(left_temp1, right_temp1);
  __m128i vtmp3 = _mm_unpackhi_epi16(left_temp1, right_temp1);

  _mm_store_si128((__m128i*)(dst_ptr + dst_step * 0), vtmp0);
  _mm_store_si128((__m128i*)(dst_ptr + dst_step * 1), vtmp1);
  _mm_store_si128((__m128i*)(dst_ptr + dst_step * 2), vtmp2);
  _mm_store_si128((__m128i*)(dst_ptr + dst_step * 3), vtmp3);
}

static void mip_upsampling_w32_ups4_hor_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref, const uint16_t dst_step, const uint8_t ref_step)
{
  const uint8_t red_pred_size = 8;
  const uint8_t ups_factor = 4; // width / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];
  const int rounding_offset = 1 << (log2_factor - 1);

  __m128i vshuf = _mm_setr_epi8(
    0x00, 0x08, 0x01, 0x09, 0x02, 0x0a, 0x03, 0x0b,
    0x04, 0x0c, 0x05, 0x0d, 0x06, 0x0e, 0x07, 0x0f
  );

  __m128i vrnd = _mm_set1_epi16(rounding_offset);

  ALIGNED(32) int16_t refs[8];
  ALIGNED(32) int16_t srcs[8];
  const uvg_pixel* ref_ptr = ref + ref_step - 1;
  const uvg_pixel* src_ptr = src;

  int step = ref_step;

  for (int i = 0; i < 8; i++) {
    for (int ref = 0; ref < 8; ++ref) {
      refs[ref] = *ref_ptr;
      srcs[ref] = *src_ptr;

      ref_ptr += step;
      src_ptr += red_pred_size;
    }

    __m128i vaccu_ref = _mm_load_si128((__m128i*)refs);
    __m128i vsub_ref = vaccu_ref;
    vaccu_ref = _mm_slli_epi16(vaccu_ref, log2_factor);

    __m128i vaccu_src = _mm_setzero_si128();
    __m128i vadd_src = _mm_load_si128((__m128i*)srcs);

    __m128i vres[4];
    for (int res = 0; res < 4; ++res) {
      vaccu_ref = _mm_sub_epi16(vaccu_ref, vsub_ref);
      vaccu_src = _mm_add_epi16(vaccu_src, vadd_src);
      vres[res] = _mm_add_epi16(vaccu_ref, vaccu_src);
      vres[res] = _mm_add_epi16(vres[res], vrnd);
      vres[res] = _mm_srli_epi16(vres[res], log2_factor);
    }

    __m128i vout0 = _mm_packus_epi16(vres[0], vres[1]);
    __m128i vout1 = _mm_packus_epi16(vres[2], vres[3]);
    vout0 = _mm_shuffle_epi8(vout0, vshuf);
    vout1 = _mm_shuffle_epi8(vout1, vshuf);

    __m128i vtmp16lo = _mm_unpacklo_epi16(vout0, vout1);
    __m128i vtmp16hi = _mm_unpackhi_epi16(vout0, vout1);

    const int dst_offset = i * 4;

    *(uint32_t*)&dst[dst_offset + dst_step * 0] = _mm_extract_epi32(vtmp16lo, 0);
    *(uint32_t*)&dst[dst_offset + dst_step * 1] = _mm_extract_epi32(vtmp16lo, 1);
    *(uint32_t*)&dst[dst_offset + dst_step * 2] = _mm_extract_epi32(vtmp16lo, 2);
    *(uint32_t*)&dst[dst_offset + dst_step * 3] = _mm_extract_epi32(vtmp16lo, 3);
    *(uint32_t*)&dst[dst_offset + dst_step * 4] = _mm_extract_epi32(vtmp16hi, 0);
    *(uint32_t*)&dst[dst_offset + dst_step * 5] = _mm_extract_epi32(vtmp16hi, 1);
    *(uint32_t*)&dst[dst_offset + dst_step * 6] = _mm_extract_epi32(vtmp16hi, 2);
    *(uint32_t*)&dst[dst_offset + dst_step * 7] = _mm_extract_epi32(vtmp16hi, 3);

    ref_ptr = src + i;
    src_ptr = src + i + 1;
    step = red_pred_size; // Switch ref step
  }
}

static void mip_upsampling_w32_ups4_hor_avx2_alt(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref, const uint16_t dst_step, const uint8_t ref_step)
{
  const uint8_t red_pred_size = 8;
  const uint8_t ups_factor = 4; // width / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];

  const uvg_pixel* ref_ptr = ref + ref_step - 1;
  const uvg_pixel* src_ptr = src;
  const uvg_pixel* dst_ptr = dst;

  int step = ref_step;
  __m256i ones = _mm256_set1_epi8(1);
  __m256i threes = _mm256_set1_epi8(3);

  __m256i permute_mask = _mm256_setr_epi32(0, 2, 4, 6, 1, 3, 5, 7);

  // This will process 4 rows at a time. Limit is always 8 rows.
  for (int i = 0; i < 2; ++i) {

    // Assign references by hand after copying sources. This will avoid the use of inserts later.
    // Before buffer length is 33 since we need to copy reference value into the first index.
    // Copying 32 samples is faster than copying 31. First indices of each 8 wide row will be replaced
    // with a reference value.
    ALIGNED(32) uint8_t before[33];
    memcpy(&before[1], src_ptr, 32);
    before[0] =  ref_ptr[ref_step * 0];
    before[8] =  ref_ptr[ref_step * 1];
    before[16] = ref_ptr[ref_step * 2];
    before[24] = ref_ptr[ref_step * 3];
    

    __m256i vbefore = _mm256_load_si256((__m256i*)before);
    __m256i vbehind = _mm256_load_si256((__m256i*)src_ptr);

    // Permute the input values to get the result in correct order.
    vbefore = _mm256_permutevar8x32_epi32(vbefore, permute_mask);
    vbehind = _mm256_permutevar8x32_epi32(vbehind, permute_mask);

    // Calculate the 3 interpolated values between before and behind, middle, left and right.
    __m256i vmiddle = _mm256_avg_epu8(vbefore, vbehind);
    __m256i vleft = _mm256_avg_epu8(vmiddle, vbefore);
    __m256i vright = _mm256_avg_epu8(vmiddle, vbehind);
    
    // Calculate the two last bits of difference between before and behind. These bits are used to determine if there will be rounding error.
    // Rounding error occurs in the left interpolated value if the two last bits of the difference between before and behind is 0b01.
    __m256i diff = _mm256_sub_epi8(vbehind, vbefore);
    diff = _mm256_and_si256(diff, threes);
    __m256i mask = _mm256_cmpeq_epi8(diff, ones); // The rounding error mask will be generated based on the calculated last bits.
    __m256i sub_amount = _mm256_blendv_epi8(_mm256_set1_epi8(0), ones, mask);

    vleft = _mm256_sub_epi8(vleft, sub_amount);

    // Same rounding error handling for right interpolated values. 
    // Error happens if the two last bits of the difference between before and behind is 0b11.
    mask = _mm256_cmpeq_epi8(diff, threes);
    sub_amount = _mm256_blendv_epi8(_mm256_set1_epi8(0), ones, mask);

    vright = _mm256_sub_epi8(vright, sub_amount);

    // Interleave results.
    __m256i left_temp0 = _mm256_unpacklo_epi8(vleft, vmiddle);
    __m256i left_temp1 = _mm256_unpackhi_epi8(vleft, vmiddle);
    __m256i right_temp0 = _mm256_unpacklo_epi8(vright, vbehind);
    __m256i right_temp1 = _mm256_unpackhi_epi8(vright, vbehind);

    __m256i vtmp0 = _mm256_unpacklo_epi16(left_temp0, right_temp0);
    __m256i vtmp1 = _mm256_unpackhi_epi16(left_temp0, right_temp0);
    __m256i vtmp2 = _mm256_unpacklo_epi16(left_temp1, right_temp1);
    __m256i vtmp3 = _mm256_unpackhi_epi16(left_temp1, right_temp1);

    _mm256_storeu_si256((__m256i*)(dst_ptr + dst_step * 0), vtmp0);
    _mm256_storeu_si256((__m256i*)(dst_ptr + dst_step * 1), vtmp1);
    _mm256_storeu_si256((__m256i*)(dst_ptr + dst_step * 2), vtmp2);
    _mm256_storeu_si256((__m256i*)(dst_ptr + dst_step * 3), vtmp3);

    src_ptr += 32;
    ref_ptr += ref_step * 4;
    dst_ptr += dst_step * 4;
  }
}

static void mip_upsampling_w32_ups8_hor_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref, const uint16_t dst_step, const uint8_t ref_step)
{
  const uint8_t red_pred_size = 4;
  const uint8_t ups_factor = 8; // width / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];
  const int rounding_offset = 1 << (log2_factor - 1);

  __m128i vshufsrc = _mm_setr_epi8(
    0xff, 0xff, 0x00, 0x01, 0x02, 0x03, 0x04, 0x05,
    0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d
  );

  __m128i vshuf0 = _mm_setr_epi8(
    0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01,
    0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01
  );

  __m128i vshuf1 = _mm_setr_epi8(
    0x02, 0x03, 0x02, 0x03, 0x02, 0x03, 0x02, 0x03,
    0x02, 0x03, 0x02, 0x03, 0x02, 0x03, 0x02, 0x03
  );

  __m128i vshuf2 = _mm_setr_epi8(
    0x04, 0x05, 0x04, 0x05, 0x04, 0x05, 0x04, 0x05,
    0x04, 0x05, 0x04, 0x05, 0x04, 0x05, 0x04, 0x05
  );

  __m128i vshuf3 = _mm_setr_epi8(
    0x06, 0x07, 0x06, 0x07, 0x06, 0x07, 0x06, 0x07,
    0x06, 0x07, 0x06, 0x07, 0x06, 0x07, 0x06, 0x07
  );

  __m128i vrnd = _mm_set1_epi16(rounding_offset);
  __m128i vmul = _mm_setr_epi16(1, 2, 3, 4, 5, 6, 7, 8);

  const uvg_pixel* ref_ptr = ref + ref_step - 1;
  const uvg_pixel* src_ptr = src;

  int step = ref_step;

  uvg_pixel* dst_ptr = dst;

  for (int i = 0; i < 8; i++) {
    // Handle input data
    int16_t before = *ref_ptr;
    __m128i vtmp = _mm_loadu_si128((__m128i*)src_ptr);
    __m128i vbehind = _mm_cvtepu8_epi16(vtmp);

    __m128i vbefore = vbehind;
    vbefore = _mm_shuffle_epi8(vbefore, vshufsrc);
    vbefore = _mm_insert_epi16(vbefore, before, 0);
    __m128i vbeforeshifted = _mm_slli_epi16(vbefore, log2_factor);

    __m128i vinterpolate = _mm_sub_epi16(vbehind, vbefore);

    // Calculate 1st 16 result chunk
    __m128i vbefore0 = _mm_shuffle_epi8(vbeforeshifted, vshuf0);
    __m128i vinterpolate0 = _mm_shuffle_epi8(vinterpolate, vshuf0);

    __m128i vmulres0 = _mm_mullo_epi16(vinterpolate0, vmul);
    vmulres0 = _mm_add_epi16(vmulres0, vbefore0);

    vmulres0 = _mm_add_epi16(vmulres0, vrnd);
    vmulres0 = _mm_srai_epi16(vmulres0, log2_factor);

    __m128i vbefore1 = _mm_shuffle_epi8(vbeforeshifted, vshuf1);
    __m128i vinterpolate1 = _mm_shuffle_epi8(vinterpolate, vshuf1);

    __m128i vmulres1 = _mm_mullo_epi16(vinterpolate1, vmul);
    vmulres1 = _mm_add_epi16(vmulres1, vbefore1);

    vmulres1 = _mm_add_epi16(vmulres1, vrnd);
    vmulres1 = _mm_srai_epi16(vmulres1, log2_factor);

    __m128i vres = _mm_packus_epi16(vmulres0, vmulres1);

    _mm_store_si128((__m128i*)(dst_ptr + 0), vres);

    // Calculate 2nd 16 result chunk
    vbefore0 = _mm_shuffle_epi8(vbeforeshifted, vshuf2);
    vinterpolate0 = _mm_shuffle_epi8(vinterpolate, vshuf2);

    vmulres0 = _mm_mullo_epi16(vinterpolate0, vmul);
    vmulres0 = _mm_add_epi16(vmulres0, vbefore0);

    vmulres0 = _mm_add_epi16(vmulres0, vrnd);
    vmulres0 = _mm_srai_epi16(vmulres0, log2_factor);

    vbefore1 = _mm_shuffle_epi8(vbeforeshifted, vshuf3);
    vinterpolate1 = _mm_shuffle_epi8(vinterpolate, vshuf3);

    vmulres1 = _mm_mullo_epi16(vinterpolate1, vmul);
    vmulres1 = _mm_add_epi16(vmulres1, vbefore1);

    vmulres1 = _mm_add_epi16(vmulres1, vrnd);
    vmulres1 = _mm_srai_epi16(vmulres1, log2_factor);

    vres = _mm_packus_epi16(vmulres0, vmulres1);

    _mm_store_si128((__m128i*)(dst_ptr + 16), vres);

    dst_ptr += dst_step;
    ref_ptr += ref_step;
    src_ptr += red_pred_size;
  }
}

static void mip_upsampling_w64_ups8_hor_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref, const uint16_t dst_step, const uint8_t ref_step)
{
  const uvg_pixel* ref_ptr = ref + ref_step - 1;
  const uvg_pixel* src_ptr = src;
  const uvg_pixel* dst_ptr = dst;
  
  const __m256i ones = _mm256_set1_epi8(1);
  const __m256i twos = _mm256_set1_epi8(2);
  const __m256i threes = _mm256_set1_epi8(3);
  const __m256i fours = _mm256_set1_epi8(4);
  const __m256i fives = _mm256_set1_epi8(5);
  const __m256i sixes = _mm256_set1_epi8(6);
  const __m256i sevens = _mm256_set1_epi8(7);
  const __m256i eights = _mm256_set1_epi8(8);

  __m256i shuffle_mask = _mm256_setr_epi8(
    0x00, 0x01, 0x04, 0x05, 0x08, 0x09, 0x0c, 0x0d,
    0x02, 0x03, 0x06, 0x07, 0x0a, 0x0b, 0x0e, 0x0f,
    0x00, 0x01, 0x04, 0x05, 0x08, 0x09, 0x0c, 0x0d,
    0x02, 0x03, 0x06, 0x07, 0x0a, 0x0b, 0x0e, 0x0f
  );


  // This will process 2 rows at a time. Limit is always 8 rows.
  for (int i = 0; i < 2; ++i) {

    // Assign references by hand after copying sources. This will avoid the use of inserts later.
    ALIGNED(32) uint8_t before[33];
    memcpy(&before[1], src_ptr, 32);
    before[0] = ref_ptr[ref_step * 0];
    before[8] = ref_ptr[ref_step * 1];
    before[16] = ref_ptr[ref_step * 2];
    before[24] = ref_ptr[ref_step * 3];

    __m256i vbefore = _mm256_load_si256((__m256i*)before);
    __m256i vbehind = _mm256_load_si256((__m256i*)src_ptr);

    // Permute the input values to get the result in correct order.
    vbefore = _mm256_shuffle_epi8(vbefore, shuffle_mask);
    vbehind = _mm256_shuffle_epi8(vbehind, shuffle_mask);
    vbefore = _mm256_permute4x64_epi64(vbefore, _MM_SHUFFLE(3, 1, 2, 0));
    vbehind = _mm256_permute4x64_epi64(vbehind, _MM_SHUFFLE(3, 1, 2, 0));

    // Calculate the 7 interpolated values between before and behind, middle, left and right.
    __m256i vmiddle = _mm256_avg_epu8(vbefore, vbehind);
    __m256i vleft_middle = _mm256_avg_epu8(vmiddle, vbefore);
    __m256i vright_middle = _mm256_avg_epu8(vmiddle, vbehind);
    __m256i vleft_left = _mm256_avg_epu8(vbefore, vleft_middle);
    __m256i vleft_right = _mm256_avg_epu8(vleft_middle, vmiddle);
    __m256i vright_left = _mm256_avg_epu8(vmiddle, vright_middle);
    __m256i vright_right = _mm256_avg_epu8(vright_middle, vbehind);

    // Calculate the three and two last bits of difference between before and behind. These bits are used to determine if there will be rounding error.
    __m256i diff = _mm256_sub_epi8(vbehind, vbefore);
    diff = _mm256_and_si256(diff, sevens);
    __m256i three_diff = _mm256_and_si256(diff, threes);

    // Right side
    __m256i mask = _mm256_cmpgt_epi8(diff, fours);  // The rounding error mask will be generated based on the calculated last bits.
    __m256i sub_amount = _mm256_blendv_epi8(_mm256_set1_epi8(0), ones, mask); // If 5, 6, 7 select one
    vright_right = _mm256_sub_epi8(vright_right, sub_amount);

    mask = _mm256_cmpeq_epi8(three_diff, threes);
    sub_amount = _mm256_blendv_epi8(_mm256_set1_epi8(0), ones, mask); // If 3 or 7 select one
    vright_middle = _mm256_sub_epi8(vright_middle, sub_amount);

    __m256i is_two = _mm256_cmpeq_epi8(diff, twos);
    __m256i is_five = _mm256_cmpeq_epi8(diff, fives);
    mask = _mm256_or_si256(mask, is_two);
    mask = _mm256_or_si256(mask, is_five);
    sub_amount = _mm256_blendv_epi8(_mm256_set1_epi8(0), ones, mask); // If 2, 3, 5, or 7 select one
    vright_left = _mm256_sub_epi8(vright_left, sub_amount);

    // Left side
    diff = _mm256_blendv_epi8(diff, eights, _mm256_cmpeq_epi8(_mm256_set1_epi8(0), diff)); // Replace zeros with eights to enable using GT
    mask = _mm256_cmpgt_epi8(diff, threes);
    sub_amount = _mm256_blendv_epi8(ones, _mm256_set1_epi8(0), mask); // If greater than three select zero
    vleft_left = _mm256_sub_epi8(vleft_left, sub_amount);
    
    mask = _mm256_cmpeq_epi8(three_diff, ones);
    sub_amount = _mm256_blendv_epi8(_mm256_set1_epi8(0), ones, mask); // If 1 or 5 select one
    vleft_middle = _mm256_sub_epi8(vleft_middle, sub_amount);

    __m256i is_three = _mm256_cmpeq_epi8(diff, threes);
    __m256i is_six = _mm256_cmpeq_epi8(diff, sixes);
    mask = _mm256_or_si256(mask, is_three);
    mask = _mm256_or_si256(mask, is_six);
    sub_amount = _mm256_blendv_epi8(_mm256_set1_epi8(0), ones, mask); // If 1, 3, 5, 6 select one
    vleft_right = _mm256_sub_epi8(vleft_right, sub_amount);

    // Interleave results.
    __m256i left_left_temp0 = _mm256_unpacklo_epi8(vleft_left, vleft_middle);
    __m256i left_left_temp1 = _mm256_unpackhi_epi8(vleft_left, vleft_middle);
    __m256i left_right_temp0 = _mm256_unpacklo_epi8(vleft_right, vmiddle);
    __m256i left_right_temp1 = _mm256_unpackhi_epi8(vleft_right, vmiddle);
    __m256i right_left_temp0 = _mm256_unpacklo_epi8(vright_left, vright_middle);
    __m256i right_left_temp1 = _mm256_unpackhi_epi8(vright_left, vright_middle);
    __m256i right_right_temp0 = _mm256_unpacklo_epi8(vright_right, vbehind);
    __m256i right_right_temp1 = _mm256_unpackhi_epi8(vright_right, vbehind);

    __m256i vleft_temp0 = _mm256_unpacklo_epi16(left_left_temp0, left_right_temp0);
    __m256i vleft_temp1 = _mm256_unpackhi_epi16(left_left_temp0, left_right_temp0);
    __m256i vleft_temp2 = _mm256_unpacklo_epi16(left_left_temp1, left_right_temp1);
    __m256i vleft_temp3 = _mm256_unpackhi_epi16(left_left_temp1, left_right_temp1);
    __m256i vright_temp0 = _mm256_unpacklo_epi16(right_left_temp0, right_right_temp0);
    __m256i vright_temp1 = _mm256_unpackhi_epi16(right_left_temp0, right_right_temp0);
    __m256i vright_temp2 = _mm256_unpacklo_epi16(right_left_temp1, right_right_temp1);
    __m256i vright_temp3 = _mm256_unpackhi_epi16(right_left_temp1, right_right_temp1);

    __m256i vtmp0 = _mm256_unpacklo_epi32(vleft_temp0, vright_temp0);
    __m256i vtmp1 = _mm256_unpackhi_epi32(vleft_temp0, vright_temp0);
    __m256i vtmp2 = _mm256_unpacklo_epi32(vleft_temp1, vright_temp1);
    __m256i vtmp3 = _mm256_unpackhi_epi32(vleft_temp1, vright_temp1);
    __m256i vtmp4 = _mm256_unpacklo_epi32(vleft_temp2, vright_temp2);
    __m256i vtmp5 = _mm256_unpackhi_epi32(vleft_temp2, vright_temp2);
    __m256i vtmp6 = _mm256_unpacklo_epi32(vleft_temp3, vright_temp3);
    __m256i vtmp7 = _mm256_unpackhi_epi32(vleft_temp3, vright_temp3);

    _mm256_store_si256((__m256i*)(dst_ptr + dst_step * 0 + 00), vtmp0);
    _mm256_store_si256((__m256i*)(dst_ptr + dst_step * 0 + 32), vtmp1);
    _mm256_store_si256((__m256i*)(dst_ptr + dst_step * 1 + 00), vtmp2);
    _mm256_store_si256((__m256i*)(dst_ptr + dst_step * 1 + 32), vtmp3);
    _mm256_store_si256((__m256i*)(dst_ptr + dst_step * 2 + 00), vtmp4);
    _mm256_store_si256((__m256i*)(dst_ptr + dst_step * 2 + 32), vtmp5);
    _mm256_store_si256((__m256i*)(dst_ptr + dst_step * 3 + 00), vtmp6);
    _mm256_store_si256((__m256i*)(dst_ptr + dst_step * 3 + 32), vtmp7);

    src_ptr += 32;
    ref_ptr += ref_step * 4;
    dst_ptr += dst_step * 4;
  }
}




static void mip_upsampling_w4_ups2_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uint8_t red_pred_size = 4;
  const uint8_t ups_factor = 2; // height / red_pred_size

  int32_t refline = *(int32_t*)ref;
  __m128i vbehind = _mm_loadu_si128((__m128i*)src);
  __m128i vbefore = vbehind;

  vbefore = _mm_slli_si128(vbefore, 4); // Shift left to make room for one 32-bit integer. This could be done with a shuffle, but there should be no performance difference.
  vbefore = _mm_insert_epi32(vbefore, refline, 0);

  __m128i vavg = _mm_avg_epu8(vbefore, vbehind);

  __m128i vres0 = _mm_unpacklo_epi32(vavg, vbehind);
  __m128i vres1 = _mm_unpackhi_epi32(vavg, vbehind);

  _mm_store_si128((__m128i*)(dst +  0), vres0);
  _mm_store_si128((__m128i*)(dst + 16), vres1);  
}

static void mip_upsampling_w4_ups4_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uint8_t red_pred_size = 4;
  const uint8_t ups_factor = 4; // height / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];
  const int rounding_offset = 1 << (log2_factor - 1);

  __m256i vrnd = _mm256_set1_epi16(rounding_offset);

  /*__m128i vshufbefore = _mm_setr_epi8(
    0xff, 0xff, 0xff, 0xff, 0x00, 0x01, 0x02, 0x03,
    0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b
  );*/

  __m256i vshufres = _mm256_setr_epi8(
    0x00, 0x01, 0x02, 0x03, 0x08, 0x09, 0x0a, 0x0b,
    0x04, 0x05, 0x06, 0x07, 0x0c, 0x0d, 0x0e, 0x0f,
    0x00, 0x01, 0x02, 0x03, 0x08, 0x09, 0x0a, 0x0b,
    0x04, 0x05, 0x06, 0x07, 0x0c, 0x0d, 0x0e, 0x0f
  );

  int32_t refline = *(int32_t*)ref;
  //__m128i vidx = _mm_setr_epi32(0, 32, 64, 96);
  //__m128i vbehind = _mm_i32gather_epi32((const int*)src, vidx, 1);
  __m128i vbehind = _mm_loadu_si128((__m128i*)src);
  __m128i vbefore = vbehind;

  vbefore = _mm_slli_si128(vbefore, 4); // Shift left to make room for one 32-bit integer. This could be done with a shuffle, but there should be no performance difference.
  vbefore = _mm_insert_epi32(vbefore, refline, 0);

  __m256i vbefore256 = _mm256_cvtepu8_epi16(vbefore);
  __m256i vbehind256 = _mm256_cvtepu8_epi16(vbehind);

  __m256i vbeforeshifted = _mm256_slli_epi16(vbefore256, log2_factor);

  // Add rounding offset
  vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

  __m256i vinterpolate = _mm256_sub_epi16(vbehind256, vbefore256);

  __m256i vrow0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
  __m256i vrow1 = _mm256_add_epi16(vrow0, vinterpolate);
  __m256i vrow2 = _mm256_add_epi16(vrow1, vinterpolate);

  /*vrow0 = _mm256_add_epi16(vrow0, vrnd);
  vrow1 = _mm256_add_epi16(vrow1, vrnd);
  vrow2 = _mm256_add_epi16(vrow2, vrnd);*/

  vrow0 = _mm256_srai_epi16(vrow0, log2_factor);
  vrow1 = _mm256_srai_epi16(vrow1, log2_factor);
  vrow2 = _mm256_srai_epi16(vrow2, log2_factor);

  __m256i vres0 = _mm256_packus_epi16(vrow0, vrow1);
  __m256i vres1 = _mm256_packus_epi16(vrow2, vbehind256);

  vres0 = _mm256_shuffle_epi8(vres0, vshufres);
  vres1 = _mm256_shuffle_epi8(vres1, vshufres);

  __m256i vlo128 = _mm256_permute2x128_si256(vres0, vres1, 0x20);
  __m256i vhi128 = _mm256_permute2x128_si256(vres0, vres1, 0x31);
  
  vres0 = _mm256_permute4x64_epi64(vlo128, _MM_SHUFFLE(3, 1, 2, 0));
  vres1 = _mm256_permute4x64_epi64(vhi128, _MM_SHUFFLE(3, 1, 2, 0));

  _mm256_store_si256((__m256i*)(dst + 0), vres0);
  _mm256_store_si256((__m256i*)(dst + 32), vres1);
  
}

static void mip_upsampling_w4_ups8_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uint8_t red_pred_size = 4;
  const uint8_t ups_factor = 8; // height / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];
  const int rounding_offset = 1 << (log2_factor - 1);

  __m256i vrnd = _mm256_set1_epi16(rounding_offset);

  /*__m128i vshufbefore = _mm_setr_epi8(
    0xff, 0xff, 0xff, 0xff, 0x00, 0x01, 0x02, 0x03,
    0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b
  );*/

  __m256i vshufres = _mm256_setr_epi8(
    0x00, 0x01, 0x02, 0x03, 0x08, 0x09, 0x0a, 0x0b,
    0x04, 0x05, 0x06, 0x07, 0x0c, 0x0d, 0x0e, 0x0f,
    0x00, 0x01, 0x02, 0x03, 0x08, 0x09, 0x0a, 0x0b,
    0x04, 0x05, 0x06, 0x07, 0x0c, 0x0d, 0x0e, 0x0f
  );

  int32_t refline = *(int32_t*)ref;
  //__m128i vidx = _mm_setr_epi32(0, 32, 64, 96);
  //__m128i vbehind = _mm_i32gather_epi32((const int*)src, vidx, 1);
  __m128i vbehind = _mm_loadu_si128((__m128i*)src);
  __m128i vbefore = vbehind;
  
  vbefore = _mm_slli_si128(vbefore, 4); // Shift left to make room for one 32-bit integer. This could be done with a shuffle, but there should be no performance difference.
  vbefore = _mm_insert_epi32(vbefore, refline, 0);

  __m256i vbefore256 = _mm256_cvtepu8_epi16(vbefore);
  __m256i vbehind256 = _mm256_cvtepu8_epi16(vbehind);

  __m256i vbeforeshifted = _mm256_slli_epi16(vbefore256, log2_factor);

  // Add rounding offset
  vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

  __m256i vinterpolate = _mm256_sub_epi16(vbehind256, vbefore256);

  __m256i vrow0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
  __m256i vrow1 = _mm256_add_epi16(vrow0, vinterpolate);
  __m256i vrow2 = _mm256_add_epi16(vrow1, vinterpolate);
  __m256i vrow3 = _mm256_add_epi16(vrow2, vinterpolate);
  __m256i vrow4 = _mm256_add_epi16(vrow3, vinterpolate);
  __m256i vrow5 = _mm256_add_epi16(vrow4, vinterpolate);
  __m256i vrow6 = _mm256_add_epi16(vrow5, vinterpolate);

  /*
  vrow0 = _mm256_add_epi16(vrow0, vrnd);
  vrow1 = _mm256_add_epi16(vrow1, vrnd);
  vrow2 = _mm256_add_epi16(vrow2, vrnd);
  vrow3 = _mm256_add_epi16(vrow3, vrnd);
  vrow4 = _mm256_add_epi16(vrow4, vrnd);
  vrow5 = _mm256_add_epi16(vrow5, vrnd);
  vrow6 = _mm256_add_epi16(vrow6, vrnd);*/

  vrow0 = _mm256_srai_epi16(vrow0, log2_factor);
  vrow1 = _mm256_srai_epi16(vrow1, log2_factor);
  vrow2 = _mm256_srai_epi16(vrow2, log2_factor);
  vrow3 = _mm256_srai_epi16(vrow3, log2_factor);
  vrow4 = _mm256_srai_epi16(vrow4, log2_factor);
  vrow5 = _mm256_srai_epi16(vrow5, log2_factor);
  vrow6 = _mm256_srai_epi16(vrow6, log2_factor);

  __m256i vres0 = _mm256_packus_epi16(vrow0, vrow1);
  __m256i vres1 = _mm256_packus_epi16(vrow2, vrow3);
  __m256i vres2 = _mm256_packus_epi16(vrow4, vrow5);
  __m256i vres3 = _mm256_packus_epi16(vrow6, vbehind256);

  vres0 = _mm256_shuffle_epi8(vres0, vshufres);
  vres1 = _mm256_shuffle_epi8(vres1, vshufres);
  vres2 = _mm256_shuffle_epi8(vres2, vshufres);
  vres3 = _mm256_shuffle_epi8(vres3, vshufres);

  __m256i vupklo0 = _mm256_unpacklo_epi64(vres0, vres1);
  __m256i vupklo1 = _mm256_unpacklo_epi64(vres2, vres3);
  __m256i vupkhi0 = _mm256_unpackhi_epi64(vres0, vres1);
  __m256i vupkhi1 = _mm256_unpackhi_epi64(vres2, vres3);

  vres0 = _mm256_permute2x128_si256(vupklo0, vupklo1, 0x20);
  vres1 = _mm256_permute2x128_si256(vupkhi0, vupkhi1, 0x20);
  vres2 = _mm256_permute2x128_si256(vupklo0, vupklo1, 0x31);
  vres3 = _mm256_permute2x128_si256(vupkhi0, vupkhi1, 0x31);

  _mm256_store_si256((__m256i*)(dst +  0), vres0);
  _mm256_store_si256((__m256i*)(dst + 32), vres1);
  _mm256_store_si256((__m256i*)(dst + 64), vres2);
  _mm256_store_si256((__m256i*)(dst + 96), vres3);
}

static void mip_upsampling_w8_ups2_h8_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  int64_t refline = *(int64_t*)ref;
  __m128i vidx0 = _mm_set_epi64x(16, 0);
  __m128i vidx1 = _mm_set_epi64x(32, 16);
  __m128i vidx2 = _mm_set_epi64x(48, 32);

  __m128i vbehind0 = _mm_i64gather_epi64((const long long*)src, vidx0, 1);
  __m128i vbefore1 = _mm_i64gather_epi64((const long long*)src, vidx1, 1);
  __m128i vbehind1 = _mm_i64gather_epi64((const long long*)src, vidx2, 1);
  
  __m128i vbefore0 = vbehind0;
  vbefore0 = _mm_slli_si128(vbefore0, 8); // Shift left to make room for one 64-bit integer. This could be done with a shuffle, but there should be no performance difference.
  vbefore0 = _mm_insert_epi64(vbefore0, refline, 0);

  __m128i vavg0 = _mm_avg_epu8(vbefore0, vbehind0);
  __m128i vavg1 = _mm_avg_epu8(vbefore1, vbehind1);

  __m128i vres0 = _mm_unpacklo_epi64(vavg0, vbehind0);
  __m128i vres1 = _mm_unpackhi_epi64(vavg0, vbehind0);
  __m128i vres2 = _mm_unpacklo_epi64(vavg1, vbehind1);
  __m128i vres3 = _mm_unpackhi_epi64(vavg1, vbehind1);
  
  _mm_store_si128((__m128i*)(dst +  0), vres0);
  _mm_store_si128((__m128i*)(dst + 16), vres1);
  _mm_store_si128((__m128i*)(dst + 32), vres2);
  _mm_store_si128((__m128i*)(dst + 48), vres3);
}

static void mip_upsampling_w8_ups2_h16_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  int64_t refline = *(int64_t*)ref;
  __m128i vbehind0 = _mm_load_si128((__m128i*)(src + 0));
  __m128i vbefore1 = _mm_load_si128((__m128i*)(src + 8));
  __m128i vbehind1 = _mm_load_si128((__m128i*)(src + 16));

  __m128i vbefore0 = vbehind0;
  vbefore0 = _mm_slli_si128(vbefore0, 8); // Shift left to make room for one 64-bit integer. This could be done with a shuffle, but there should be no performance difference.
  vbefore0 = _mm_insert_epi64(vbefore0, refline, 0);

  __m128i vavg0 = _mm_avg_epu8(vbefore0, vbehind0);
  __m128i vavg1 = _mm_avg_epu8(vbefore1, vbehind1);

  __m128i vres0 = _mm_unpacklo_epi64(vavg0, vbehind0);
  __m128i vres1 = _mm_unpackhi_epi64(vavg0, vbehind0);
  __m128i vres2 = _mm_unpacklo_epi64(vavg1, vbehind1);
  __m128i vres3 = _mm_unpackhi_epi64(vavg1, vbehind1);

  _mm_store_si128((__m128i*)(dst + 0), vres0);
  _mm_store_si128((__m128i*)(dst + 16), vres1);
  _mm_store_si128((__m128i*)(dst + 32), vres2);
  _mm_store_si128((__m128i*)(dst + 48), vres3);

  vbefore0 = _mm_load_si128((__m128i*)(src + 24));
  vbehind0 = _mm_load_si128((__m128i*)(src + 32));
  vbefore1 = _mm_load_si128((__m128i*)(src + 40));
  vbehind1 = _mm_load_si128((__m128i*)(src + 48));

  vavg0 = _mm_avg_epu8(vbefore0, vbehind0);
  vavg1 = _mm_avg_epu8(vbefore1, vbehind1);

  vres0 = _mm_unpacklo_epi64(vavg0, vbehind0);
  vres1 = _mm_unpackhi_epi64(vavg0, vbehind0);
  vres2 = _mm_unpacklo_epi64(vavg1, vbehind1);
  vres3 = _mm_unpackhi_epi64(vavg1, vbehind1);

  _mm_store_si128((__m128i*)(dst +  64), vres0);
  _mm_store_si128((__m128i*)(dst +  80), vres1);
  _mm_store_si128((__m128i*)(dst +  96), vres2);
  _mm_store_si128((__m128i*)(dst + 112), vres3);
}
//
static void mip_upsampling_w8_ups4_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uint8_t red_pred_size = 8;
  const uint8_t ups_factor = 4; // height / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];
  const int rounding_offset = 1 << (log2_factor - 1);
  __m256i vrnd = _mm256_set1_epi16(rounding_offset);

  int64_t refline = *(int64_t*)ref;
  __m128i vbehind = _mm_loadu_si128((__m128i*)src);
  __m128i vbefore = vbehind;

  vbefore = _mm_slli_si128(vbefore, 8); // Shift left to make room for one 64-bit integer. This could be done with a shuffle, but there should be no performance difference.
  vbefore = _mm_insert_epi64(vbefore, refline, 0);

  __m256i vbefore256 = _mm256_cvtepu8_epi16(vbefore);
  __m256i vbehind256 = _mm256_cvtepu8_epi16(vbehind);

  __m256i vbeforeshifted = _mm256_slli_epi16(vbefore256, log2_factor);

  // Add rounding offset
  vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

  __m256i vinterpolate = _mm256_sub_epi16(vbehind256, vbefore256);

  __m256i vrow0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
  __m256i vrow1 = _mm256_add_epi16(vrow0, vinterpolate);
  __m256i vrow2 = _mm256_add_epi16(vrow1, vinterpolate);

  vrow0 = _mm256_srai_epi16(vrow0, log2_factor);
  vrow1 = _mm256_srai_epi16(vrow1, log2_factor);
  vrow2 = _mm256_srai_epi16(vrow2, log2_factor);

  __m256i vres0 = _mm256_packus_epi16(vrow0, vrow1);
  __m256i vres1 = _mm256_packus_epi16(vrow2, vbehind256);

  __m256i vlo128 = _mm256_permute2x128_si256(vres0, vres1, 0x20);
  __m256i vhi128 = _mm256_permute2x128_si256(vres0, vres1, 0x31);

  _mm256_store_si256((__m256i*)(dst + 0), vlo128);
  _mm256_store_si256((__m256i*)(dst + 32), vhi128);

  for (int i = 1; i < 4; ++i) {
    vbefore = _mm_loadu_si128((__m128i*)(src + (i * 16 - 8)));
    vbehind = _mm_loadu_si128((__m128i*)(src + (i * 16)));
    vbefore256 = _mm256_cvtepu8_epi16(vbefore);
    vbehind256 = _mm256_cvtepu8_epi16(vbehind);

    vbeforeshifted = _mm256_slli_epi16(vbefore256, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    vinterpolate = _mm256_sub_epi16(vbehind256, vbefore256);

    vrow0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    vrow1 = _mm256_add_epi16(vrow0, vinterpolate);
    vrow2 = _mm256_add_epi16(vrow1, vinterpolate);

    vrow0 = _mm256_srai_epi16(vrow0, log2_factor);
    vrow1 = _mm256_srai_epi16(vrow1, log2_factor);
    vrow2 = _mm256_srai_epi16(vrow2, log2_factor);

    vres0 = _mm256_packus_epi16(vrow0, vrow1);
    vres1 = _mm256_packus_epi16(vrow2, vbehind256);

    vlo128 = _mm256_permute2x128_si256(vres0, vres1, 0x20);
    vhi128 = _mm256_permute2x128_si256(vres0, vres1, 0x31);

    _mm256_store_si256((__m256i*)(dst + (i * 64) +  0), vlo128);
    _mm256_store_si256((__m256i*)(dst + (i * 64) + 32), vhi128);
  }
}

static void mip_upsampling_w8_ups8_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uint8_t red_pred_size = 8;
  const uint8_t ups_factor = 8; // height / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];
  const int rounding_offset = 1 << (log2_factor - 1);
  __m256i vrnd = _mm256_set1_epi16(rounding_offset);

  int64_t refline = *(int64_t*)ref;
  __m128i vbehind = _mm_loadu_si128((__m128i*)src);
  __m128i vbefore = vbehind;

  vbefore = _mm_slli_si128(vbefore, 8); // Shift left to make room for one 64-bit integer. This could be done with a shuffle, but there should be no performance difference.
  vbefore = _mm_insert_epi64(vbefore, refline, 0);

  __m256i vbefore256 = _mm256_cvtepu8_epi16(vbefore);
  __m256i vbehind256 = _mm256_cvtepu8_epi16(vbehind);

  __m256i vbeforeshifted = _mm256_slli_epi16(vbefore256, log2_factor);

  // Add rounding offset
  vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

  __m256i vinterpolate = _mm256_sub_epi16(vbehind256, vbefore256);

  __m256i vrow0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
  __m256i vrow1 = _mm256_add_epi16(vrow0, vinterpolate);
  __m256i vrow2 = _mm256_add_epi16(vrow1, vinterpolate);
  __m256i vrow3 = _mm256_add_epi16(vrow2, vinterpolate);
  __m256i vrow4 = _mm256_add_epi16(vrow3, vinterpolate);
  __m256i vrow5 = _mm256_add_epi16(vrow4, vinterpolate);
  __m256i vrow6 = _mm256_add_epi16(vrow5, vinterpolate);

  vrow0 = _mm256_srai_epi16(vrow0, log2_factor);
  vrow1 = _mm256_srai_epi16(vrow1, log2_factor);
  vrow2 = _mm256_srai_epi16(vrow2, log2_factor);
  vrow3 = _mm256_srai_epi16(vrow3, log2_factor);
  vrow4 = _mm256_srai_epi16(vrow4, log2_factor);
  vrow5 = _mm256_srai_epi16(vrow5, log2_factor);
  vrow6 = _mm256_srai_epi16(vrow6, log2_factor);

  __m256i vres0 = _mm256_packus_epi16(vrow0, vrow1);
  __m256i vres1 = _mm256_packus_epi16(vrow2, vrow3);
  __m256i vres2 = _mm256_packus_epi16(vrow4, vrow5);
  __m256i vres3 = _mm256_packus_epi16(vrow6, vbehind256);

  __m256i vlo128a = _mm256_permute2x128_si256(vres0, vres1, 0x20);
  __m256i vlo128b = _mm256_permute2x128_si256(vres2, vres3, 0x20);
  __m256i vhi128a = _mm256_permute2x128_si256(vres0, vres1, 0x31);
  __m256i vhi128b = _mm256_permute2x128_si256(vres2, vres3, 0x31);

  _mm256_store_si256((__m256i*)(dst + 0),  vlo128a);
  _mm256_store_si256((__m256i*)(dst + 32), vlo128b);
  _mm256_store_si256((__m256i*)(dst + 64), vhi128a);
  _mm256_store_si256((__m256i*)(dst + 96), vhi128b);

  for (int i = 1; i < 4; ++i) {
    vbefore = _mm_loadu_si128((__m128i*)(src + (i * 16 - 8)));
    vbehind = _mm_loadu_si128((__m128i*)(src + (i * 16)));
    vbefore256 = _mm256_cvtepu8_epi16(vbefore);
    vbehind256 = _mm256_cvtepu8_epi16(vbehind);

    vbeforeshifted = _mm256_slli_epi16(vbefore256, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    vinterpolate = _mm256_sub_epi16(vbehind256, vbefore256);

    vrow0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    vrow1 = _mm256_add_epi16(vrow0, vinterpolate);
    vrow2 = _mm256_add_epi16(vrow1, vinterpolate);
    vrow3 = _mm256_add_epi16(vrow2, vinterpolate);
    vrow4 = _mm256_add_epi16(vrow3, vinterpolate);
    vrow5 = _mm256_add_epi16(vrow4, vinterpolate);
    vrow6 = _mm256_add_epi16(vrow5, vinterpolate);

    vrow0 = _mm256_srai_epi16(vrow0, log2_factor);
    vrow1 = _mm256_srai_epi16(vrow1, log2_factor);
    vrow2 = _mm256_srai_epi16(vrow2, log2_factor);
    vrow3 = _mm256_srai_epi16(vrow3, log2_factor);
    vrow4 = _mm256_srai_epi16(vrow4, log2_factor);
    vrow5 = _mm256_srai_epi16(vrow5, log2_factor);
    vrow6 = _mm256_srai_epi16(vrow6, log2_factor);

    vres0 = _mm256_packus_epi16(vrow0, vrow1);
    vres1 = _mm256_packus_epi16(vrow2, vrow3);
    vres2 = _mm256_packus_epi16(vrow4, vrow5);
    vres3 = _mm256_packus_epi16(vrow6, vbehind256);

    vlo128a = _mm256_permute2x128_si256(vres0, vres1, 0x20);
    vlo128b = _mm256_permute2x128_si256(vres2, vres3, 0x20);
    vhi128a = _mm256_permute2x128_si256(vres0, vres1, 0x31);
    vhi128b = _mm256_permute2x128_si256(vres2, vres3, 0x31);

    _mm256_store_si256((__m256i*)(dst + (i * 128) +  0), vlo128a);
    _mm256_store_si256((__m256i*)(dst + (i * 128) + 32), vlo128b);
    _mm256_store_si256((__m256i*)(dst + (i * 128) + 64), vhi128a);
    _mm256_store_si256((__m256i*)(dst + (i * 128) + 96), vhi128b);
  }
}

static void mip_upsampling_w16_ups2_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  __m128i vbehind0 = _mm_loadu_si128((__m128i*)(src + 0));
  __m128i vbehind1 = _mm_loadu_si128((__m128i*)(src + 32));
  __m128i vbehind2 = _mm_loadu_si128((__m128i*)(src + 64));
  __m128i vbehind3 = _mm_loadu_si128((__m128i*)(src + 96));
  __m128i vbehind4 = _mm_loadu_si128((__m128i*)(src + 128));
  __m128i vbehind5 = _mm_loadu_si128((__m128i*)(src + 160));
  __m128i vbehind6 = _mm_loadu_si128((__m128i*)(src + 192));
  __m128i vbehind7 = _mm_loadu_si128((__m128i*)(src + 224));

  __m128i vbefore0 = _mm_loadu_si128((__m128i*)ref);
  __m128i vbefore1 = vbehind0;
  __m128i vbefore2 = vbehind1;
  __m128i vbefore3 = vbehind2;
  __m128i vbefore4 = vbehind3;
  __m128i vbefore5 = vbehind4;
  __m128i vbefore6 = vbehind5;
  __m128i vbefore7 = vbehind6;

  __m128i vavg0 = _mm_avg_epu8(vbefore0, vbehind0);
  __m128i vavg1 = _mm_avg_epu8(vbefore1, vbehind1);
  __m128i vavg2 = _mm_avg_epu8(vbefore2, vbehind2);
  __m128i vavg3 = _mm_avg_epu8(vbefore3, vbehind3);
  __m128i vavg4 = _mm_avg_epu8(vbefore4, vbehind4);
  __m128i vavg5 = _mm_avg_epu8(vbefore5, vbehind5);
  __m128i vavg6 = _mm_avg_epu8(vbefore6, vbehind6);
  __m128i vavg7 = _mm_avg_epu8(vbefore7, vbehind7);

  _mm_store_si128((__m128i*)(dst +   0), vavg0);
  _mm_store_si128((__m128i*)(dst +  32), vavg1);
  _mm_store_si128((__m128i*)(dst +  64), vavg2);
  _mm_store_si128((__m128i*)(dst +  96), vavg3);
  _mm_store_si128((__m128i*)(dst + 128), vavg4);
  _mm_store_si128((__m128i*)(dst + 160), vavg5);
  _mm_store_si128((__m128i*)(dst + 192), vavg6);
  _mm_store_si128((__m128i*)(dst + 224), vavg7);
}

static void mip_upsampling_w16_ups4_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uint8_t red_pred_size = 8;
  const uint8_t ups_factor = 4; // height / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];
  const int rounding_offset = 1 << (log2_factor - 1);
  __m256i vrnd = _mm256_set1_epi16(rounding_offset);

  __m256i vbefore256;
  __m256i vbehind256;

  __m128i vbefore = _mm_load_si128((__m128i*)ref);
  vbefore256 = _mm256_cvtepu8_epi16(vbefore);
  
  for (int i = 0; i < 8; ++i) {
    __m128i vbehind = _mm_loadu_si128((__m128i*)(src + (i * 64)));
    vbehind256 = _mm256_cvtepu8_epi16(vbehind);

    __m256i vbeforeshifted = _mm256_slli_epi16(vbefore256, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    __m256i vinterpolate = _mm256_sub_epi16(vbehind256, vbefore256);

    __m256i vrow0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    __m256i vrow1 = _mm256_add_epi16(vrow0, vinterpolate);
    __m256i vrow2 = _mm256_add_epi16(vrow1, vinterpolate);

    vrow0 = _mm256_srai_epi16(vrow0, log2_factor);
    vrow1 = _mm256_srai_epi16(vrow1, log2_factor);
    vrow2 = _mm256_srai_epi16(vrow2, log2_factor);

    __m256i vres0 = _mm256_packus_epi16(vrow0, vrow1);
    __m256i vres1 = _mm256_packus_epi16(vrow2, vbehind256);

    vres0 = _mm256_permute4x64_epi64(vres0, _MM_SHUFFLE(3, 1, 2, 0));
    vres1 = _mm256_permute4x64_epi64(vres1, _MM_SHUFFLE(3, 1, 2, 0));

    _mm256_store_si256((__m256i*)(dst + (i * 64) +  0), vres0);
    _mm256_store_si256((__m256i*)(dst + (i * 64) + 32), vres1);

    vbefore256 = vbehind256;
  }
}

static void mip_upsampling_w16_ups4_ver_avx2_alt(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uvg_pixel* src_ptr = src;
  const uvg_pixel* dst_ptr = dst;

  __m128i vbefore = _mm_load_si128((__m128i*)ref);

  const __m128i zeros  = _mm_setzero_si128();
  const __m128i ones   = _mm_set1_epi8(1);
  const __m128i threes = _mm_set1_epi8(3);

  for (int i = 0; i < 8; ++i) {
    __m128i vbehind = _mm_load_si128((__m128i*)src_ptr);

    // Calculate the 3 interpolated lines between before and behind. Top row, middle row and bottom row.
    __m128i vmiddle = _mm_avg_epu8(vbefore, vbehind);
    __m128i vtop    = _mm_avg_epu8(vbefore, vmiddle);
    __m128i vbottom = _mm_avg_epu8(vmiddle, vbehind);

    // Calculate the two last bits of difference between before and behind. These bits are used to determine if there will be rounding error.
    // Rounding error occurs in the left interpolated value if the two last bits of the difference between before and behind is 0b01.
    __m128i diff       = _mm_sub_epi8(vbehind, vbefore);
            diff       = _mm_and_si128(diff, threes);
    __m128i mask       = _mm_cmpeq_epi8(diff, ones); // The rounding error mask will be generated based on the calculated last bits.
    __m128i sub_amount = _mm_blendv_epi8(zeros, ones, mask);

    vtop = _mm_sub_epi8(vtop, sub_amount);

    // Same rounding error handling for bottom interpolated values. 
    // Error happens if the two last bits of the difference between before and behind is 0b11.
    mask       = _mm_cmpeq_epi8(diff, threes);
    sub_amount = _mm_blendv_epi8(zeros, ones, mask);

    vbottom = _mm_sub_epi8(vbottom, sub_amount);

    // Store results
    _mm_store_si128((__m128i*)(dst_ptr +  0), vtop);
    _mm_store_si128((__m128i*)(dst_ptr + 16), vmiddle);
    _mm_store_si128((__m128i*)(dst_ptr + 32), vbottom);

    vbefore = vbehind;
    src_ptr += 64;
    dst_ptr += 64;
  }
}

static void mip_upsampling_w16_ups8_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uint8_t red_pred_size = 8;
  const uint8_t ups_factor = 8; // height / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];
  const int rounding_offset = 1 << (log2_factor - 1);
  __m256i vrnd = _mm256_set1_epi16(rounding_offset);

  __m256i vbefore256;
  __m256i vbehind256;

  __m128i vbefore = _mm_load_si128((__m128i*)ref);
  vbefore256 = _mm256_cvtepu8_epi16(vbefore);

  for (int i = 0; i < 8; ++i) {
    __m128i vbehind = _mm_loadu_si128((__m128i*)(src + (i * 128)));
    vbehind256 = _mm256_cvtepu8_epi16(vbehind);

    __m256i vbeforeshifted = _mm256_slli_epi16(vbefore256, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    __m256i vinterpolate = _mm256_sub_epi16(vbehind256, vbefore256);

    __m256i vrow0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    __m256i vrow1 = _mm256_add_epi16(vrow0, vinterpolate);
    __m256i vrow2 = _mm256_add_epi16(vrow1, vinterpolate);
    __m256i vrow3 = _mm256_add_epi16(vrow2, vinterpolate);
    __m256i vrow4 = _mm256_add_epi16(vrow3, vinterpolate);
    __m256i vrow5 = _mm256_add_epi16(vrow4, vinterpolate);
    __m256i vrow6 = _mm256_add_epi16(vrow5, vinterpolate);

    vrow0 = _mm256_srai_epi16(vrow0, log2_factor);
    vrow1 = _mm256_srai_epi16(vrow1, log2_factor);
    vrow2 = _mm256_srai_epi16(vrow2, log2_factor);
    vrow3 = _mm256_srai_epi16(vrow3, log2_factor);
    vrow4 = _mm256_srai_epi16(vrow4, log2_factor);
    vrow5 = _mm256_srai_epi16(vrow5, log2_factor);
    vrow6 = _mm256_srai_epi16(vrow6, log2_factor);

    __m256i vres0 = _mm256_packus_epi16(vrow0, vrow1);
    __m256i vres1 = _mm256_packus_epi16(vrow2, vrow3);
    __m256i vres2 = _mm256_packus_epi16(vrow4, vrow5);
    __m256i vres3 = _mm256_packus_epi16(vrow6, vbehind256);

    vres0 = _mm256_permute4x64_epi64(vres0, _MM_SHUFFLE(3, 1, 2, 0));
    vres1 = _mm256_permute4x64_epi64(vres1, _MM_SHUFFLE(3, 1, 2, 0));
    vres2 = _mm256_permute4x64_epi64(vres2, _MM_SHUFFLE(3, 1, 2, 0));
    vres3 = _mm256_permute4x64_epi64(vres3, _MM_SHUFFLE(3, 1, 2, 0));

    _mm256_store_si256((__m256i*)(dst + (i * 128) +  0), vres0);
    _mm256_store_si256((__m256i*)(dst + (i * 128) + 32), vres1);
    _mm256_store_si256((__m256i*)(dst + (i * 128) + 64), vres2);
    _mm256_store_si256((__m256i*)(dst + (i * 128) + 96), vres3);

    vbefore256 = vbehind256;
  }
}

// Note: this alternate version is slower than the original version. It is kept here for reference.
static void mip_upsampling_w16_ups8_ver_avx2_alt(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uvg_pixel* src_ptr = src;
  const uvg_pixel* dst_ptr = dst;

  const __m128i zeros  = _mm_setzero_si128();
  const __m128i ones   = _mm_set1_epi8(1);
  const __m128i twos   = _mm_set1_epi8(2);
  const __m128i threes = _mm_set1_epi8(3);
  const __m128i fours  = _mm_set1_epi8(4);
  const __m128i fives  = _mm_set1_epi8(5);
  const __m128i sixes  = _mm_set1_epi8(6);
  const __m128i sevens = _mm_set1_epi8(7);
  const __m128i eights = _mm_set1_epi8(8);

  __m128i vbefore = _mm_load_si128((__m128i*)ref);

  for (int i = 0; i < 8; ++i) {
    __m128i vbehind = _mm_loadu_si128((__m128i*)src_ptr);

    // Calculate the 7 interpolated lines between before and behind. Ordered by number from top to bottom.
    __m128i vrow3 = _mm_avg_epu8(vbefore, vbehind); // Middle
    __m128i vrow1 = _mm_avg_epu8(vrow3, vbefore);   // Top middle
    __m128i vrow5 = _mm_avg_epu8(vrow3, vbehind);   // Bottom middle
    __m128i vrow0 = _mm_avg_epu8(vbefore, vrow1);   // Top middle top
    __m128i vrow2 = _mm_avg_epu8(vrow1, vrow3);     // Top middle bottom
    __m128i vrow4 = _mm_avg_epu8(vrow3, vrow5);     // Bottom middle top
    __m128i vrow6 = _mm_avg_epu8(vrow5, vbehind);   // Bottom middle bottom

    // Calculate the three and two last bits of difference between before and behind. These bits are used to determine if there will be rounding error.
    __m128i diff       = _mm_sub_epi8(vbehind, vbefore);
            diff       = _mm_and_si128(diff, sevens);
    __m128i three_diff = _mm_and_si128(diff, threes);

    // Bottom side
    __m128i mask       = _mm_cmpgt_epi8(diff, fours);  // The rounding error mask will be generated based on the calculated last bits.
    __m128i sub_amount = _mm_blendv_epi8(zeros, ones, mask); // If 5, 6, 7 select one
            vrow6      = _mm_sub_epi8(vrow6, sub_amount);

    mask       = _mm_cmpeq_epi8(three_diff, threes);
    sub_amount = _mm_blendv_epi8(zeros, ones, mask); // If 3 or 7 select one
    vrow5      = _mm_sub_epi8(vrow5, sub_amount);

    __m128i is_two     = _mm_cmpeq_epi8(diff, twos);
    __m128i is_five    = _mm_cmpeq_epi8(diff, fives);
            mask       = _mm_or_si128(mask, is_two);
            mask       = _mm_or_si128(mask, is_five);
            sub_amount = _mm_blendv_epi8(zeros, ones, mask); // If 2, 3, 5, or 7 select one
            vrow4      = _mm_sub_epi8(vrow4, sub_amount);

    // Top side
    diff       = _mm_blendv_epi8(diff, eights, _mm_cmpeq_epi8(zeros, diff)); // Replace zeros with eights to enable using GT
    mask       = _mm_cmpgt_epi8(diff, threes);
    sub_amount = _mm_blendv_epi8(ones, zeros, mask); // If greater than three select zero
    vrow0      = _mm_sub_epi8(vrow0, sub_amount);

    mask       = _mm_cmpeq_epi8(three_diff, ones);
    sub_amount = _mm_blendv_epi8(zeros, ones, mask); // If 1 or 5 select one
    vrow1      = _mm_sub_epi8(vrow1, sub_amount);

    __m128i is_three   = _mm_cmpeq_epi8(diff, threes);
    __m128i is_six     = _mm_cmpeq_epi8(diff, sixes);
            mask       = _mm_or_si128(mask, is_three);
            mask       = _mm_or_si128(mask, is_six);
            sub_amount = _mm_blendv_epi8(zeros, ones, mask); // If 1, 3, 5, 6 select one
            vrow2      = _mm_sub_epi8(vrow2, sub_amount);
    
    // Store results
    _mm_store_si128((__m128i*)(dst_ptr +  0), vrow0);
    _mm_store_si128((__m128i*)(dst_ptr + 16), vrow1);
    _mm_store_si128((__m128i*)(dst_ptr + 32), vrow2);
    _mm_store_si128((__m128i*)(dst_ptr + 48), vrow3);
    _mm_store_si128((__m128i*)(dst_ptr + 64), vrow4);
    _mm_store_si128((__m128i*)(dst_ptr + 80), vrow5);
    _mm_store_si128((__m128i*)(dst_ptr + 96), vrow6);

    vbefore = vbehind;
    src_ptr += 128;
    dst_ptr += 128;
  }
}

static void mip_upsampling_w32_ups2_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  __m256i vbefore = _mm256_load_si256((__m256i*)ref);

  for (int i = 0; i < 8; ++i) {
    __m256i vbehind = _mm256_load_si256((__m256i*)(src + (i * 64)));
    __m256i vavg = _mm256_avg_epu8(vbefore, vbehind);

    _mm256_store_si256((__m256i*)(dst + (i * 64)), vavg);

    vbefore = vbehind;
  }
}

static void mip_upsampling_w32_ups4_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uint8_t red_pred_size = 8;
  const uint8_t ups_factor = 4; // height / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];
  const int rounding_offset = 1 << (log2_factor - 1);
  __m256i vrnd = _mm256_set1_epi16(rounding_offset);

  __m256i vbefore256a;
  __m256i vbehind256a;

  __m256i vbefore256b;
  __m256i vbehind256b;

  __m128i vbeforea = _mm_load_si128((__m128i*)(ref + 0));
  __m128i vbeforeb = _mm_load_si128((__m128i*)(ref + 16));
  vbefore256a = _mm256_cvtepu8_epi16(vbeforea);
  vbefore256b = _mm256_cvtepu8_epi16(vbeforeb);

  for (int i = 0; i < 8; ++i) {
    __m128i vbehinda = _mm_loadu_si128((__m128i*)(src + (i * 128) + 0));
    __m128i vbehindb = _mm_loadu_si128((__m128i*)(src + (i * 128) + 16));
    vbehind256a = _mm256_cvtepu8_epi16(vbehinda);
    vbehind256b = _mm256_cvtepu8_epi16(vbehindb);

    // Calculate left side of 32 wide lane
    __m256i vbeforeshifted = _mm256_slli_epi16(vbefore256a, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    __m256i vinterpolate = _mm256_sub_epi16(vbehind256a, vbefore256a);

    __m256i vrowleft0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    __m256i vrowleft1 = _mm256_add_epi16(vrowleft0, vinterpolate);
    __m256i vrowleft2 = _mm256_add_epi16(vrowleft1, vinterpolate);

    vrowleft0 = _mm256_srai_epi16(vrowleft0, log2_factor);
    vrowleft1 = _mm256_srai_epi16(vrowleft1, log2_factor);
    vrowleft2 = _mm256_srai_epi16(vrowleft2, log2_factor);


    // Calculate right side of 32 wide lane
    vbeforeshifted = _mm256_slli_epi16(vbefore256b, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    vinterpolate = _mm256_sub_epi16(vbehind256b, vbefore256b);

    __m256i vrowright0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    __m256i vrowright1 = _mm256_add_epi16(vrowright0, vinterpolate);
    __m256i vrowright2 = _mm256_add_epi16(vrowright1, vinterpolate);

    vrowright0 = _mm256_srai_epi16(vrowright0, log2_factor);
    vrowright1 = _mm256_srai_epi16(vrowright1, log2_factor);
    vrowright2 = _mm256_srai_epi16(vrowright2, log2_factor);


    // Store results
    __m256i vres0 = _mm256_packus_epi16(vrowleft0, vrowright0);
    __m256i vres1 = _mm256_packus_epi16(vrowleft1, vrowright1);
    __m256i vres2 = _mm256_packus_epi16(vrowleft2, vrowright2);

    vres0 = _mm256_permute4x64_epi64(vres0, _MM_SHUFFLE(3, 1, 2, 0));
    vres1 = _mm256_permute4x64_epi64(vres1, _MM_SHUFFLE(3, 1, 2, 0));
    vres2 = _mm256_permute4x64_epi64(vres2, _MM_SHUFFLE(3, 1, 2, 0));

    _mm256_store_si256((__m256i*)(dst + (i * 128) +  0), vres0);
    _mm256_store_si256((__m256i*)(dst + (i * 128) + 32), vres1);
    _mm256_store_si256((__m256i*)(dst + (i * 128) + 64), vres2);

    vbefore256a = vbehind256a;
    vbefore256b = vbehind256b;
  }
}

static void mip_upsampling_w32_ups4_ver_avx2_alt(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uvg_pixel* src_ptr = src;
  const uvg_pixel* dst_ptr = dst;

  __m256i vbefore = _mm256_loadu_si256((__m256i*)ref);

  const __m256i zeros = _mm256_setzero_si256();
  const __m256i ones = _mm256_set1_epi8(1);
  const __m256i threes = _mm256_set1_epi8(3);

  for (int i = 0; i < 8; ++i) {
    __m256i vbehind = _mm256_loadu_si256((__m256i*)src_ptr);

    // Calculate the 3 interpolated lines between before and behind. Top row, middle row and bottom row.
    __m256i vmiddle = _mm256_avg_epu8(vbefore, vbehind);
    __m256i vtop    = _mm256_avg_epu8(vbefore, vmiddle);
    __m256i vbottom = _mm256_avg_epu8(vmiddle, vbehind);

    // Calculate the two last bits of difference between before and behind. These bits are used to determine if there will be rounding error.
    // Rounding error occurs in the left interpolated value if the two last bits of the difference between before and behind is 0b01.
    __m256i diff = _mm256_sub_epi8(vbehind, vbefore);
    diff = _mm256_and_si256(diff, threes);
    __m256i mask = _mm256_cmpeq_epi8(diff, ones); // The rounding error mask will be generated based on the calculated last bits.
    __m256i sub_amount = _mm256_blendv_epi8(zeros, ones, mask);

    vtop = _mm256_sub_epi8(vtop, sub_amount);

    // Same rounding error handling for bottom interpolated values. 
    // Error happens if the two last bits of the difference between before and behind is 0b11.
    mask = _mm256_cmpeq_epi8(diff, threes);
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask);

    vbottom = _mm256_sub_epi8(vbottom, sub_amount);

    // Store results
    _mm256_storeu_si256((__m256i*)(dst_ptr +  0), vtop);
    _mm256_storeu_si256((__m256i*)(dst_ptr + 32), vmiddle);
    _mm256_storeu_si256((__m256i*)(dst_ptr + 64), vbottom);

    vbefore = vbehind;
    src_ptr += 128;
    dst_ptr += 128;
  }
}

static void mip_upsampling_w32_ups8_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uint8_t red_pred_size = 8;
  const uint8_t ups_factor = 8; // height / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];
  const int rounding_offset = 1 << (log2_factor - 1);
  __m256i vrnd = _mm256_set1_epi16(rounding_offset);

  __m256i vbefore256a;
  __m256i vbehind256a;

  __m256i vbefore256b;
  __m256i vbehind256b;

  __m128i vbeforea = _mm_load_si128((__m128i*)(ref + 0));
  __m128i vbeforeb = _mm_load_si128((__m128i*)(ref + 16));
  vbefore256a = _mm256_cvtepu8_epi16(vbeforea);
  vbefore256b = _mm256_cvtepu8_epi16(vbeforeb);

  for (int i = 0; i < 8; ++i) {
    __m128i vbehinda = _mm_loadu_si128((__m128i*)(src + (i * 256) + 0));
    __m128i vbehindb = _mm_loadu_si128((__m128i*)(src + (i * 256) + 16));
    vbehind256a = _mm256_cvtepu8_epi16(vbehinda);
    vbehind256b = _mm256_cvtepu8_epi16(vbehindb);

    // Calculate left side of 32 wide lane
    __m256i vbeforeshifted = _mm256_slli_epi16(vbefore256a, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    __m256i vinterpolate = _mm256_sub_epi16(vbehind256a, vbefore256a);

    __m256i vrowleft0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    __m256i vrowleft1 = _mm256_add_epi16(vrowleft0, vinterpolate);
    __m256i vrowleft2 = _mm256_add_epi16(vrowleft1, vinterpolate);
    __m256i vrowleft3 = _mm256_add_epi16(vrowleft2, vinterpolate);
    __m256i vrowleft4 = _mm256_add_epi16(vrowleft3, vinterpolate);
    __m256i vrowleft5 = _mm256_add_epi16(vrowleft4, vinterpolate);
    __m256i vrowleft6 = _mm256_add_epi16(vrowleft5, vinterpolate);

    vrowleft0 = _mm256_srai_epi16(vrowleft0, log2_factor);
    vrowleft1 = _mm256_srai_epi16(vrowleft1, log2_factor);
    vrowleft2 = _mm256_srai_epi16(vrowleft2, log2_factor);
    vrowleft3 = _mm256_srai_epi16(vrowleft3, log2_factor);
    vrowleft4 = _mm256_srai_epi16(vrowleft4, log2_factor);
    vrowleft5 = _mm256_srai_epi16(vrowleft5, log2_factor);
    vrowleft6 = _mm256_srai_epi16(vrowleft6, log2_factor);


    // Calculate right side of 32 wide lane
    vbeforeshifted = _mm256_slli_epi16(vbefore256b, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    vinterpolate = _mm256_sub_epi16(vbehind256b, vbefore256b);

    __m256i vrowright0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    __m256i vrowright1 = _mm256_add_epi16(vrowright0, vinterpolate);
    __m256i vrowright2 = _mm256_add_epi16(vrowright1, vinterpolate);
    __m256i vrowright3 = _mm256_add_epi16(vrowright2, vinterpolate);
    __m256i vrowright4 = _mm256_add_epi16(vrowright3, vinterpolate);
    __m256i vrowright5 = _mm256_add_epi16(vrowright4, vinterpolate);
    __m256i vrowright6 = _mm256_add_epi16(vrowright5, vinterpolate);

    vrowright0 = _mm256_srai_epi16(vrowright0, log2_factor);
    vrowright1 = _mm256_srai_epi16(vrowright1, log2_factor);
    vrowright2 = _mm256_srai_epi16(vrowright2, log2_factor);
    vrowright3 = _mm256_srai_epi16(vrowright3, log2_factor);
    vrowright4 = _mm256_srai_epi16(vrowright4, log2_factor);
    vrowright5 = _mm256_srai_epi16(vrowright5, log2_factor);
    vrowright6 = _mm256_srai_epi16(vrowright6, log2_factor);

    
    // Store results
    __m256i vres0 = _mm256_packus_epi16(vrowleft0, vrowright0);
    __m256i vres1 = _mm256_packus_epi16(vrowleft1, vrowright1);
    __m256i vres2 = _mm256_packus_epi16(vrowleft2, vrowright2);
    __m256i vres3 = _mm256_packus_epi16(vrowleft3, vrowright3);
    __m256i vres4 = _mm256_packus_epi16(vrowleft4, vrowright4);
    __m256i vres5 = _mm256_packus_epi16(vrowleft5, vrowright5);
    __m256i vres6 = _mm256_packus_epi16(vrowleft6, vrowright6);

    vres0 = _mm256_permute4x64_epi64(vres0, _MM_SHUFFLE(3, 1, 2, 0));
    vres1 = _mm256_permute4x64_epi64(vres1, _MM_SHUFFLE(3, 1, 2, 0));
    vres2 = _mm256_permute4x64_epi64(vres2, _MM_SHUFFLE(3, 1, 2, 0));
    vres3 = _mm256_permute4x64_epi64(vres3, _MM_SHUFFLE(3, 1, 2, 0));
    vres4 = _mm256_permute4x64_epi64(vres4, _MM_SHUFFLE(3, 1, 2, 0));
    vres5 = _mm256_permute4x64_epi64(vres5, _MM_SHUFFLE(3, 1, 2, 0));
    vres6 = _mm256_permute4x64_epi64(vres6, _MM_SHUFFLE(3, 1, 2, 0));

    _mm256_store_si256((__m256i*)(dst + (i * 256) +   0), vres0);
    _mm256_store_si256((__m256i*)(dst + (i * 256) +  32), vres1);
    _mm256_store_si256((__m256i*)(dst + (i * 256) +  64), vres2);
    _mm256_store_si256((__m256i*)(dst + (i * 256) +  96), vres3);
    _mm256_store_si256((__m256i*)(dst + (i * 256) + 128), vres4);
    _mm256_store_si256((__m256i*)(dst + (i * 256) + 160), vres5);
    _mm256_store_si256((__m256i*)(dst + (i * 256) + 192), vres6);

    vbefore256a = vbehind256a;
    vbefore256b = vbehind256b;
  }
}

static void mip_upsampling_w32_ups8_ver_avx2_alt(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uvg_pixel* src_ptr = src;
  const uvg_pixel* dst_ptr = dst;

  const __m256i zeros = _mm256_setzero_si256();
  const __m256i ones = _mm256_set1_epi8(1);
  const __m256i twos = _mm256_set1_epi8(2);
  const __m256i threes = _mm256_set1_epi8(3);
  const __m256i fours = _mm256_set1_epi8(4);
  const __m256i fives = _mm256_set1_epi8(5);
  const __m256i sixes = _mm256_set1_epi8(6);
  const __m256i sevens = _mm256_set1_epi8(7);
  const __m256i eights = _mm256_set1_epi8(8);

  __m256i vbefore = _mm256_load_si256((__m256i*)(ref + 0));

  for (int i = 0; i < 8; ++i) {
    __m256i vbehind = _mm256_load_si256((__m256i*)src_ptr);

    // Calculate the 7 interpolated lines between before and behind. Ordered by number from top to bottom.
    __m256i vrow3 = _mm256_avg_epu8(vbefore, vbehind); // Middle
    __m256i vrow1 = _mm256_avg_epu8(vrow3, vbefore);   // Top middle
    __m256i vrow5 = _mm256_avg_epu8(vrow3, vbehind);   // Bottom middle
    __m256i vrow0 = _mm256_avg_epu8(vbefore, vrow1);   // Top middle top
    __m256i vrow2 = _mm256_avg_epu8(vrow1, vrow3);     // Top middle bottom
    __m256i vrow4 = _mm256_avg_epu8(vrow3, vrow5);     // Bottom middle top
    __m256i vrow6 = _mm256_avg_epu8(vrow5, vbehind);   // Bottom middle bottom

    // Calculate the three and two last bits of difference between before and behind. These bits are used to determine if there will be rounding error.
    __m256i diff = _mm256_sub_epi8(vbehind, vbefore);
    diff = _mm256_and_si256(diff, sevens);
    __m256i three_diff = _mm256_and_si256(diff, threes);

    // Bottom side
    __m256i mask = _mm256_cmpgt_epi8(diff, fours);  // The rounding error mask will be generated based on the calculated last bits.
    __m256i sub_amount = _mm256_blendv_epi8(zeros, ones, mask); // If 5, 6, 7 select one
    vrow6 = _mm256_sub_epi8(vrow6, sub_amount);

    mask = _mm256_cmpeq_epi8(three_diff, threes);
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask); // If 3 or 7 select one
    vrow5 = _mm256_sub_epi8(vrow5, sub_amount);

    __m256i is_two = _mm256_cmpeq_epi8(diff, twos);
    __m256i is_five = _mm256_cmpeq_epi8(diff, fives);
    mask = _mm256_or_si256(mask, is_two);
    mask = _mm256_or_si256(mask, is_five);
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask); // If 2, 3, 5, or 7 select one
    vrow4 = _mm256_sub_epi8(vrow4, sub_amount);

    // Top side
    diff = _mm256_blendv_epi8(diff, eights, _mm256_cmpeq_epi8(zeros, diff)); // Replace zeros with eights to enable using GT
    mask = _mm256_cmpgt_epi8(diff, threes);
    sub_amount = _mm256_blendv_epi8(ones, zeros, mask); // If greater than three select zero
    vrow0 = _mm256_sub_epi8(vrow0, sub_amount);

    mask = _mm256_cmpeq_epi8(three_diff, ones);
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask); // If 1 or 5 select one
    vrow1 = _mm256_sub_epi8(vrow1, sub_amount);

    __m256i is_three = _mm256_cmpeq_epi8(diff, threes);
    __m256i is_six = _mm256_cmpeq_epi8(diff, sixes);
    mask = _mm256_or_si256(mask, is_three);
    mask = _mm256_or_si256(mask, is_six);
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask); // If 1, 3, 5, 6 select one
    vrow2 = _mm256_sub_epi8(vrow2, sub_amount);

    // Store results
    _mm256_store_si256((__m256i*)(dst_ptr +   0), vrow0);
    _mm256_store_si256((__m256i*)(dst_ptr +  32), vrow1);
    _mm256_store_si256((__m256i*)(dst_ptr +  64), vrow2);
    _mm256_store_si256((__m256i*)(dst_ptr +  96), vrow3);
    _mm256_store_si256((__m256i*)(dst_ptr + 128), vrow4);
    _mm256_store_si256((__m256i*)(dst_ptr + 160), vrow5);
    _mm256_store_si256((__m256i*)(dst_ptr + 192), vrow6);

    vbefore = vbehind;
    src_ptr += 256;
    dst_ptr += 256;
  }
}

static void mip_upsampling_w64_ups2_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  __m256i vbeforeleft  = _mm256_load_si256((__m256i*)(ref + 0));
  __m256i vbeforeright = _mm256_load_si256((__m256i*)(ref + 32));

  for (int i = 0; i < 8; ++i) {
    __m256i vbehindleft  = _mm256_load_si256((__m256i*)(src + (i * 128) +  0));
    __m256i vbehindright = _mm256_load_si256((__m256i*)(src + (i * 128) + 32));
    __m256i vavgleft = _mm256_avg_epu8(vbeforeleft, vbehindleft);
    __m256i vavgright = _mm256_avg_epu8(vbeforeright, vbehindright);

    _mm256_store_si256((__m256i*)(dst + (i * 128) +  0), vavgleft);
    _mm256_store_si256((__m256i*)(dst + (i * 128) + 32), vavgright);

    vbeforeleft = vbehindleft;
    vbeforeright = vbehindright;
  }
}

static void mip_upsampling_w64_ups4_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uint8_t red_pred_size = 8;
  const uint8_t ups_factor = 4; // height / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];
  const int rounding_offset = 1 << (log2_factor - 1);
  __m256i vrnd = _mm256_set1_epi16(rounding_offset);

  __m256i vbefore256a;
  __m256i vbehind256a;

  __m256i vbefore256b;
  __m256i vbehind256b;

  __m256i vbefore256c;
  __m256i vbehind256c;

  __m256i vbefore256d;
  __m256i vbehind256d;

  __m128i vbeforea = _mm_load_si128((__m128i*)(ref + 0));
  __m128i vbeforeb = _mm_load_si128((__m128i*)(ref + 16));
  __m128i vbeforec = _mm_load_si128((__m128i*)(ref + 32));
  __m128i vbefored = _mm_load_si128((__m128i*)(ref + 48));
  vbefore256a = _mm256_cvtepu8_epi16(vbeforea);
  vbefore256b = _mm256_cvtepu8_epi16(vbeforeb);
  vbefore256c = _mm256_cvtepu8_epi16(vbeforec);
  vbefore256d = _mm256_cvtepu8_epi16(vbefored);

  for (int i = 0; i < 8; ++i) {
    __m128i vbehinda = _mm_loadu_si128((__m128i*)(src + (i * 256) + 0));
    __m128i vbehindb = _mm_loadu_si128((__m128i*)(src + (i * 256) + 16));
    __m128i vbehindc = _mm_loadu_si128((__m128i*)(src + (i * 256) + 32));
    __m128i vbehindd = _mm_loadu_si128((__m128i*)(src + (i * 256) + 48));
    vbehind256a = _mm256_cvtepu8_epi16(vbehinda);
    vbehind256b = _mm256_cvtepu8_epi16(vbehindb);
    vbehind256c = _mm256_cvtepu8_epi16(vbehindc);
    vbehind256d = _mm256_cvtepu8_epi16(vbehindd);

    // Calculate 1/4 part of 64 wide lane
    __m256i vbeforeshifted = _mm256_slli_epi16(vbefore256a, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    __m256i vinterpolate = _mm256_sub_epi16(vbehind256a, vbefore256a);

    __m256i vrowa0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    __m256i vrowa1 = _mm256_add_epi16(vrowa0, vinterpolate);
    __m256i vrowa2 = _mm256_add_epi16(vrowa1, vinterpolate);

    vrowa0 = _mm256_srai_epi16(vrowa0, log2_factor);
    vrowa1 = _mm256_srai_epi16(vrowa1, log2_factor);
    vrowa2 = _mm256_srai_epi16(vrowa2, log2_factor);


    // Calculate 2/4 part of 64 wide lane
    vbeforeshifted = _mm256_slli_epi16(vbefore256b, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    vinterpolate = _mm256_sub_epi16(vbehind256b, vbefore256b);

    __m256i vrowb0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    __m256i vrowb1 = _mm256_add_epi16(vrowb0, vinterpolate);
    __m256i vrowb2 = _mm256_add_epi16(vrowb1, vinterpolate);

    vrowb0 = _mm256_srai_epi16(vrowb0, log2_factor);
    vrowb1 = _mm256_srai_epi16(vrowb1, log2_factor);
    vrowb2 = _mm256_srai_epi16(vrowb2, log2_factor);


    // Calculate 3/4 part of 64 wide lane
    vbeforeshifted = _mm256_slli_epi16(vbefore256c, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    vinterpolate = _mm256_sub_epi16(vbehind256c, vbefore256c);

    __m256i vrowc0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    __m256i vrowc1 = _mm256_add_epi16(vrowc0, vinterpolate);
    __m256i vrowc2 = _mm256_add_epi16(vrowc1, vinterpolate);

    vrowc0 = _mm256_srai_epi16(vrowc0, log2_factor);
    vrowc1 = _mm256_srai_epi16(vrowc1, log2_factor);
    vrowc2 = _mm256_srai_epi16(vrowc2, log2_factor);


    // Calculate 3/4 part of 64 wide lane
    vbeforeshifted = _mm256_slli_epi16(vbefore256d, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    vinterpolate = _mm256_sub_epi16(vbehind256d, vbefore256d);

    __m256i vrowd0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    __m256i vrowd1 = _mm256_add_epi16(vrowd0, vinterpolate);
    __m256i vrowd2 = _mm256_add_epi16(vrowd1, vinterpolate);

    vrowd0 = _mm256_srai_epi16(vrowd0, log2_factor);
    vrowd1 = _mm256_srai_epi16(vrowd1, log2_factor);
    vrowd2 = _mm256_srai_epi16(vrowd2, log2_factor);


    // Store results
    __m256i vres0left  = _mm256_packus_epi16(vrowa0, vrowb0);
    __m256i vres0right = _mm256_packus_epi16(vrowc0, vrowd0);
    __m256i vres1left  = _mm256_packus_epi16(vrowa1, vrowb1);
    __m256i vres1right = _mm256_packus_epi16(vrowc1, vrowd1);
    __m256i vres2left  = _mm256_packus_epi16(vrowa2, vrowb2);
    __m256i vres2right = _mm256_packus_epi16(vrowc2, vrowd2);

    /*vres0 = _mm256_permute4x64_epi64(vres0, _MM_SHUFFLE(3, 1, 2, 0));
    vres1 = _mm256_permute4x64_epi64(vres1, _MM_SHUFFLE(3, 1, 2, 0));
    vres2 = _mm256_permute4x64_epi64(vres2, _MM_SHUFFLE(3, 1, 2, 0));*/

    vres0left  = _mm256_permute4x64_epi64(vres0left,  _MM_SHUFFLE(3, 1, 2, 0));
    vres0right = _mm256_permute4x64_epi64(vres0right, _MM_SHUFFLE(3, 1, 2, 0));
    vres1left  = _mm256_permute4x64_epi64(vres1left,  _MM_SHUFFLE(3, 1, 2, 0));
    vres1right = _mm256_permute4x64_epi64(vres1right, _MM_SHUFFLE(3, 1, 2, 0));
    vres2left  = _mm256_permute4x64_epi64(vres2left,  _MM_SHUFFLE(3, 1, 2, 0));
    vres2right = _mm256_permute4x64_epi64(vres2right, _MM_SHUFFLE(3, 1, 2, 0));

    /*_mm256_store_si256((__m256i*)(dst + (i * 128) +  0), vres0);
    _mm256_store_si256((__m256i*)(dst + (i * 128) + 32), vres1);
    _mm256_store_si256((__m256i*)(dst + (i * 128) + 64), vres2);*/

    _mm256_store_si256((__m256i*)(dst + (i * 256) +   0), vres0left);
    _mm256_store_si256((__m256i*)(dst + (i * 256) +  32), vres0right);
    _mm256_store_si256((__m256i*)(dst + (i * 256) +  64), vres1left);
    _mm256_store_si256((__m256i*)(dst + (i * 256) +  96), vres1right);
    _mm256_store_si256((__m256i*)(dst + (i * 256) + 128), vres2left);
    _mm256_store_si256((__m256i*)(dst + (i * 256) + 160), vres2right);

    vbefore256a = vbehind256a;
    vbefore256b = vbehind256b;
    vbefore256c = vbehind256c;
    vbefore256d = vbehind256d;
  }
}

static void mip_upsampling_w64_ups4_ver_avx2_alt(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uvg_pixel* src_ptr = src;
  const uvg_pixel* dst_ptr = dst;

  __m256i vbeforeleft = _mm256_load_si256((__m256i*)(ref + 0));
  __m256i vbeforeright = _mm256_load_si256((__m256i*)(ref + 32));

  const __m256i zeros = _mm256_setzero_si256();
  const __m256i ones = _mm256_set1_epi8(1);
  const __m256i threes = _mm256_set1_epi8(3);

  for (int i = 0; i < 8; ++i) {
    // Calculate 4 lines at a time
    __m256i vbehindleft  = _mm256_load_si256((__m256i*)(src_ptr + 0));
    __m256i vbehindright = _mm256_load_si256((__m256i*)(src_ptr + 32));

    // Calculate left side of 64 wide lane
    // Calculate the 3 interpolated lines between before and behind. Top row, middle row and bottom row.
    __m256i vmiddleleft = _mm256_avg_epu8(vbeforeleft, vbehindleft);
    __m256i vtopleft = _mm256_avg_epu8(vbeforeleft, vmiddleleft);
    __m256i vbottomleft = _mm256_avg_epu8(vmiddleleft, vbehindleft);

    // Calculate the two last bits of difference between before and behind. These bits are used to determine if there will be rounding error.
    // Rounding error occurs in the left interpolated value if the two last bits of the difference between before and behind is 0b01.
    __m256i diff = _mm256_sub_epi8(vbehindleft, vbeforeleft);
    diff = _mm256_and_si256(diff, threes);
    __m256i mask = _mm256_cmpeq_epi8(diff, ones); // The rounding error mask will be generated based on the calculated last bits.
    __m256i sub_amount = _mm256_blendv_epi8(zeros, ones, mask);

    vtopleft = _mm256_sub_epi8(vtopleft, sub_amount);

    // Same rounding error handling for bottom interpolated values. 
    // Error happens if the two last bits of the difference between before and behind is 0b11.
    mask = _mm256_cmpeq_epi8(diff, threes);
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask);

    vbottomleft = _mm256_sub_epi8(vbottomleft, sub_amount);


    // Calculate right side of 64 wide lane
    // Calculate the 3 interpolated lines between before and behind. Top row, middle row and bottom row.
    __m256i vmiddleright = _mm256_avg_epu8(vbeforeright, vbehindright);
    __m256i vtopright = _mm256_avg_epu8(vbeforeright, vmiddleright);
    __m256i vbottomright = _mm256_avg_epu8(vmiddleright, vbehindright);

    // Calculate the two last bits of difference between before and behind. These bits are used to determine if there will be rounding error.
    // Rounding error occurs in the right interpolated value if the two last bits of the difference between before and behind is 0b01.
    diff = _mm256_sub_epi8(vbehindright, vbeforeright);
    diff = _mm256_and_si256(diff, threes);
    mask = _mm256_cmpeq_epi8(diff, ones); // The rounding error mask will be generated based on the calculated last bits.
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask);

    vtopright = _mm256_sub_epi8(vtopright, sub_amount);

    // Same rounding error handling for bottom interpolated values. 
    // Error happens if the two last bits of the difference between before and behind is 0b11.
    mask = _mm256_cmpeq_epi8(diff, threes);
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask);

    vbottomright = _mm256_sub_epi8(vbottomright, sub_amount);

    // Store results
    _mm256_store_si256((__m256i*)(dst_ptr + 0), vtopleft);
    _mm256_store_si256((__m256i*)(dst_ptr + 32), vtopright);
    _mm256_store_si256((__m256i*)(dst_ptr + 64), vmiddleleft);
    _mm256_store_si256((__m256i*)(dst_ptr + 96), vmiddleright);
    _mm256_store_si256((__m256i*)(dst_ptr + 128), vbottomleft);
    _mm256_store_si256((__m256i*)(dst_ptr + 160), vbottomright);
    // No need to store the last line of the 4 lines as it is already present in the result array and it was not modified in any way.

    vbeforeleft = vbehindleft;
    vbeforeright = vbehindright;

    dst_ptr += 256;
    src_ptr += 256;
  }
}

static void mip_upsampling_w64_ups8_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uint8_t red_pred_size = 8;
  const uint8_t ups_factor = 8; // height / red_pred_size

  const int log2_factor = uvg_g_convert_to_log2[ups_factor];
  const int rounding_offset = 1 << (log2_factor - 1);
  __m256i vrnd = _mm256_set1_epi16(rounding_offset);

  __m256i vbefore256a;
  __m256i vbehind256a;

  __m256i vbefore256b;
  __m256i vbehind256b;

  __m256i vbefore256c;
  __m256i vbehind256c;

  __m256i vbefore256d;
  __m256i vbehind256d;

  __m128i vbeforea = _mm_load_si128((__m128i*)(ref + 0));
  __m128i vbeforeb = _mm_load_si128((__m128i*)(ref + 16));
  __m128i vbeforec = _mm_load_si128((__m128i*)(ref + 32));
  __m128i vbefored = _mm_load_si128((__m128i*)(ref + 48));
  vbefore256a = _mm256_cvtepu8_epi16(vbeforea);
  vbefore256b = _mm256_cvtepu8_epi16(vbeforeb);
  vbefore256c = _mm256_cvtepu8_epi16(vbeforec);
  vbefore256d = _mm256_cvtepu8_epi16(vbefored);

  for (int i = 0; i < 8; ++i) {
    __m128i vbehinda = _mm_loadu_si128((__m128i*)(src + (i * 512) + 0));
    __m128i vbehindb = _mm_loadu_si128((__m128i*)(src + (i * 512) + 16));
    __m128i vbehindc = _mm_loadu_si128((__m128i*)(src + (i * 512) + 32));
    __m128i vbehindd = _mm_loadu_si128((__m128i*)(src + (i * 512) + 48));
    vbehind256a = _mm256_cvtepu8_epi16(vbehinda);
    vbehind256b = _mm256_cvtepu8_epi16(vbehindb);
    vbehind256c = _mm256_cvtepu8_epi16(vbehindc);
    vbehind256d = _mm256_cvtepu8_epi16(vbehindd);

    // Calculate 1/4 part of 64 wide lane
    __m256i vbeforeshifted = _mm256_slli_epi16(vbefore256a, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    __m256i vinterpolate = _mm256_sub_epi16(vbehind256a, vbefore256a);

    __m256i vrowa0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    __m256i vrowa1 = _mm256_add_epi16(vrowa0, vinterpolate);
    __m256i vrowa2 = _mm256_add_epi16(vrowa1, vinterpolate);
    __m256i vrowa3 = _mm256_add_epi16(vrowa2, vinterpolate);
    __m256i vrowa4 = _mm256_add_epi16(vrowa3, vinterpolate);
    __m256i vrowa5 = _mm256_add_epi16(vrowa4, vinterpolate);
    __m256i vrowa6 = _mm256_add_epi16(vrowa5, vinterpolate);

    vrowa0 = _mm256_srai_epi16(vrowa0, log2_factor);
    vrowa1 = _mm256_srai_epi16(vrowa1, log2_factor);
    vrowa2 = _mm256_srai_epi16(vrowa2, log2_factor);
    vrowa3 = _mm256_srai_epi16(vrowa3, log2_factor);
    vrowa4 = _mm256_srai_epi16(vrowa4, log2_factor);
    vrowa5 = _mm256_srai_epi16(vrowa5, log2_factor);
    vrowa6 = _mm256_srai_epi16(vrowa6, log2_factor);


    // Calculate 2/4 part of 64 wide lane
    vbeforeshifted = _mm256_slli_epi16(vbefore256b, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    vinterpolate = _mm256_sub_epi16(vbehind256b, vbefore256b);

    __m256i vrowb0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    __m256i vrowb1 = _mm256_add_epi16(vrowb0, vinterpolate);
    __m256i vrowb2 = _mm256_add_epi16(vrowb1, vinterpolate);
    __m256i vrowb3 = _mm256_add_epi16(vrowb2, vinterpolate);
    __m256i vrowb4 = _mm256_add_epi16(vrowb3, vinterpolate);
    __m256i vrowb5 = _mm256_add_epi16(vrowb4, vinterpolate);
    __m256i vrowb6 = _mm256_add_epi16(vrowb5, vinterpolate);

    vrowb0 = _mm256_srai_epi16(vrowb0, log2_factor);
    vrowb1 = _mm256_srai_epi16(vrowb1, log2_factor);
    vrowb2 = _mm256_srai_epi16(vrowb2, log2_factor);
    vrowb3 = _mm256_srai_epi16(vrowb3, log2_factor);
    vrowb4 = _mm256_srai_epi16(vrowb4, log2_factor);
    vrowb5 = _mm256_srai_epi16(vrowb5, log2_factor);
    vrowb6 = _mm256_srai_epi16(vrowb6, log2_factor);


    // Calculate 3/4 part of 64 wide lane
    vbeforeshifted = _mm256_slli_epi16(vbefore256c, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    vinterpolate = _mm256_sub_epi16(vbehind256c, vbefore256c);

    __m256i vrowc0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    __m256i vrowc1 = _mm256_add_epi16(vrowc0, vinterpolate);
    __m256i vrowc2 = _mm256_add_epi16(vrowc1, vinterpolate);
    __m256i vrowc3 = _mm256_add_epi16(vrowc2, vinterpolate);
    __m256i vrowc4 = _mm256_add_epi16(vrowc3, vinterpolate);
    __m256i vrowc5 = _mm256_add_epi16(vrowc4, vinterpolate);
    __m256i vrowc6 = _mm256_add_epi16(vrowc5, vinterpolate);

    vrowc0 = _mm256_srai_epi16(vrowc0, log2_factor);
    vrowc1 = _mm256_srai_epi16(vrowc1, log2_factor);
    vrowc2 = _mm256_srai_epi16(vrowc2, log2_factor);
    vrowc3 = _mm256_srai_epi16(vrowc3, log2_factor);
    vrowc4 = _mm256_srai_epi16(vrowc4, log2_factor);
    vrowc5 = _mm256_srai_epi16(vrowc5, log2_factor);
    vrowc6 = _mm256_srai_epi16(vrowc6, log2_factor);


    // Calculate 3/4 part of 64 wide lane
    vbeforeshifted = _mm256_slli_epi16(vbefore256d, log2_factor);

    // Add rounding offset
    vbeforeshifted = _mm256_add_epi16(vbeforeshifted, vrnd);

    vinterpolate = _mm256_sub_epi16(vbehind256d, vbefore256d);

    __m256i vrowd0 = _mm256_add_epi16(vbeforeshifted, vinterpolate);
    __m256i vrowd1 = _mm256_add_epi16(vrowd0, vinterpolate);
    __m256i vrowd2 = _mm256_add_epi16(vrowd1, vinterpolate);
    __m256i vrowd3 = _mm256_add_epi16(vrowd2, vinterpolate);
    __m256i vrowd4 = _mm256_add_epi16(vrowd3, vinterpolate);
    __m256i vrowd5 = _mm256_add_epi16(vrowd4, vinterpolate);
    __m256i vrowd6 = _mm256_add_epi16(vrowd5, vinterpolate);

    vrowd0 = _mm256_srai_epi16(vrowd0, log2_factor);
    vrowd1 = _mm256_srai_epi16(vrowd1, log2_factor);
    vrowd2 = _mm256_srai_epi16(vrowd2, log2_factor);
    vrowd3 = _mm256_srai_epi16(vrowd3, log2_factor);
    vrowd4 = _mm256_srai_epi16(vrowd4, log2_factor);
    vrowd5 = _mm256_srai_epi16(vrowd5, log2_factor);
    vrowd6 = _mm256_srai_epi16(vrowd6, log2_factor);


    // Store results
    __m256i vres00 = _mm256_packus_epi16(vrowa0, vrowb0);
    __m256i vres01 = _mm256_packus_epi16(vrowc0, vrowd0);
    
    __m256i vres10 = _mm256_packus_epi16(vrowa1, vrowb1);
    __m256i vres11 = _mm256_packus_epi16(vrowc1, vrowd1);

    __m256i vres20 = _mm256_packus_epi16(vrowa2, vrowb2);
    __m256i vres21 = _mm256_packus_epi16(vrowc2, vrowd2);

    __m256i vres30 = _mm256_packus_epi16(vrowa3, vrowb3);
    __m256i vres31 = _mm256_packus_epi16(vrowc3, vrowd3);

    __m256i vres40 = _mm256_packus_epi16(vrowa4, vrowb4);
    __m256i vres41 = _mm256_packus_epi16(vrowc4, vrowd4);

    __m256i vres50 = _mm256_packus_epi16(vrowa5, vrowb5);
    __m256i vres51 = _mm256_packus_epi16(vrowc5, vrowd5);

    __m256i vres60 = _mm256_packus_epi16(vrowa6, vrowb6);
    __m256i vres61 = _mm256_packus_epi16(vrowc6, vrowd6);

    
    vres00 = _mm256_permute4x64_epi64(vres00, _MM_SHUFFLE(3, 1, 2, 0));
    vres01 = _mm256_permute4x64_epi64(vres01, _MM_SHUFFLE(3, 1, 2, 0));
    vres10 = _mm256_permute4x64_epi64(vres10, _MM_SHUFFLE(3, 1, 2, 0));
    vres11 = _mm256_permute4x64_epi64(vres11, _MM_SHUFFLE(3, 1, 2, 0));
    vres20 = _mm256_permute4x64_epi64(vres20, _MM_SHUFFLE(3, 1, 2, 0));
    vres21 = _mm256_permute4x64_epi64(vres21, _MM_SHUFFLE(3, 1, 2, 0));
    vres30 = _mm256_permute4x64_epi64(vres30, _MM_SHUFFLE(3, 1, 2, 0));
    vres31 = _mm256_permute4x64_epi64(vres31, _MM_SHUFFLE(3, 1, 2, 0));
    vres40 = _mm256_permute4x64_epi64(vres40, _MM_SHUFFLE(3, 1, 2, 0));
    vres41 = _mm256_permute4x64_epi64(vres41, _MM_SHUFFLE(3, 1, 2, 0));
    vres50 = _mm256_permute4x64_epi64(vres50, _MM_SHUFFLE(3, 1, 2, 0));
    vres51 = _mm256_permute4x64_epi64(vres51, _MM_SHUFFLE(3, 1, 2, 0));
    vres60 = _mm256_permute4x64_epi64(vres60, _MM_SHUFFLE(3, 1, 2, 0));
    vres61 = _mm256_permute4x64_epi64(vres61, _MM_SHUFFLE(3, 1, 2, 0));


    _mm256_store_si256((__m256i*)(dst + (i * 512) +   0), vres00);
    _mm256_store_si256((__m256i*)(dst + (i * 512) +  32), vres01);
    _mm256_store_si256((__m256i*)(dst + (i * 512) +  64), vres10);
    _mm256_store_si256((__m256i*)(dst + (i * 512) +  96), vres11);
    _mm256_store_si256((__m256i*)(dst + (i * 512) + 128), vres20);
    _mm256_store_si256((__m256i*)(dst + (i * 512) + 160), vres21);
    _mm256_store_si256((__m256i*)(dst + (i * 512) + 192), vres30);
    _mm256_store_si256((__m256i*)(dst + (i * 512) + 224), vres31);
    _mm256_store_si256((__m256i*)(dst + (i * 512) + 256), vres40);
    _mm256_store_si256((__m256i*)(dst + (i * 512) + 288), vres41);
    _mm256_store_si256((__m256i*)(dst + (i * 512) + 320), vres50);
    _mm256_store_si256((__m256i*)(dst + (i * 512) + 352), vres51);
    _mm256_store_si256((__m256i*)(dst + (i * 512) + 384), vres60);
    _mm256_store_si256((__m256i*)(dst + (i * 512) + 416), vres61);
    

    vbefore256a = vbehind256a;
    vbefore256b = vbehind256b;
    vbefore256c = vbehind256c;
    vbefore256d = vbehind256d;
  }
}

static void mip_upsampling_w64_ups8_ver_avx2_alt(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  const uvg_pixel* src_ptr = src;
  const uvg_pixel* dst_ptr = dst;

  const __m256i zeros = _mm256_setzero_si256();
  const __m256i ones = _mm256_set1_epi8(1);
  const __m256i twos = _mm256_set1_epi8(2);
  const __m256i threes = _mm256_set1_epi8(3);
  const __m256i fours = _mm256_set1_epi8(4);
  const __m256i fives = _mm256_set1_epi8(5);
  const __m256i sixes = _mm256_set1_epi8(6);
  const __m256i sevens = _mm256_set1_epi8(7);
  const __m256i eights = _mm256_set1_epi8(8);

  __m256i vbeforeleft = _mm256_load_si256((__m256i*)(ref + 0));
  __m256i vbeforeright = _mm256_load_si256((__m256i*)(ref + 32));

  for (int i = 0; i < 8; ++i) {
    __m256i vbehindleft  = _mm256_load_si256((__m256i*)(src_ptr + 0));
    __m256i vbehindright = _mm256_load_si256((__m256i*)(src_ptr + 32));

    // Calculate left side of 64 wide lane.
    // Calculate the 7 interpolated lines between before and behind. Ordered by number from top to bottom.
    __m256i vleft3 = _mm256_avg_epu8(vbeforeleft, vbehindleft); // Middle
    __m256i vleft1 = _mm256_avg_epu8(vleft3, vbeforeleft);      // Top middle
    __m256i vleft5 = _mm256_avg_epu8(vleft3, vbehindleft);      // Bottom middle
    __m256i vleft0 = _mm256_avg_epu8(vbeforeleft, vleft1);      // Top middle top
    __m256i vleft2 = _mm256_avg_epu8(vleft1, vleft3);           // Top middle bottom
    __m256i vleft4 = _mm256_avg_epu8(vleft3, vleft5);           // Bottom middle top
    __m256i vleft6 = _mm256_avg_epu8(vleft5, vbehindleft);      // Bottom middle bottom

    // Calculate the three and two last bits of difference between before and behind. These bits are used to determine if there will be rounding error.
    __m256i diff = _mm256_sub_epi8(vbehindleft, vbeforeleft);
    diff = _mm256_and_si256(diff, sevens);
    __m256i three_diff = _mm256_and_si256(diff, threes);

    // Bottom side
    __m256i mask = _mm256_cmpgt_epi8(diff, fours);  // The rounding error mask will be generated based on the calculated last bits.
    __m256i sub_amount = _mm256_blendv_epi8(zeros, ones, mask); // If 5, 6, 7 select one
    vleft6 = _mm256_sub_epi8(vleft6, sub_amount);

    mask = _mm256_cmpeq_epi8(three_diff, threes);
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask); // If 3 or 7 select one
    vleft5 = _mm256_sub_epi8(vleft5, sub_amount);

    __m256i is_two = _mm256_cmpeq_epi8(diff, twos);
    __m256i is_five = _mm256_cmpeq_epi8(diff, fives);
    mask = _mm256_or_si256(mask, is_two);
    mask = _mm256_or_si256(mask, is_five);
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask); // If 2, 3, 5, or 7 select one
    vleft4 = _mm256_sub_epi8(vleft4, sub_amount);

    // Top side
    diff = _mm256_blendv_epi8(diff, eights, _mm256_cmpeq_epi8(zeros, diff)); // Replace zeros with eights to enable using GT
    mask = _mm256_cmpgt_epi8(diff, threes);
    sub_amount = _mm256_blendv_epi8(ones, zeros, mask); // If greater than three select zero
    vleft0 = _mm256_sub_epi8(vleft0, sub_amount);

    mask = _mm256_cmpeq_epi8(three_diff, ones);
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask); // If 1 or 5 select one
    vleft1 = _mm256_sub_epi8(vleft1, sub_amount);

    __m256i is_three = _mm256_cmpeq_epi8(diff, threes);
    __m256i is_six = _mm256_cmpeq_epi8(diff, sixes);
    mask = _mm256_or_si256(mask, is_three);
    mask = _mm256_or_si256(mask, is_six);
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask); // If 1, 3, 5, 6 select one
    vleft2 = _mm256_sub_epi8(vleft2, sub_amount);

    
    // Calculate right side of 64 wide lane.
    // Calculate the 7 interpolated lines between before and behind. Ordered by number from top to bottom.
    __m256i vright3 = _mm256_avg_epu8(vbeforeright, vbehindright); // Middle
    __m256i vright1 = _mm256_avg_epu8(vright3, vbeforeright);      // Top middle
    __m256i vright5 = _mm256_avg_epu8(vright3, vbehindright);      // Bottom middle
    __m256i vright0 = _mm256_avg_epu8(vbeforeright, vright1);      // Top middle top
    __m256i vright2 = _mm256_avg_epu8(vright1, vright3);           // Top middle bottom
    __m256i vright4 = _mm256_avg_epu8(vright3, vright5);           // Bottom middle top
    __m256i vright6 = _mm256_avg_epu8(vright5, vbehindright);      // Bottom middle bottom

    // Calculate the three and two last bits of difference between before and behind. These bits are used to determine if there will be rounding error.
    diff = _mm256_sub_epi8(vbehindright, vbeforeright);
    diff = _mm256_and_si256(diff, sevens);
    three_diff = _mm256_and_si256(diff, threes);

    // Bottom side
    mask = _mm256_cmpgt_epi8(diff, fours);  // The rounding error mask will be generated based on the calculated last bits.
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask); // If 5, 6, 7 select one
    vright6 = _mm256_sub_epi8(vright6, sub_amount);

    mask = _mm256_cmpeq_epi8(three_diff, threes);
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask); // If 3 or 7 select one
    vright5 = _mm256_sub_epi8(vright5, sub_amount);

    is_two = _mm256_cmpeq_epi8(diff, twos);
    is_five = _mm256_cmpeq_epi8(diff, fives);
    mask = _mm256_or_si256(mask, is_two);
    mask = _mm256_or_si256(mask, is_five);
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask); // If 2, 3, 5, or 7 select one
    vright4 = _mm256_sub_epi8(vright4, sub_amount);

    // Top side
    diff = _mm256_blendv_epi8(diff, eights, _mm256_cmpeq_epi8(zeros, diff)); // Replace zeros with eights to enable using GT
    mask = _mm256_cmpgt_epi8(diff, threes);
    sub_amount = _mm256_blendv_epi8(ones, zeros, mask); // If greater than three select zero
    vright0 = _mm256_sub_epi8(vright0, sub_amount);

    mask = _mm256_cmpeq_epi8(three_diff, ones);
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask); // If 1 or 5 select one
    vright1 = _mm256_sub_epi8(vright1, sub_amount);

    is_three = _mm256_cmpeq_epi8(diff, threes);
    is_six = _mm256_cmpeq_epi8(diff, sixes);
    mask = _mm256_or_si256(mask, is_three);
    mask = _mm256_or_si256(mask, is_six);
    sub_amount = _mm256_blendv_epi8(zeros, ones, mask); // If 1, 3, 5, 6 select one
    vright2 = _mm256_sub_epi8(vright2, sub_amount);


    // Store results
    _mm256_store_si256((__m256i*)(dst_ptr +   0), vleft0);
    _mm256_store_si256((__m256i*)(dst_ptr +  32), vright0);
    _mm256_store_si256((__m256i*)(dst_ptr +  64), vleft1);
    _mm256_store_si256((__m256i*)(dst_ptr +  96), vright1);
    _mm256_store_si256((__m256i*)(dst_ptr + 128), vleft2);
    _mm256_store_si256((__m256i*)(dst_ptr + 160), vright2);
    _mm256_store_si256((__m256i*)(dst_ptr + 192), vleft3);
    _mm256_store_si256((__m256i*)(dst_ptr + 224), vright3);
    _mm256_store_si256((__m256i*)(dst_ptr + 256), vleft4);
    _mm256_store_si256((__m256i*)(dst_ptr + 288), vright4);
    _mm256_store_si256((__m256i*)(dst_ptr + 320), vleft5);
    _mm256_store_si256((__m256i*)(dst_ptr + 352), vright5);
    _mm256_store_si256((__m256i*)(dst_ptr + 384), vleft6);
    _mm256_store_si256((__m256i*)(dst_ptr + 416), vright6);

    vbeforeleft = vbehindleft;
    vbeforeright = vbehindright;

    dst_ptr += 512;
    src_ptr += 512;
  }
}

/** \brief Matrix weighted intra prediction.
*/
static void mip_predict_avx2(
  //const encoder_state_t* const state,
  const uvg_intra_references* const refs,
  const uint16_t pred_block_width,
  const uint16_t pred_block_height,
  uvg_pixel* dst,
  const int mip_mode,
  const bool mip_transp)
{
  // MIP prediction uses int values instead of uvg_pixel as some temp values may be negative

  //uvg_pixel* out = dst;
  //uvg_pixel result[64 * 64] = { 0 };
  uvg_pixel* result = dst;
  const int mode_idx = mip_mode;

  // *** INPUT PREP ***

  // Initialize prediction parameters START
  uint16_t width = pred_block_width;
  uint16_t height = pred_block_height;

  int size_id; // Prediction block type
  if (width == 4 && height == 4) {
    size_id = 0;
  }
  else if (width == 4 || height == 4 || (width == 8 && height == 8)) {
    size_id = 1;
  }
  else {
    size_id = 2;
  }

  // Reduced boundary and prediction sizes
  int red_bdry_size = (size_id == 0) ? 2 : 4;
  int red_pred_size = (size_id < 2) ? 4 : 8;

  // Upsampling factors
  uint16_t ups_hor_factor = width / red_pred_size;
  uint16_t ups_ver_factor = height / red_pred_size;
  // Initialize prediction parameters END

  const uvg_pixel* ref_samples_top = &refs->ref.top[1];
  const uvg_pixel* ref_samples_left = &refs->ref.left[1];

  // Compute reduced boundary with Haar-downsampling
  const int input_size = 2 * red_bdry_size;

  uvg_pixel red_bdry[MIP_MAX_INPUT_SIZE];
  uvg_pixel red_bdry_trans[MIP_MAX_INPUT_SIZE];
  int16_t red_bdry16[MIP_MAX_INPUT_SIZE];
  int16_t red_bdry_trans16[MIP_MAX_INPUT_SIZE];

  uvg_pixel* const top_reduced = &red_bdry[0];
  uvg_pixel* const left_reduced = &red_bdry[red_bdry_size];

  if (width == 4 && height == 4) {
    // 4 to 2 downsampling for both dimensions
    mip_ref_downsampling_4x4_4to2_avx2(top_reduced, ref_samples_top, ref_samples_left);
  }
  else if (width == 8 && height == 8) {
    // 8 to 4 downsampling for both dimensions
    mip_ref_downsampling_8x8_8to4_avx2(top_reduced, ref_samples_top, ref_samples_left);
  }
  else {
    // Horizontal downsampling
    switch (width) {
      case 4:
        // 4x4 case handled elsewhere.
        // No horizontal downsampling needed. Copy pixels.
        memcpy(top_reduced, ref_samples_top, 4 * sizeof(uvg_pixel));
        break;
      case 8:  mip_ref_downsampling_1D_8to4_avx2(top_reduced, ref_samples_top);  break; // 8x8 case handled elsewhere.
      case 16: mip_ref_downsampling_1D_16to4_avx2(top_reduced, ref_samples_top); break;
      case 32: mip_ref_downsampling_1D_32to4_avx2(top_reduced, ref_samples_top); break;
      case 64: mip_ref_downsampling_1D_64to4_avx2(top_reduced, ref_samples_top); break;
      default:
        assert(false && "MIP horizontal downsampling. Invalid width.\n");
        break;
    }

    // Vertical downsampling
    switch (height) {
      case 4:
        // 4x4 case handled elsewhere.
        // No vertical downsampling needed. Copy pixels.
        memcpy(left_reduced, ref_samples_left, 4 * sizeof(uvg_pixel));
        break;
      case 8:  mip_ref_downsampling_1D_8to4_avx2(left_reduced, ref_samples_left);  break; // 8x8 case handled elsewhere.
      case 16: mip_ref_downsampling_1D_16to4_avx2(left_reduced, ref_samples_left); break;
      case 32: mip_ref_downsampling_1D_32to4_avx2(left_reduced, ref_samples_left); break;
      case 64: mip_ref_downsampling_1D_64to4_avx2(left_reduced, ref_samples_left); break;
      default:
        assert(false && "MIP vertical downsampling. Invalid height.\n");
        break;
    }
  }
  

  // Transposed reduced boundaries
  uvg_pixel* const left_reduced_trans = &red_bdry_trans[0];
  uvg_pixel* const top_reduced_trans = &red_bdry_trans[red_bdry_size];

  for (int x = 0; x < red_bdry_size; x++) {
    top_reduced_trans[x] = top_reduced[x];
  }
  for (int y = 0; y < red_bdry_size; y++) {
    left_reduced_trans[y] = left_reduced[y];
  }

  uvg_pixel input_offset = red_bdry[0];
  uvg_pixel input_offset_trans = red_bdry_trans[0];

  const bool has_first_col = (size_id < 2);
  // First column of matrix not needed for large blocks
  // These can potentially fail with uvg_pixel
  red_bdry16[0] = has_first_col ? ((1 << (UVG_BIT_DEPTH - 1)) - input_offset) : 0;
  red_bdry_trans16[0] = has_first_col ? ((1 << (UVG_BIT_DEPTH - 1)) - input_offset_trans) : 0;

  // This fails with uvg_pixel, here at least int16_t is needed
  for (int i = 1; i < input_size; ++i) {
    red_bdry16[i] = red_bdry[i] - input_offset;
    red_bdry_trans16[i] = red_bdry_trans[i] - input_offset_trans;
  }

  // *** INPUT PREP *** END

  // *** BLOCK PREDICT ***

  const bool need_upsampling = (ups_hor_factor > 1) || (ups_ver_factor > 1);
  const bool transpose = mip_transp;

  const uint8_t* matrix = 0;
  const uint16_t* matrix16 = 0;
  switch (size_id) {
  case 0:
    matrix16 = &uvg_mip_sid0_weights[mode_idx][0][0];
    break;
  case 1:
    matrix16 = &uvg_mip_sid1_weights[mode_idx * 128];
    break;
  case 2:
    //matrix = &uvg_mip_matrix_16x16[mode_idx][0][0];
    matrix16 = &uvg_mip_sid2_weights[mode_idx * 512];
    break;
  default:
    assert(false && "Invalid MIP size id.");
  }

  // Max possible size is red_pred_size * red_pred_size, red_pred_size can be either 4 or 8
  uvg_pixel red_pred_buffer[8 * 8];
  uvg_pixel* const reduced_pred = need_upsampling ? red_pred_buffer : result;

  const uvg_pixel* const reduced_bdry = transpose ? red_bdry_trans : red_bdry;
  const int16_t* const reduced_bdry16 = transpose ? red_bdry_trans16 : red_bdry16;

  switch (size_id) {
    case 0: mip_reduced_pred_sid0_avx2(reduced_pred, reduced_bdry16, matrix16, transpose, input_offset, input_offset_trans); break;
    case 1: mip_reduced_pred_sid1_avx2(reduced_pred, reduced_bdry16, matrix16, transpose, input_offset, input_offset_trans); break;
    case 2: mip_reduced_pred_sid2_avx2(reduced_pred, reduced_bdry16, matrix16, transpose, input_offset, input_offset_trans); break;
    default:
      assert(false && "Intra MIP: invalid size id.\n");
      break;
  }
  if (need_upsampling) {
    const uvg_pixel* ver_src = reduced_pred;
    uint16_t ver_src_step = width;

    if (ups_hor_factor > 1) {
      uvg_pixel* const hor_dst = result + (ups_ver_factor - 1) * width;
      ver_src = hor_dst;
      ver_src_step *= ups_ver_factor;

      switch (width) {
        // Case 4 does not exist. There is no need for horizontal upsampling when width is 4.
        case 8: 
          // This will only get called for 8x8 blocks.
          mip_upsampling_w8_ups2_hor_avx2_alt(hor_dst, reduced_pred, ref_samples_left, ver_src_step, ups_ver_factor);
          break;
        case 16: 
          if (red_pred_size == 4) {
            mip_upsampling_w16_ups4_hor_avx2(hor_dst, reduced_pred, ref_samples_left, ver_src_step, ups_ver_factor);
          }
          else {
            mip_upsampling_w16_ups2_hor_avx2(hor_dst, reduced_pred, ref_samples_left, ver_src_step, ups_ver_factor); // Works for height 8, 16, 32 and 64. Upsamples 1 to 2.
          }
          break;
        case 32:
          if (red_pred_size == 4) {
            mip_upsampling_w32_ups8_hor_avx2(hor_dst, reduced_pred, ref_samples_left, ver_src_step, ups_ver_factor);
          }
          else {
            mip_upsampling_w32_ups4_hor_avx2_alt(hor_dst, reduced_pred, ref_samples_left, ver_src_step, ups_ver_factor); // Works for height 8, 16, 32 and 64. Upsamples 1 to 4.
          }
          break;
        case 64: 
          mip_upsampling_w64_ups8_hor_avx2(hor_dst, reduced_pred, ref_samples_left, ver_src_step, ups_ver_factor); // Works for height 8, 16, 32 and 64. Upsamples 1 to 8.
          break;
        default:
          assert(false && "Invalid MIP width.\n");
          break;
      }
    }

    //uvg_pixel tmp[64 * 64] = {0};
    if (ups_ver_factor > 1) {
      switch (width) {
        case 4: 
          if (ups_ver_factor == 2) {
            mip_upsampling_w4_ups2_ver_avx2(result, ver_src, ref_samples_top);
          }
          else if (ups_ver_factor == 4) {
            mip_upsampling_w4_ups4_ver_avx2(result, ver_src, ref_samples_top);
          }
          else {
            mip_upsampling_w4_ups8_ver_avx2(result, ver_src, ref_samples_top);
          }
          break;
        
        case 8: 
          if (ups_ver_factor == 2) {
            if (height == 8) {
              mip_upsampling_w8_ups2_h8_ver_avx2(result, ver_src, ref_samples_top);
            }
            else { // Height == 16
              mip_upsampling_w8_ups2_h16_ver_avx2(result, ver_src, ref_samples_top);
            }
          }
          else if (ups_ver_factor == 4) {
            mip_upsampling_w8_ups4_ver_avx2(result, ver_src, ref_samples_top);
          }
          else {
            mip_upsampling_w8_ups8_ver_avx2(result, ver_src, ref_samples_top);
          }
          break;

        case 16: 
          if (ups_ver_factor == 2) {
            mip_upsampling_w16_ups2_ver_avx2(result, ver_src, ref_samples_top);
          }
          else if (ups_ver_factor == 4) {
            mip_upsampling_w16_ups4_ver_avx2(result, ver_src, ref_samples_top);
          }
          else {
            mip_upsampling_w16_ups8_ver_avx2(result, ver_src, ref_samples_top);
          }
          break;

        case 32: 
          if (ups_ver_factor == 2) {
            mip_upsampling_w32_ups2_ver_avx2(result, ver_src, ref_samples_top);
          }
          else if (ups_ver_factor == 4) {
            mip_upsampling_w32_ups4_ver_avx2_alt(result, ver_src, ref_samples_top);
          }
          else {
            mip_upsampling_w32_ups8_ver_avx2_alt(result, ver_src, ref_samples_top);
          }
          break;
          
        case 64: 
          if (ups_ver_factor == 2) {
            mip_upsampling_w64_ups2_ver_avx2(result, ver_src, ref_samples_top);
          }
          else if (ups_ver_factor == 4) {
            mip_upsampling_w64_ups4_ver_avx2_alt(result, ver_src, ref_samples_top);
          }
          else {
            mip_upsampling_w64_ups8_ver_avx2_alt(result, ver_src, ref_samples_top);
          }
          break;

        default:
          assert(false && "Invalid MIP width.\n");
          break;
      }
    }
  }
  // *** BLOCK PREDICT *** END
}


#endif // UVG_BIT_DEPTH == 8

#endif // COMPILE_INTEL_AVX2 && defined X86_64

int uvg_strategy_register_intra_avx2(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if COMPILE_INTEL_AVX2 && defined X86_64
#if UVG_BIT_DEPTH == 8
  if (bitdepth == 8) {
    success &= uvg_strategyselector_register(opaque, "angular_pred", "avx2", 40, &uvg_angular_pred_avx2);
    success &= uvg_strategyselector_register(opaque, "intra_pred_planar", "avx2", 40, &uvg_intra_pred_planar_avx2);
    success &= uvg_strategyselector_register(opaque, "intra_pred_filtered_dc", "avx2", 40, &uvg_intra_pred_filtered_dc_avx2);
    success &= uvg_strategyselector_register(opaque, "pdpc_planar_dc", "avx2", 40, &uvg_pdpc_planar_dc_avx2);
    success &= uvg_strategyselector_register(opaque, "mip_predict", "avx2", 40, &mip_predict_avx2);
  }
#endif //UVG_BIT_DEPTH == 8
#endif //COMPILE_INTEL_AVX2 && defined X86_64
  return success;
}

