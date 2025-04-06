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


ALIGNED(32) static const int16_t cubic_filter[32][4] =
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
ALIGNED(32) static const int8_t cubic_filter_8bit_c[32][4] =
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
ALIGNED(32) static const int8_t cubic_filter_8bit_g[32][4] = 
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


static void angular_pred_w4_ver_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t* delta_int, const int16_t* delta_fract, const int height, const int8_t(*filter)[4])
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
    0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01,
    0x04, 0x05, 0x04, 0x05, 0x04, 0x05, 0x04, 0x05,
    0x08, 0x09, 0x08, 0x09, 0x08, 0x09, 0x08, 0x09,
    0x0c, 0x0d, 0x0c, 0x0d, 0x0c, 0x0d, 0x0c, 0x0d
  );

  const __m256i w_shuf_23 = _mm256_setr_epi8(
    0x02, 0x03, 0x02, 0x03, 0x02, 0x03, 0x02, 0x03,
    0x06, 0x07, 0x06, 0x07, 0x06, 0x07, 0x06, 0x07,
    0x0a, 0x0b, 0x0a, 0x0b, 0x0a, 0x0b, 0x0a, 0x0b,
    0x0e, 0x0f, 0x0e, 0x0f, 0x0e, 0x0f, 0x0e, 0x0f
  );

  // Do 4-tap intra interpolation filtering
  // For a 4 width block, height must be at least 4. Handle 4 lines at once
  for (int y = 0; y < height; y += 4) {
    // Copy the filter to local memory
    __m128i vdfract = _mm_loadu_si128((__m128i*)&delta_fract[y]);
    __m128i vidxw = _mm_cvtepi16_epi32(vdfract);
    __m128i all_weights = _mm_i32gather_epi32((const int32_t*)filter, vidxw, 4);

    __m256i weights256 = _mm256_insertf128_si256(_mm256_castsi128_si256(all_weights), all_weights, 1);

    // Shuffle the interpolation weights into place.
    __m256i w01 = _mm256_shuffle_epi8(weights256, w_shuf_01);
    __m256i w23 = _mm256_shuffle_epi8(weights256, w_shuf_23);

    // This solution assumes the delta int values to be 64-bit
    // Cast from 16-bit to 64-bit.
    __m128i vdelta_int = _mm_loadu_si128((__m128i*)&delta_int[y]);
    __m256i vidx = _mm256_cvtepi16_epi64(vdelta_int);

    __m256i vp = _mm256_i64gather_epi64((const long long int*)ref_main, vidx, 1);
    __m256i vp_01 = _mm256_shuffle_epi8(vp, p_shuf_01);
    __m256i vp_23 = _mm256_shuffle_epi8(vp, p_shuf_23);

    __m256i dot_01 = _mm256_maddubs_epi16(vp_01, w01);
    __m256i dot_23 = _mm256_maddubs_epi16(vp_23, w23);
    __m256i sum = _mm256_add_epi16(dot_01, dot_23);
    sum = _mm256_add_epi16(sum, _mm256_set1_epi16(32));
    sum = _mm256_srai_epi16(sum, 6);

    __m128i lo = _mm256_castsi256_si128(sum);
    __m128i hi = _mm256_extracti128_si256(sum, 1);
    __m128i packed = _mm_packus_epi16(lo, hi);

    _mm_store_si128((__m128i*)dst, packed);
    dst += 16;
  }
}

static void angular_pred_w8_h2_ver_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t* delta_int, const int16_t* delta_fract, const int height, const int8_t(*filter)[4])
{
  //const int width = 8;

  const __m128i p_shuf_01 = _mm_setr_epi8(
    0x00, 0x01, 0x01, 0x02, 0x02, 0x03, 0x03, 0x04,
    0x04, 0x05, 0x05, 0x06, 0x06, 0x07, 0x07, 0x08
  );

  const __m128i p_shuf_23 = _mm_setr_epi8(
    0x02, 0x03, 0x03, 0x04, 0x04, 0x05, 0x05, 0x06,
    0x06, 0x07, 0x07, 0x08, 0x08, 0x09, 0x09, 0x0a
  );

  const __m256i w_shuf_01 = _mm256_setr_epi8(
    0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01,
    0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01,
    0x04, 0x05, 0x04, 0x05, 0x04, 0x05, 0x04, 0x05,
    0x04, 0x05, 0x04, 0x05, 0x04, 0x05, 0x04, 0x05
  );

  const __m256i w_shuf_23 = _mm256_setr_epi8(
    0x02, 0x03, 0x02, 0x03, 0x02, 0x03, 0x02, 0x03,
    0x02, 0x03, 0x02, 0x03, 0x02, 0x03, 0x02, 0x03,
    0x06, 0x07, 0x06, 0x07, 0x06, 0x07, 0x06, 0x07,
    0x06, 0x07, 0x06, 0x07, 0x06, 0x07, 0x06, 0x07
  );

  // Do 4-tap intra interpolation filtering
  // For a 8 width block, height must be at least 2. Handle 2 lines at once
  for (int y = 0; y < height; y += 2) {
    
    // Load and shuffle filter weights
    __m128i vidxw = _mm_loadu_si128((__m128i*)&delta_fract[y]);
    __m128i vidxw32 = _mm_cvtepi16_epi32(vidxw);
    __m128i all_weights = _mm_i32gather_epi32((const int32_t*)filter, vidxw32, 4);
    __m256i aw256 = _mm256_inserti128_si256(_mm256_castsi128_si256(all_weights), all_weights, 1);

    __m256i w01 = _mm256_shuffle_epi8(aw256, w_shuf_01);
    __m256i w23 = _mm256_shuffle_epi8(aw256, w_shuf_23);

    // Load and shuffle reference pixels
    __m128i vp0 = _mm_loadu_si128((__m128i*)(ref_main + delta_int[y + 0]));
    __m128i vp1 = _mm_loadu_si128((__m128i*)(ref_main + delta_int[y + 1]));

    __m256i vp_01 = _mm256_castsi128_si256(_mm_shuffle_epi8(vp0, p_shuf_01));
    vp_01 = _mm256_inserti128_si256(vp_01, _mm_shuffle_epi8(vp1, p_shuf_01), 1);

    __m256i vp_23 = _mm256_castsi128_si256(_mm_shuffle_epi8(vp0, p_shuf_23));
    vp_23 = _mm256_inserti128_si256(vp_23, _mm_shuffle_epi8(vp1, p_shuf_23), 1);

    __m256i vmadd01 = _mm256_maddubs_epi16(vp_01, w01);
    __m256i vmadd23 = _mm256_maddubs_epi16(vp_23, w23);
    __m256i sum = _mm256_add_epi16(vmadd01, vmadd23);
    sum = _mm256_add_epi16(sum, _mm256_set1_epi16(32));
    sum = _mm256_srai_epi16(sum, 6);

    __m128i lo = _mm256_castsi256_si128(sum);
    __m128i hi = _mm256_extracti128_si256(sum, 1);
    __m128i packed = _mm_packus_epi16(lo, hi);

    _mm_store_si128((__m128i*)dst, packed);
    dst += 16;
  }
}

static void angular_pred_w8_ver_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t* delta_int, const int16_t* delta_fract, const int height, const int8_t(*filter)[4])
{
  //const int width = 8;

  const __m128i p_shuf_01 = _mm_setr_epi8(
    0x00, 0x01, 0x01, 0x02, 0x02, 0x03, 0x03, 0x04,
    0x04, 0x05, 0x05, 0x06, 0x06, 0x07, 0x07, 0x08
  );

  const __m128i p_shuf_23 = _mm_setr_epi8(
    0x02, 0x03, 0x03, 0x04, 0x04, 0x05, 0x05, 0x06,
    0x06, 0x07, 0x07, 0x08, 0x08, 0x09, 0x09, 0x0a
  );

  const __m256i w_shuf_01_row01 = _mm256_setr_epi8(
    0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01,
    0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01,
    0x04, 0x05, 0x04, 0x05, 0x04, 0x05, 0x04, 0x05,
    0x04, 0x05, 0x04, 0x05, 0x04, 0x05, 0x04, 0x05
  );

  const __m256i w_shuf_23_row01 = _mm256_setr_epi8(
    0x02, 0x03, 0x02, 0x03, 0x02, 0x03, 0x02, 0x03,
    0x02, 0x03, 0x02, 0x03, 0x02, 0x03, 0x02, 0x03,
    0x06, 0x07, 0x06, 0x07, 0x06, 0x07, 0x06, 0x07,
    0x06, 0x07, 0x06, 0x07, 0x06, 0x07, 0x06, 0x07
  );

  const __m256i w_shuf_01_row23 = _mm256_setr_epi8(
    0x08, 0x09, 0x08, 0x09, 0x08, 0x09, 0x08, 0x09,
    0x08, 0x09, 0x08, 0x09, 0x08, 0x09, 0x08, 0x09,
    0x0c, 0x0d, 0x0c, 0x0d, 0x0c, 0x0d, 0x0c, 0x0d,
    0x0c, 0x0d, 0x0c, 0x0d, 0x0c, 0x0d, 0x0c, 0x0d
  );

  const __m256i w_shuf_23_row23 = _mm256_setr_epi8(
    0x0a, 0x0b, 0x0a, 0x0b, 0x0a, 0x0b, 0x0a, 0x0b,
    0x0a, 0x0b, 0x0a, 0x0b, 0x0a, 0x0b, 0x0a, 0x0b,
    0x0e, 0x0f, 0x0e, 0x0f, 0x0e, 0x0f, 0x0e, 0x0f,
    0x0e, 0x0f, 0x0e, 0x0f, 0x0e, 0x0f, 0x0e, 0x0f
  );

  // Do 4-tap intra interpolation filtering
  // For a 8 width block, height must be at least 2. This version handles 4 lines at once to minimize vidx loads.
  // No need to check height 2 cases, other function handles that.
  for (int y = 0; y < height; y += 4) {

    // Load and shuffle filter weights
    __m128i vidxw = _mm_loadu_si128((__m128i*) &delta_fract[y]);
    __m128i vidxw32 = _mm_cvtepi16_epi32(vidxw);
    __m128i all_weights = _mm_i32gather_epi32((const int32_t*)filter, vidxw32, 4);
    __m256i aw256 = _mm256_inserti128_si256(_mm256_castsi128_si256(all_weights), all_weights, 1);

    __m256i w01_row01 = _mm256_shuffle_epi8(aw256, w_shuf_01_row01);
    __m256i w23_row01 = _mm256_shuffle_epi8(aw256, w_shuf_23_row01);
    __m256i w01_row23 = _mm256_shuffle_epi8(aw256, w_shuf_01_row23);
    __m256i w23_row23 = _mm256_shuffle_epi8(aw256, w_shuf_23_row23);

    // Load and shuffle reference pixels
    __m128i vp0 = _mm_loadu_si128((__m128i*)(ref_main + delta_int[y + 0]));
    __m128i vp1 = _mm_loadu_si128((__m128i*)(ref_main + delta_int[y + 1]));
    __m128i vp2 = _mm_loadu_si128((__m128i*)(ref_main + delta_int[y + 2]));
    __m128i vp3 = _mm_loadu_si128((__m128i*)(ref_main + delta_int[y + 3]));

    __m256i vp_01_row01 = _mm256_castsi128_si256(_mm_shuffle_epi8(vp0, p_shuf_01));
    vp_01_row01 = _mm256_inserti128_si256(vp_01_row01, _mm_shuffle_epi8(vp1, p_shuf_01), 1);

    __m256i vp_23_row01 = _mm256_castsi128_si256(_mm_shuffle_epi8(vp0, p_shuf_23));
    vp_23_row01 = _mm256_inserti128_si256(vp_23_row01, _mm_shuffle_epi8(vp1, p_shuf_23), 1);

    __m256i vp_01_row23 = _mm256_castsi128_si256(_mm_shuffle_epi8(vp2, p_shuf_01));
    vp_01_row23 = _mm256_inserti128_si256(vp_01_row23, _mm_shuffle_epi8(vp3, p_shuf_01), 1);

    __m256i vp_23_row23 = _mm256_castsi128_si256(_mm_shuffle_epi8(vp2, p_shuf_23));
    vp_23_row23 = _mm256_inserti128_si256(vp_23_row23, _mm_shuffle_epi8(vp3, p_shuf_23), 1);

    __m256i vmadd01_row01 = _mm256_maddubs_epi16(vp_01_row01, w01_row01);
    __m256i vmadd23_row01 = _mm256_maddubs_epi16(vp_23_row01, w23_row01);
    __m256i vmadd01_row23 = _mm256_maddubs_epi16(vp_01_row23, w01_row23);
    __m256i vmadd23_row23 = _mm256_maddubs_epi16(vp_23_row23, w23_row23);


    __m256i sum01 = _mm256_add_epi16(vmadd01_row01, vmadd23_row01);
    __m256i sum23 = _mm256_add_epi16(vmadd01_row23, vmadd23_row23);
    sum01 = _mm256_add_epi16(sum01, _mm256_set1_epi16(32));
    sum23 = _mm256_add_epi16(sum23, _mm256_set1_epi16(32));
    sum01 = _mm256_srai_epi16(sum01, 6);
    sum23 = _mm256_srai_epi16(sum23, 6);

    __m128i lo01 = _mm256_castsi256_si128(sum01);
    __m128i hi01 = _mm256_extracti128_si256(sum01, 1);
    __m128i lo23 = _mm256_castsi256_si128(sum23);
    __m128i hi23 = _mm256_extracti128_si256(sum23, 1);

    __m128i packed01 = _mm_packus_epi16(lo01, hi01);
    __m128i packed23 = _mm_packus_epi16(lo23, hi23);
    //__m256i packed = _mm256_inserti128_si256(_mm256_castsi128_si256(packed01), packed23, 1);

    //_mm256_store_si256((__m256i*)dst, packed);
    _mm_store_si128((__m128i*)(dst + 0), packed01);
    _mm_store_si128((__m128i*)(dst + 16), packed23);
    dst += 32;
  }
}

NO_ASAN
static void angular_pred_w16_ver_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t* delta_int, const int16_t* delta_fract, const int width, const int height, const int8_t(*filter)[4])
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
    0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01,
    0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01,
    0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01,
    0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01
  );

  const __m256i w_shuf_23 = _mm256_setr_epi8(
    0x02, 0x03, 0x02, 0x03, 0x02, 0x03, 0x02, 0x03,
    0x02, 0x03, 0x02, 0x03, 0x02, 0x03, 0x02, 0x03,
    0x02, 0x03, 0x02, 0x03, 0x02, 0x03, 0x02, 0x03,
    0x02, 0x03, 0x02, 0x03, 0x02, 0x03, 0x02, 0x03
  );

  // Do 4-tap intra interpolation filtering
  // For a 16 width block, height can be 1.
  for (int y = 0; y < height; ++y) {

    // Load and shuffle filter weights
    // This load can read beyond the end of the filter table, however the values
    // are not used in the shuffle operation.
    __m128i vweights = _mm_loadu_si128((__m128i*)&filter[delta_fract[y]]);
    __m256i vw256 = _mm256_inserti128_si256(_mm256_castsi128_si256(vweights), vweights, 1);

    __m256i w01 = _mm256_shuffle_epi8(vw256, w_shuf_01);
    __m256i w23 = _mm256_shuffle_epi8(vw256, w_shuf_23);

    for (int x = 0; x < width; x += 16) {
      __m256i vp = _mm256_loadu_si256((__m256i*)(ref_main + x + delta_int[y]));

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


static void angular_pred_w4_hor_high_angle_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t* delta_int, const int16_t* delta_fract, const int height, const int8_t(*filter)[4])
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
    0x00, 0x01, 0x04, 0x05, 0x08, 0x09, 0x0c, 0x0d,
    0x00, 0x01, 0x04, 0x05, 0x08, 0x09, 0x0c, 0x0d,
    0x00, 0x01, 0x04, 0x05, 0x08, 0x09, 0x0c, 0x0d,
    0x00, 0x01, 0x04, 0x05, 0x08, 0x09, 0x0c, 0x0d
  );

  const __m256i w_shuf_23 = _mm256_setr_epi8(
    0x02, 0x03, 0x06, 0x07, 0x0a, 0x0b, 0x0e, 0x0f,
    0x02, 0x03, 0x06, 0x07, 0x0a, 0x0b, 0x0e, 0x0f,
    0x02, 0x03, 0x06, 0x07, 0x0a, 0x0b, 0x0e, 0x0f,
    0x02, 0x03, 0x06, 0x07, 0x0a, 0x0b, 0x0e, 0x0f
  );

  // Copy the filter to local memory
  __m128i vdfract = _mm_load_si128((__m128i*)delta_fract);
  __m128i vidx = _mm_cvtepi16_epi32(vdfract);
  __m128i all_weights = _mm_i32gather_epi32((const int32_t*)filter, vidx, 4);

  __m256i weights256 = _mm256_insertf128_si256(_mm256_castsi128_si256(all_weights), all_weights, 1);
  // Shuffle the interpolation weights into place.
  __m256i w01 = _mm256_shuffle_epi8(weights256, w_shuf_01);
  __m256i w23 = _mm256_shuffle_epi8(weights256, w_shuf_23);

  // For a 4 width block, height must be at least 4. Handle 4 lines at once
  for (int y = 0; y < height; y += 4) {
    // This solution assumes the delta int values to be 64-bit
    // Cast from 16-bit to 64-bit.
    __m128i vidx = _mm_loadu_si128((__m128i*)delta_int);
    __m256i vidx256 = _mm256_cvtepu16_epi64(vidx);
    
    __m256i vp = _mm256_i64gather_epi64((const long long int*)&ref_main[y], vidx256, 1);

    __m256i vp_01 = _mm256_shuffle_epi8(vp, p_shuf_01);
    __m256i vp_23 = _mm256_shuffle_epi8(vp, p_shuf_23);

    vp_01 = _mm256_permute4x64_epi64(vp_01, _MM_SHUFFLE(3, 1, 2, 0));
    vp_23 = _mm256_permute4x64_epi64(vp_23, _MM_SHUFFLE(3, 1, 2, 0));
    vp_01 = _mm256_shuffle_epi32(vp_01, _MM_SHUFFLE(3, 1, 2, 0));
    vp_23 = _mm256_shuffle_epi32(vp_23, _MM_SHUFFLE(3, 1, 2, 0));

    __m256i vmadd01 = _mm256_maddubs_epi16(vp_01, w01);
    __m256i vmadd23 = _mm256_maddubs_epi16(vp_23, w23);
    __m256i sum = _mm256_add_epi16(vmadd01, vmadd23);
    sum = _mm256_add_epi16(sum, _mm256_set1_epi16(32));
    sum = _mm256_srai_epi16(sum, 6);

    __m128i lo = _mm256_castsi256_si128(sum);
    __m128i hi = _mm256_extracti128_si256(sum, 1);
    __m128i packed = _mm_packus_epi16(lo, hi);

    _mm_store_si128((__m128i*)dst, packed);
    dst += 16;
  }
}

static void angular_pred_w4_hor_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t pred_mode, const int16_t multi_ref_line, const int16_t* delta_int, const int16_t* delta_fract, const int height, const int8_t (*filter)[4])
{
  // const int width = 4;

  const __m256i w_shuf_01 = _mm256_setr_epi8(
    0x00, 0x01, 0x04, 0x05, 0x08, 0x09, 0x0c, 0x0d,
    0x00, 0x01, 0x04, 0x05, 0x08, 0x09, 0x0c, 0x0d,
    0x00, 0x01, 0x04, 0x05, 0x08, 0x09, 0x0c, 0x0d,
    0x00, 0x01, 0x04, 0x05, 0x08, 0x09, 0x0c, 0x0d
  );

  const __m256i w_shuf_23 = _mm256_setr_epi8(
    0x02, 0x03, 0x06, 0x07, 0x0a, 0x0b, 0x0e, 0x0f,
    0x02, 0x03, 0x06, 0x07, 0x0a, 0x0b, 0x0e, 0x0f,
    0x02, 0x03, 0x06, 0x07, 0x0a, 0x0b, 0x0e, 0x0f,
    0x02, 0x03, 0x06, 0x07, 0x0a, 0x0b, 0x0e, 0x0f
  );

  const int mode_idx = pred_mode <= 34 ? pred_mode + 12 : 80 - pred_mode; // Considers also wide angle modes.
  const int table_offset = mode_idx * 192 + multi_ref_line * 64;

  const __m256i vpshuf0 = _mm256_load_si256((__m256i*) &intra_luma_interpolation_shuffle_vectors_w4_hor[table_offset + 0]);
  const __m256i vpshuf1 = _mm256_load_si256((__m256i*) &intra_luma_interpolation_shuffle_vectors_w4_hor[table_offset + 32]);
  
  int ref_offset = MIN(delta_int[0], delta_int[3]);

  // Copy the filter to local memory
  __m128i vdfract = _mm_loadu_si128((__m128i*)delta_fract);
  __m128i vidx = _mm_cvtepi16_epi32(vdfract);
  __m128i all_weights = _mm_i32gather_epi32((const int32_t*)filter, vidx, 4);
  
  __m256i weights256 = _mm256_insertf128_si256(_mm256_castsi128_si256(all_weights), all_weights, 1);

  // Shuffle the interpolation weights into place.
  __m256i w01 = _mm256_shuffle_epi8(weights256, w_shuf_01);
  __m256i w23 = _mm256_shuffle_epi8(weights256, w_shuf_23);

  // 4-tap interpolation filtering.
  // For a 4 width block, height must be at least 4. Handle 4 lines at once
  for (int y = 0; y < height; y += 4) {
    // Load 16 samples and shuffle into place
    __m128i vref = _mm_loadu_si128((__m128i*)&ref_main[y + ref_offset]);
    __m256i vp = _mm256_insertf128_si256(_mm256_castsi128_si256(vref), vref, 1);
      
    __m256i vp_01 = _mm256_shuffle_epi8(vp, vpshuf0);
    __m256i vp_23 = _mm256_shuffle_epi8(vp, vpshuf1);

    __m256i vmadd01 = _mm256_maddubs_epi16(vp_01, w01);
    __m256i vmadd23 = _mm256_maddubs_epi16(vp_23, w23);
    __m256i sum = _mm256_add_epi16(vmadd01, vmadd23);
    sum = _mm256_add_epi16(sum, _mm256_set1_epi16(32));
    sum = _mm256_srai_epi16(sum, 6);

    __m128i lo = _mm256_castsi256_si128(sum);
    __m128i hi = _mm256_extracti128_si256(sum, 1);
    __m128i packed = _mm_packus_epi16(lo, hi);

    _mm_store_si128((__m128i*)dst, packed);
    dst += 16;
  }
}

static void angular_pred_w8_hor_high_angle_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t* delta_int, const int16_t* delta_fract, const int height, const int8_t(*filter)[4])
{
  const int width = 8;

  __m128i tmp = _mm_loadu_si128((__m128i*)delta_int);
  __m256i vidx = _mm256_cvtepi16_epi32(tmp);
  // Load weights
  tmp = _mm_load_si128((__m128i*)delta_fract);
  __m256i vidxw = _mm256_cvtepi16_epi32(tmp);
  __m256i vweights = _mm256_i32gather_epi32((const int32_t*)filter, vidxw, 4);

  for (int y = 0; y < height; y += 2) {

    // Do 4-tap intra interpolation filtering
    uvg_pixel* p = (uvg_pixel*)(ref_main + y);
    __m256i vp0 = _mm256_i32gather_epi32((const int*)&ref_main[y + 0], vidx, 1);
    __m256i vp1 = _mm256_i32gather_epi32((const int*)&ref_main[y + 1], vidx, 1);

    __m256i vmadd0 = _mm256_maddubs_epi16(vp0, vweights);
    __m256i vmadd1 = _mm256_maddubs_epi16(vp1, vweights);
    __m256i sum = _mm256_hadd_epi16(vmadd0, vmadd1);
    sum = _mm256_add_epi16(sum, _mm256_set1_epi16(32));
    sum = _mm256_srai_epi16(sum, 6);

    __m128i lo = _mm256_castsi256_si128(sum);
    __m128i hi = _mm256_extracti128_si256(sum, 1);
    __m128i packed = _mm_packus_epi16(lo, hi);
    packed = _mm_shuffle_epi32(packed, _MM_SHUFFLE(3, 1, 2, 0));
      
    _mm_store_si128((__m128i*)dst, packed);
 
    dst += 16;
  }
}

static void angular_pred_w8_hor_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t pred_mode, const int16_t multi_ref_line, const int16_t* delta_int, const int16_t* delta_fract, const int height, const int8_t (*filter)[4])
{
  // const int width = 8;
  
  __m256i vwshuf01 = _mm256_setr_epi8(
    0x00, 0x01, 0x04, 0x05, 0x08, 0x09, 0x0c, 0x0d,
    0x00, 0x01, 0x04, 0x05, 0x08, 0x09, 0x0c, 0x0d,
    0x00, 0x01, 0x04, 0x05, 0x08, 0x09, 0x0c, 0x0d,
    0x00, 0x01, 0x04, 0x05, 0x08, 0x09, 0x0c, 0x0d
  );

  __m256i vwshuf23 = _mm256_setr_epi8(
    0x02, 0x03, 0x06, 0x07, 0x0a, 0x0b, 0x0e, 0x0f,
    0x02, 0x03, 0x06, 0x07, 0x0a, 0x0b, 0x0e, 0x0f,
    0x02, 0x03, 0x06, 0x07, 0x0a, 0x0b, 0x0e, 0x0f,
    0x02, 0x03, 0x06, 0x07, 0x0a, 0x0b, 0x0e, 0x0f
  );

  int ref_offset = MIN(delta_int[0], delta_int[7]);
  const __m256i v32s = _mm256_set1_epi16(32);

  // Load weights
  __m128i tmp = _mm_loadu_si128((__m128i*)delta_fract);
  __m256i vidxw = _mm256_cvtepi16_epi32(tmp);
  __m256i vweights = _mm256_i32gather_epi32((const int32_t*)filter, vidxw, 4);
  
  __m256i vw01 = _mm256_shuffle_epi8(vweights, vwshuf01);
  __m256i vw23 = _mm256_shuffle_epi8(vweights, vwshuf23);

  vw01 = _mm256_permute4x64_epi64(vw01, _MM_SHUFFLE(3, 1, 2, 0));
  vw23 = _mm256_permute4x64_epi64(vw23, _MM_SHUFFLE(3, 1, 2, 0));

  const int mode_idx = pred_mode <= 34 ? pred_mode + 12 : 80 - pred_mode; // Considers also wide angle modes.
  const int table_offset = mode_idx * 192 + multi_ref_line * 64;

  const __m256i vpshuf01 = _mm256_loadu_si256((__m256i*) &intra_luma_interpolation_shuffle_vectors_w8_hor[table_offset + 0]);
  const __m256i vpshuf23 = _mm256_loadu_si256((__m256i*) &intra_luma_interpolation_shuffle_vectors_w8_hor[table_offset + 32]);

  // 4-tap interpolation filtering.
  // For a 8 width block, height must be at least 2. Handle 2 lines at once.
  for (int y = 0; y < height; y += 2) {
    // Load samples and shuffle into place
    __m128i vp = _mm_loadu_si128((__m128i*)&ref_main[y + ref_offset]);
    __m256i vp256 = _mm256_inserti128_si256(_mm256_castsi128_si256(vp), vp, 1);

    __m256i vp01 = _mm256_shuffle_epi8(vp256, vpshuf01);
    __m256i vp23 = _mm256_shuffle_epi8(vp256, vpshuf23);
    
    __m256i vmadd01 = _mm256_maddubs_epi16(vp01, vw01);
    __m256i vmadd23 = _mm256_maddubs_epi16(vp23, vw23);
    __m256i sum = _mm256_add_epi16(vmadd01, vmadd23);
    sum = _mm256_add_epi16(sum, v32s);
    sum = _mm256_srai_epi16(sum, 6);

    __m128i lo = _mm256_castsi256_si128(sum);
    __m128i hi = _mm256_extracti128_si256(sum, 1);
    __m128i packed = _mm_packus_epi16(lo, hi);

    _mm_store_si128((__m128i*)dst, packed);

    dst += 16;
  }
}

static void angular_pred_w16_hor_high_angle_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t* delta_int, const int16_t* delta_fract, const int width, const int height, const int8_t(*filter)[4])
{
  __m256i vw0[4];
  __m256i vw1[4];
  for (int x = 0, i = 0; x < width; x += 16, ++i) {
    __m128i tmp0 = _mm_loadu_si128((__m128i*) &delta_fract[x + 0]);
    __m128i tmp1 = _mm_loadu_si128((__m128i*) &delta_fract[x + 8]);

    __m256i vidx0 = _mm256_cvtepi16_epi32(tmp0);
    __m256i vidx1 = _mm256_cvtepi16_epi32(tmp1);

    vw0[i] = _mm256_i32gather_epi32((const int32_t*)filter, vidx0, 4);
    vw1[i] = _mm256_i32gather_epi32((const int32_t*)filter, vidx1, 4);
  }
  
  for (int x = 0, vi = 0; x < width; x += 16, ++vi) {
    __m128i tmp0 = _mm_loadu_si128((__m128i*)&delta_int[x]);
    __m128i tmp1 = _mm_loadu_si128((__m128i*)&delta_int[x + 8]);
    __m256i vidx0 = _mm256_cvtepi16_epi32(tmp0);
    __m256i vidx1 = _mm256_cvtepi16_epi32(tmp1);

    // Width 16, handle one row at a time
    for (int y = 0; y < height; ++y) {
      // Do 4-tap intra interpolation filtering
      __m256i vp0 = _mm256_i32gather_epi32((const int*)&ref_main[y], vidx0, 1);
      __m256i vp1 = _mm256_i32gather_epi32((const int*)&ref_main[y], vidx1, 1);

      __m256i vmadd0 = _mm256_maddubs_epi16(vp0, vw0[vi]);
      __m256i vmadd1 = _mm256_maddubs_epi16(vp1, vw1[vi]);
      __m256i sum = _mm256_hadd_epi16(vmadd0, vmadd1);
      sum = _mm256_add_epi16(sum, _mm256_set1_epi16(32));
      sum = _mm256_srai_epi16(sum, 6);

      __m128i lo = _mm256_castsi256_si128(sum);
      __m128i hi = _mm256_extracti128_si256(sum, 1);
      __m128i packed = _mm_packus_epi16(lo, hi);
      packed = _mm_shuffle_epi32(packed, _MM_SHUFFLE(3, 1, 2, 0));

      _mm_store_si128((__m128i*)(dst + (y * width + x)), packed);
    }
  }
}

static void angular_pred_w16_hor_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t pred_mode, const int16_t multi_ref_line, const int16_t* delta_int, const int16_t* delta_fract, const int height, const int8_t(*filter)[4])
{
  const int width = 16;
  const int ref_offset = MIN(delta_int[0], delta_int[15]);
  const __m256i v32s = _mm256_set1_epi16(32);
  
  __m128i tmp0 = _mm_loadu_si128((__m128i*) &delta_fract[0]);
  __m128i tmp1 = _mm_loadu_si128((__m128i*) &delta_fract[8]);

  __m256i vidx0 = _mm256_cvtepi16_epi32(tmp0);
  __m256i vidx1 = _mm256_cvtepi16_epi32(tmp1);

  __m256i vw0 = _mm256_i32gather_epi32((const int32_t*)filter, vidx0, 4);
  __m256i vw1 = _mm256_i32gather_epi32((const int32_t*)filter, vidx1, 4);

  // Unused modes are pruned from the table and it starts from mode 5. Offset mode 5 to zero index.
  const int mode_idx = pred_mode - 5;
  const int table_offset = mode_idx * 768 + multi_ref_line * 256; // mode_idx * (3 * 256) + mrl * 256

  const __m256i vpshuf0 = _mm256_loadu_si256((__m256i*) &intra_luma_interpolation_shuffle_vectors_w64_hor[table_offset + 0]);
  const __m256i vpshuf1 = _mm256_loadu_si256((__m256i*) &intra_luma_interpolation_shuffle_vectors_w64_hor[table_offset + 32]);

  // Width 16, handle one row at a time
  for (int y = 0; y < height; ++y) {
    // Do 4-tap intra interpolation filtering
    __m128i vp = _mm_loadu_si128((__m128i*)&ref_main[y + ref_offset]);
    __m256i vp256 = _mm256_inserti128_si256(_mm256_castsi128_si256(vp), vp, 1);

    __m256i vp0 = _mm256_shuffle_epi8(vp256, vpshuf0);
    __m256i vp1 = _mm256_shuffle_epi8(vp256, vpshuf1);

    __m256i vmadd0 = _mm256_maddubs_epi16(vp0, vw0);
    __m256i vmadd1 = _mm256_maddubs_epi16(vp1, vw1);
    __m256i sum = _mm256_hadd_epi16(vmadd0, vmadd1);
    sum = _mm256_add_epi16(sum, v32s);
    sum = _mm256_srai_epi16(sum, 6);

    __m128i lo = _mm256_castsi256_si128(sum);
    __m128i hi = _mm256_extracti128_si256(sum, 1);
    __m128i packed = _mm_packus_epi16(lo, hi);
    packed = _mm_shuffle_epi32(packed, _MM_SHUFFLE(3, 1, 2, 0));

    _mm_store_si128((__m128i*)dst, packed);
    dst += 16;
  }
}

// Note: use this same function also for w64. w16 could use this, but it was slightly faster without the for loop overheads
static void angular_pred_w32_hor_avx2(uvg_pixel* dst, const uvg_pixel* ref_main, const int16_t pred_mode, const int16_t multi_ref_line, const int16_t* delta_int, const int16_t* delta_fract, const int width, const int height, const int8_t(*filter)[4])
{
  const __m256i v32s = _mm256_set1_epi16(32);

  // Unused modes are pruned from the table and it starts from mode 5. Offset mode 5 to zero index.
  const int mode_idx = pred_mode - 5;
  const int table_offset = mode_idx * 768 + multi_ref_line * 256; // mode_idx * (3 * 256) + mrl * 256

  for (int x = 0, shuf = table_offset; x < width; x += 16, shuf += 64) {
    const int ref_offset = MIN(delta_int[x], delta_int[x + 15]);

    __m128i tmp0 = _mm_loadu_si128((__m128i*)&delta_fract[x]);
    __m128i tmp1 = _mm_loadu_si128((__m128i*)&delta_fract[x + 8]);

    __m256i vidx0 = _mm256_cvtepi16_epi32(tmp0);
    __m256i vidx1 = _mm256_cvtepi16_epi32(tmp1);

    __m256i vw0 = _mm256_i32gather_epi32((const int32_t*)filter, vidx0, 4);
    __m256i vw1 = _mm256_i32gather_epi32((const int32_t*)filter, vidx1, 4);

    __m256i vpshuf0 = _mm256_loadu_si256((__m256i*) &intra_luma_interpolation_shuffle_vectors_w64_hor[shuf + 0]);
    __m256i vpshuf1 = _mm256_loadu_si256((__m256i*) &intra_luma_interpolation_shuffle_vectors_w64_hor[shuf + 32]);

    // Width 16, handle one row at a time
    for (int y = 0; y < height; ++y) {
      // Do 4-tap intra interpolation filtering
      __m128i vp = _mm_loadu_si128((__m128i*) &ref_main[y + ref_offset]);
      __m256i vp256 = _mm256_inserti128_si256(_mm256_castsi128_si256(vp), vp, 1);

      __m256i vp0 = _mm256_shuffle_epi8(vp256, vpshuf0);
      __m256i vp1 = _mm256_shuffle_epi8(vp256, vpshuf1);

      __m256i vmadd0 = _mm256_maddubs_epi16(vp0, vw0);
      __m256i vmadd1 = _mm256_maddubs_epi16(vp1, vw1);
      __m256i sum = _mm256_hadd_epi16(vmadd0, vmadd1);
      sum = _mm256_add_epi16(sum, v32s);
      sum = _mm256_srai_epi16(sum, 6);

      __m128i lo = _mm256_castsi256_si128(sum);
      __m128i hi = _mm256_extracti128_si256(sum, 1);
      __m128i packed = _mm_packus_epi16(lo, hi);
      packed = _mm_shuffle_epi32(packed, _MM_SHUFFLE(3, 1, 2, 0));

      _mm_store_si128((__m128i*)(dst + (y * width + x)), packed);
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

// Used for angles which do not require interpolation.
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

// Horizontal pixel copy for prediction mode 2.
static void angular_pred_non_fractional_angle_pxl_copy_w4_mode2_hor_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int multi_ref_offset)
{
  // const int width = 4;
  
  const __m128i vrefshuf0 = _mm_setr_epi8(
    0x00, 0x01, 0x02, 0x03, 0x01, 0x02, 0x03, 0x04,
    0x02, 0x03, 0x04, 0x05, 0x03, 0x04, 0x05, 0x06
  );

  const __m128i vrefshuf1 = _mm_setr_epi8(
    0x04, 0x05, 0x06, 0x07, 0x05, 0x06, 0x07, 0x08,
    0x06, 0x07, 0x08, 0x09, 0x07, 0x08, 0x09, 0x0a
  );

  // Handle as 4x4 blocks. There is no case where height < 4.
  if (height == 4) {
    // Offset indices by one since index 0 is top left and plus one since delta_int[0] for mode 2 is 1.
    __m128i vref = _mm_loadu_si128((__m128i*)(&ref[2] + multi_ref_offset));
    vref = _mm_shuffle_epi8(vref, vrefshuf0);

    _mm_store_si128((__m128i*)dst, vref);
  }
  else {
    // Can handle 8 rows at once
    for (int y = 0; y < height; y += 8) {
      
      __m128i vref = _mm_loadu_si128((__m128i*)(ref + 2 + y + multi_ref_offset));

      __m128i vres0 = _mm_shuffle_epi8(vref, vrefshuf0);
      __m128i vres1 = _mm_shuffle_epi8(vref, vrefshuf1);

      _mm_store_si128((__m128i*)(dst + 0), vres0);
      _mm_store_si128((__m128i*)(dst + 16), vres1);
      dst += 32;
    }
  }
}

static void angular_pred_non_fractional_angle_pxl_copy_w8_mode2_hor_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int multi_ref_offset)
{

  const __m128i vrefshuf0 = _mm_setr_epi8(
    0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
    0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08
  );

  const __m128i vrefshuf1 = _mm_setr_epi8(
    0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09,
    0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a
  );

  const __m128i vrefshuf2 = _mm_setr_epi8(
    0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b,
    0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c
  );

  const __m128i vrefshuf3 = _mm_setr_epi8(
    0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d,
    0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e
  );

  // Can handle 8 rows at once. There is no case for height 2 and 4, this function is not reached in those cases.
  for (int y = 0; y < height; y += 8) {
    // Offset indices by one since index 0 is top left and plus one since delta_int[0] for mode 2 is 1.
    __m128i vref = _mm_loadu_si128((__m128i*)(ref + 2 + y + multi_ref_offset));
    _mm_store_si128((__m128i*)(dst + 0), _mm_shuffle_epi8(vref, vrefshuf0));
    _mm_store_si128((__m128i*)(dst + 16), _mm_shuffle_epi8(vref, vrefshuf1));
    _mm_store_si128((__m128i*)(dst + 32), _mm_shuffle_epi8(vref, vrefshuf2));
    _mm_store_si128((__m128i*)(dst + 48), _mm_shuffle_epi8(vref, vrefshuf3));
    dst += 64;
  }
}

// Horizontal pixel copy for wide angles modes.
static void angular_pred_non_fractional_angle_pxl_copy_w4_wide_angle_hor_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int16_t* delta_int)
{
  // const int width = 4;

  const __m128i vtranspose = _mm_setr_epi8(
    0x00, 0x04, 0x08, 0x0c, 0x01, 0x05, 0x09, 0x0d,
    0x02, 0x06, 0x0a, 0x0e, 0x03, 0x07, 0x0b, 0x0f
  );

  //__m128i vidx = _mm_setr_epi32(delta_int[0], delta_int[1], delta_int[2], delta_int[3]);
  __m128i vidx = _mm_loadu_si128((__m128i*)delta_int);
  vidx = _mm_cvtepi16_epi32(vidx);

  // Handle as 4x4 blocks. There is no case where height < 4.
  for (int y = 0; y < height; y += 4) {
    // Offset indices by one since index 0 is top left.

    __m128i vref = _mm_i32gather_epi32((const int*)(ref + y + 1), vidx, 1);

    vref = _mm_shuffle_epi8(vref, vtranspose);

    _mm_store_si128((__m128i*)dst, vref);
    dst += 16;
  }
}

static void angular_pred_non_fractional_angle_pxl_copy_w8_wide_angle_hor_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int16_t* delta_int)
{
  // const int width = 8;

  // Place the next 4 16-bit delta int values in the lower half of the register.
  const __m128i vidxshuf = _mm_setr_epi8(
    0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff  // Don't care.
  );

  // 1st step of the transpose
  const __m256i vtranspose0 = _mm256_setr_epi8(
    0x00, 0x08, 0x01, 0x09, 0x02, 0x0a, 0x03, 0x0b,
    0x04, 0x0c, 0x05, 0x0d, 0x06, 0x0e, 0x07, 0x0f,
    0x00, 0x08, 0x01, 0x09, 0x02, 0x0a, 0x03, 0x0b,
    0x04, 0x0c, 0x05, 0x0d, 0x06, 0x0e, 0x07, 0x0f
  );

  // 3rd step of the transpose, after permute4x64_epi64
  const __m256i vtranspose1 = _mm256_setr_epi8(
    0x00, 0x01, 0x08, 0x09, 0x02, 0x03, 0x0a, 0x0b,
    0x04, 0x05, 0x0c, 0x0d, 0x06, 0x07, 0x0e, 0x0f,
    0x00, 0x01, 0x08, 0x09, 0x02, 0x03, 0x0a, 0x0b,
    0x04, 0x05, 0x0c, 0x0d, 0x06, 0x07, 0x0e, 0x0f
  );

  const __m128i vidx = _mm_loadu_si128((__m128i*)delta_int);
  const __m256i vidx0 = _mm256_cvtepi16_epi64(vidx);
  const __m256i vidx1 = _mm256_cvtepi16_epi64(_mm_shuffle_epi8(vidx, vidxshuf));

  // Handle as 8x8 blocks.
  for (int y = 0; y < height; y += 8) {
    __m256i vref0 = _mm256_i64gather_epi64((const long long*)&ref[y + 1], vidx0, 1);
    __m256i vref1 = _mm256_i64gather_epi64((const long long*)&ref[y + 1], vidx1, 1);

    // Transpose the 8x8 block
    vref0 = _mm256_shuffle_epi8(vref0, vtranspose0);
    vref1 = _mm256_shuffle_epi8(vref1, vtranspose0);
    
    vref0 = _mm256_permute4x64_epi64(vref0, _MM_SHUFFLE(3, 1, 2, 0));
    vref1 = _mm256_permute4x64_epi64(vref1, _MM_SHUFFLE(3, 1, 2, 0));

    vref0 = _mm256_shuffle_epi8(vref0, vtranspose1);
    vref1 = _mm256_shuffle_epi8(vref1, vtranspose1);

    __m256i vlo32 = _mm256_unpacklo_epi32(vref0, vref1);
    __m256i vhi32 = _mm256_unpackhi_epi32(vref0, vref1);

    __m256i vfinal0 = _mm256_permute2x128_si256(vlo32, vhi32, 0x20);
    __m256i vfinal1 = _mm256_permute2x128_si256(vlo32, vhi32, 0x31);

    _mm256_store_si256((__m256i*)(dst + 0),  vfinal0);
    _mm256_store_si256((__m256i*)(dst + 32), vfinal1);

    dst += 64;
  }
}

static void angular_pred_non_fractional_angle_pxl_copy_w16_wide_angle_hor_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int16_t* delta_int)
{
  // const int width = 16;
  
  // Handle as 16x16 blocks. This function can handle widths from 16 onwards.
  for (int y = 0; y < height; y += 16) {
    // Offset indices by one since ref[0] is top left.
    __m128i vref0 = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x00]));
    __m128i vref1 = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x01]));
    __m128i vref2 = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x02]));
    __m128i vref3 = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x03]));
    __m128i vref4 = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x04]));
    __m128i vref5 = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x05]));
    __m128i vref6 = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x06]));
    __m128i vref7 = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x07]));

    __m128i vref8 = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x08]));
    __m128i vref9 = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x09]));
    __m128i vrefa = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x0a]));
    __m128i vrefb = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x0b]));
    __m128i vrefc = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x0c]));
    __m128i vrefd = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x0d]));
    __m128i vrefe = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x0e]));
    __m128i vreff = _mm_loadu_si128((__m128i*)(ref + y + 1 + delta_int[0x0f]));

    // The result is just a transpose of the 16x16 block.
    __m128i vlo8_0 = _mm_unpacklo_epi8(vref0, vref1);
    __m128i vlo8_1 = _mm_unpacklo_epi8(vref2, vref3);
    __m128i vlo8_2 = _mm_unpacklo_epi8(vref4, vref5);
    __m128i vlo8_3 = _mm_unpacklo_epi8(vref6, vref7);
    __m128i vlo8_4 = _mm_unpacklo_epi8(vref8, vref9);
    __m128i vlo8_5 = _mm_unpacklo_epi8(vrefa, vrefb);
    __m128i vlo8_6 = _mm_unpacklo_epi8(vrefc, vrefd);
    __m128i vlo8_7 = _mm_unpacklo_epi8(vrefe, vreff);

    __m128i vhi8_0 = _mm_unpackhi_epi8(vref0, vref1);
    __m128i vhi8_1 = _mm_unpackhi_epi8(vref2, vref3);
    __m128i vhi8_2 = _mm_unpackhi_epi8(vref4, vref5);
    __m128i vhi8_3 = _mm_unpackhi_epi8(vref6, vref7);
    __m128i vhi8_4 = _mm_unpackhi_epi8(vref8, vref9);
    __m128i vhi8_5 = _mm_unpackhi_epi8(vrefa, vrefb);
    __m128i vhi8_6 = _mm_unpackhi_epi8(vrefc, vrefd);
    __m128i vhi8_7 = _mm_unpackhi_epi8(vrefe, vreff);

    __m128i vlo16_0 = _mm_unpacklo_epi16(vlo8_0, vlo8_1);
    __m128i vlo16_1 = _mm_unpacklo_epi16(vlo8_2, vlo8_3);
    __m128i vlo16_2 = _mm_unpacklo_epi16(vhi8_0, vhi8_1);
    __m128i vlo16_3 = _mm_unpacklo_epi16(vhi8_2, vhi8_3);
    __m128i vlo16_4 = _mm_unpacklo_epi16(vlo8_4, vlo8_5);
    __m128i vlo16_5 = _mm_unpacklo_epi16(vlo8_6, vlo8_7);
    __m128i vlo16_6 = _mm_unpacklo_epi16(vhi8_4, vhi8_5);
    __m128i vlo16_7 = _mm_unpacklo_epi16(vhi8_6, vhi8_7);


    __m128i vhi16_0 = _mm_unpackhi_epi16(vlo8_0, vlo8_1);
    __m128i vhi16_1 = _mm_unpackhi_epi16(vlo8_2, vlo8_3);
    __m128i vhi16_2 = _mm_unpackhi_epi16(vhi8_0, vhi8_1);
    __m128i vhi16_3 = _mm_unpackhi_epi16(vhi8_2, vhi8_3);
    __m128i vhi16_4 = _mm_unpackhi_epi16(vlo8_4, vlo8_5);
    __m128i vhi16_5 = _mm_unpackhi_epi16(vlo8_6, vlo8_7);
    __m128i vhi16_6 = _mm_unpackhi_epi16(vhi8_4, vhi8_5);
    __m128i vhi16_7 = _mm_unpackhi_epi16(vhi8_6, vhi8_7);

    __m128i vlo32_0 = _mm_unpacklo_epi32(vlo16_0, vlo16_1);
    __m128i vlo32_1 = _mm_unpacklo_epi32(vlo16_2, vlo16_3);
    __m128i vlo32_2 = _mm_unpacklo_epi32(vhi16_0, vhi16_1);
    __m128i vlo32_3 = _mm_unpacklo_epi32(vhi16_2, vhi16_3);
    __m128i vlo32_4 = _mm_unpacklo_epi32(vlo16_4, vlo16_5);
    __m128i vlo32_5 = _mm_unpacklo_epi32(vlo16_6, vlo16_7);
    __m128i vlo32_6 = _mm_unpacklo_epi32(vhi16_4, vhi16_5);
    __m128i vlo32_7 = _mm_unpacklo_epi32(vhi16_6, vhi16_7);

    __m128i vhi32_0 = _mm_unpackhi_epi32(vlo16_0, vlo16_1);
    __m128i vhi32_1 = _mm_unpackhi_epi32(vlo16_2, vlo16_3);
    __m128i vhi32_2 = _mm_unpackhi_epi32(vhi16_0, vhi16_1);
    __m128i vhi32_3 = _mm_unpackhi_epi32(vhi16_2, vhi16_3);
    __m128i vhi32_4 = _mm_unpackhi_epi32(vlo16_4, vlo16_5);
    __m128i vhi32_5 = _mm_unpackhi_epi32(vlo16_6, vlo16_7);
    __m128i vhi32_6 = _mm_unpackhi_epi32(vhi16_4, vhi16_5);
    __m128i vhi32_7 = _mm_unpackhi_epi32(vhi16_6, vhi16_7);

    __m128i vrow0 = _mm_unpacklo_epi64(vlo32_0, vlo32_4);
    __m128i vrow1 = _mm_unpackhi_epi64(vlo32_0, vlo32_4);
    __m128i vrow2 = _mm_unpacklo_epi64(vhi32_0, vhi32_4);
    __m128i vrow3 = _mm_unpackhi_epi64(vhi32_0, vhi32_4);
    __m128i vrow4 = _mm_unpacklo_epi64(vlo32_2, vlo32_6);
    __m128i vrow5 = _mm_unpackhi_epi64(vlo32_2, vlo32_6);
    __m128i vrow6 = _mm_unpacklo_epi64(vhi32_2, vhi32_6);
    __m128i vrow7 = _mm_unpackhi_epi64(vhi32_2, vhi32_6);

    __m128i vrow8 = _mm_unpacklo_epi64(vlo32_1, vlo32_5);
    __m128i vrwo9 = _mm_unpackhi_epi64(vlo32_1, vlo32_5);
    __m128i vrowa = _mm_unpacklo_epi64(vhi32_1, vhi32_5);
    __m128i vrowb = _mm_unpackhi_epi64(vhi32_1, vhi32_5);
    __m128i vrowc = _mm_unpacklo_epi64(vlo32_3, vlo32_7);
    __m128i vrowd = _mm_unpackhi_epi64(vlo32_3, vlo32_7);
    __m128i vrowe = _mm_unpacklo_epi64(vhi32_3, vhi32_7);
    __m128i vrowf = _mm_unpackhi_epi64(vhi32_3, vhi32_7);

    _mm_store_si128((__m128i*)(dst + 0), vrow0);
    _mm_store_si128((__m128i*)(dst + 16), vrow1);
    _mm_store_si128((__m128i*)(dst + 32), vrow2);
    _mm_store_si128((__m128i*)(dst + 48), vrow3);
    _mm_store_si128((__m128i*)(dst + 64), vrow4);
    _mm_store_si128((__m128i*)(dst + 80), vrow5);
    _mm_store_si128((__m128i*)(dst + 96), vrow6);
    _mm_store_si128((__m128i*)(dst + 112), vrow7);

    _mm_store_si128((__m128i*)(dst + 128), vrow8);
    _mm_store_si128((__m128i*)(dst + 144), vrwo9);
    _mm_store_si128((__m128i*)(dst + 160), vrowa);
    _mm_store_si128((__m128i*)(dst + 176), vrowb);
    _mm_store_si128((__m128i*)(dst + 192), vrowc);
    _mm_store_si128((__m128i*)(dst + 208), vrowd);
    _mm_store_si128((__m128i*)(dst + 224), vrowe);
    _mm_store_si128((__m128i*)(dst + 240), vrowf);

    dst += 256;
  }
}

static void angular_pred_non_fractional_angle_pxl_copy_w32_wide_angle_hor_avx2(uvg_pixel* dst, uvg_pixel* ref, const int height, const int16_t* delta_int)
{
  // const int width = 32;
  // Handle as 32x32 blocks. Similarly to the w16 version, this is also just a transpose of the 32x32 block.
  // TODO: if this is too slow, consider doing it in 16x16 blocks. There will be a lot of moving data between registers in this solution.
  for (int y = 0; y < height; y += 32) {
    // Offset indices by one since ref[0] is top left.
    __m256i vref00 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x00]));
    __m256i vref01 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x01]));
    __m256i vref02 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x02]));
    __m256i vref03 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x03]));
    __m256i vref04 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x04]));
    __m256i vref05 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x05]));
    __m256i vref06 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x06]));
    __m256i vref07 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x07]));

    __m256i vref08 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x08]));
    __m256i vref09 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x09]));
    __m256i vref0a = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x0a]));
    __m256i vref0b = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x0b]));
    __m256i vref0c = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x0c]));
    __m256i vref0d = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x0d]));
    __m256i vref0e = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x0e]));
    __m256i vref0f = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x0f]));

    __m256i vref10 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x10]));
    __m256i vref11 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x11]));
    __m256i vref12 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x12]));
    __m256i vref13 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x13]));
    __m256i vref14 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x14]));
    __m256i vref15 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x15]));
    __m256i vref16 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x16]));
    __m256i vref17 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x17]));

    __m256i vref18 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x18]));
    __m256i vref19 = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x19]));
    __m256i vref1a = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x1a]));
    __m256i vref1b = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x1b]));
    __m256i vref1c = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x1c]));
    __m256i vref1d = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x1d]));
    __m256i vref1e = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x1e]));
    __m256i vref1f = _mm256_loadu_si256((__m256i*)(ref + y + 1 + delta_int[0x1f]));

    __m256i vlo8_0 = _mm256_unpacklo_epi8(vref00, vref01);
    __m256i vlo8_1 = _mm256_unpacklo_epi8(vref02, vref03);
    __m256i vlo8_2 = _mm256_unpacklo_epi8(vref04, vref05);
    __m256i vlo8_3 = _mm256_unpacklo_epi8(vref06, vref07);
    __m256i vlo8_4 = _mm256_unpacklo_epi8(vref08, vref09);
    __m256i vlo8_5 = _mm256_unpacklo_epi8(vref0a, vref0b);
    __m256i vlo8_6 = _mm256_unpacklo_epi8(vref0c, vref0d);
    __m256i vlo8_7 = _mm256_unpacklo_epi8(vref0e, vref0f);
    __m256i vlo8_8 = _mm256_unpacklo_epi8(vref10, vref11);
    __m256i vlo8_9 = _mm256_unpacklo_epi8(vref12, vref13);
    __m256i vlo8_a = _mm256_unpacklo_epi8(vref14, vref15);
    __m256i vlo8_b = _mm256_unpacklo_epi8(vref16, vref17);
    __m256i vlo8_c = _mm256_unpacklo_epi8(vref18, vref19);
    __m256i vlo8_d = _mm256_unpacklo_epi8(vref1a, vref1b);
    __m256i vlo8_e = _mm256_unpacklo_epi8(vref1c, vref1d);
    __m256i vlo8_f = _mm256_unpacklo_epi8(vref1e, vref1f);

    __m256i vhi8_0 = _mm256_unpackhi_epi8(vref00, vref01);
    __m256i vhi8_1 = _mm256_unpackhi_epi8(vref02, vref03);
    __m256i vhi8_2 = _mm256_unpackhi_epi8(vref04, vref05);
    __m256i vhi8_3 = _mm256_unpackhi_epi8(vref06, vref07);
    __m256i vhi8_4 = _mm256_unpackhi_epi8(vref08, vref09);
    __m256i vhi8_5 = _mm256_unpackhi_epi8(vref0a, vref0b);
    __m256i vhi8_6 = _mm256_unpackhi_epi8(vref0c, vref0d);
    __m256i vhi8_7 = _mm256_unpackhi_epi8(vref0e, vref0f);
    __m256i vhi8_8 = _mm256_unpackhi_epi8(vref10, vref11);
    __m256i vhi8_9 = _mm256_unpackhi_epi8(vref12, vref13);
    __m256i vhi8_a = _mm256_unpackhi_epi8(vref14, vref15);
    __m256i vhi8_b = _mm256_unpackhi_epi8(vref16, vref17);
    __m256i vhi8_c = _mm256_unpackhi_epi8(vref18, vref19);
    __m256i vhi8_d = _mm256_unpackhi_epi8(vref1a, vref1b);
    __m256i vhi8_e = _mm256_unpackhi_epi8(vref1c, vref1d);
    __m256i vhi8_f = _mm256_unpackhi_epi8(vref1e, vref1f);

    __m256i vlo16_0 = _mm256_unpacklo_epi16(vlo8_0, vlo8_1);
    __m256i vlo16_1 = _mm256_unpacklo_epi16(vlo8_2, vlo8_3);
    __m256i vlo16_2 = _mm256_unpacklo_epi16(vlo8_4, vlo8_5);
    __m256i vlo16_3 = _mm256_unpacklo_epi16(vlo8_6, vlo8_7);
    __m256i vlo16_4 = _mm256_unpacklo_epi16(vlo8_8, vlo8_9);
    __m256i vlo16_5 = _mm256_unpacklo_epi16(vlo8_a, vlo8_b);
    __m256i vlo16_6 = _mm256_unpacklo_epi16(vlo8_c, vlo8_d);
    __m256i vlo16_7 = _mm256_unpacklo_epi16(vlo8_e, vlo8_f);
    __m256i vlo16_8 = _mm256_unpacklo_epi16(vhi8_0, vhi8_1);
    __m256i vlo16_9 = _mm256_unpacklo_epi16(vhi8_2, vhi8_3);
    __m256i vlo16_a = _mm256_unpacklo_epi16(vhi8_4, vhi8_5);
    __m256i vlo16_b = _mm256_unpacklo_epi16(vhi8_6, vhi8_7);
    __m256i vlo16_c = _mm256_unpacklo_epi16(vhi8_8, vhi8_9);
    __m256i vlo16_d = _mm256_unpacklo_epi16(vhi8_a, vhi8_b);
    __m256i vlo16_e = _mm256_unpacklo_epi16(vhi8_c, vhi8_d);
    __m256i vlo16_f = _mm256_unpacklo_epi16(vhi8_e, vhi8_f);

    __m256i vhi16_0 = _mm256_unpackhi_epi16(vlo8_0, vlo8_1);
    __m256i vhi16_1 = _mm256_unpackhi_epi16(vlo8_2, vlo8_3);
    __m256i vhi16_2 = _mm256_unpackhi_epi16(vlo8_4, vlo8_5);
    __m256i vhi16_3 = _mm256_unpackhi_epi16(vlo8_6, vlo8_7);
    __m256i vhi16_4 = _mm256_unpackhi_epi16(vlo8_8, vlo8_9);
    __m256i vhi16_5 = _mm256_unpackhi_epi16(vlo8_a, vlo8_b);
    __m256i vhi16_6 = _mm256_unpackhi_epi16(vlo8_c, vlo8_d);
    __m256i vhi16_7 = _mm256_unpackhi_epi16(vlo8_e, vlo8_f);
    __m256i vhi16_8 = _mm256_unpackhi_epi16(vhi8_0, vhi8_1);
    __m256i vhi16_9 = _mm256_unpackhi_epi16(vhi8_2, vhi8_3);
    __m256i vhi16_a = _mm256_unpackhi_epi16(vhi8_4, vhi8_5);
    __m256i vhi16_b = _mm256_unpackhi_epi16(vhi8_6, vhi8_7);
    __m256i vhi16_c = _mm256_unpackhi_epi16(vhi8_8, vhi8_9);
    __m256i vhi16_d = _mm256_unpackhi_epi16(vhi8_a, vhi8_b);
    __m256i vhi16_e = _mm256_unpackhi_epi16(vhi8_c, vhi8_d);
    __m256i vhi16_f = _mm256_unpackhi_epi16(vhi8_e, vhi8_f);

    __m256i vlo32_0 = _mm256_unpacklo_epi32(vlo16_0, vlo16_1);
    __m256i vlo32_1 = _mm256_unpacklo_epi32(vlo16_2, vlo16_3);
    __m256i vlo32_2 = _mm256_unpacklo_epi32(vlo16_4, vlo16_5);
    __m256i vlo32_3 = _mm256_unpacklo_epi32(vlo16_6, vlo16_7);
    __m256i vlo32_4 = _mm256_unpacklo_epi32(vhi16_0, vhi16_1);
    __m256i vlo32_5 = _mm256_unpacklo_epi32(vhi16_2, vhi16_3);
    __m256i vlo32_6 = _mm256_unpacklo_epi32(vhi16_4, vhi16_5);
    __m256i vlo32_7 = _mm256_unpacklo_epi32(vhi16_6, vhi16_7);
    __m256i vlo32_8 = _mm256_unpacklo_epi32(vlo16_8, vlo16_9);
    __m256i vlo32_9 = _mm256_unpacklo_epi32(vlo16_a, vlo16_b);
    __m256i vlo32_a = _mm256_unpacklo_epi32(vlo16_c, vlo16_d);
    __m256i vlo32_b = _mm256_unpacklo_epi32(vlo16_e, vlo16_f);
    __m256i vlo32_c = _mm256_unpacklo_epi32(vhi16_8, vhi16_9);
    __m256i vlo32_d = _mm256_unpacklo_epi32(vhi16_a, vhi16_b);
    __m256i vlo32_e = _mm256_unpacklo_epi32(vhi16_c, vhi16_d);
    __m256i vlo32_f = _mm256_unpacklo_epi32(vhi16_e, vhi16_f);

    __m256i vhi32_0 = _mm256_unpackhi_epi32(vlo16_0, vlo16_1);
    __m256i vhi32_1 = _mm256_unpackhi_epi32(vlo16_2, vlo16_3);
    __m256i vhi32_2 = _mm256_unpackhi_epi32(vlo16_4, vlo16_5);
    __m256i vhi32_3 = _mm256_unpackhi_epi32(vlo16_6, vlo16_7);
    __m256i vhi32_4 = _mm256_unpackhi_epi32(vhi16_0, vhi16_1);
    __m256i vhi32_5 = _mm256_unpackhi_epi32(vhi16_2, vhi16_3);
    __m256i vhi32_6 = _mm256_unpackhi_epi32(vhi16_4, vhi16_5);
    __m256i vhi32_7 = _mm256_unpackhi_epi32(vhi16_6, vhi16_7);
    __m256i vhi32_8 = _mm256_unpackhi_epi32(vlo16_8, vlo16_9);
    __m256i vhi32_9 = _mm256_unpackhi_epi32(vlo16_a, vlo16_b);
    __m256i vhi32_a = _mm256_unpackhi_epi32(vlo16_c, vlo16_d);
    __m256i vhi32_b = _mm256_unpackhi_epi32(vlo16_e, vlo16_f);
    __m256i vhi32_c = _mm256_unpackhi_epi32(vhi16_8, vhi16_9);
    __m256i vhi32_d = _mm256_unpackhi_epi32(vhi16_a, vhi16_b);
    __m256i vhi32_e = _mm256_unpackhi_epi32(vhi16_c, vhi16_d);
    __m256i vhi32_f = _mm256_unpackhi_epi32(vhi16_e, vhi16_f);

    __m256i vlo64_0 = _mm256_unpacklo_epi64(vlo32_0, vlo32_1);
    __m256i vlo64_1 = _mm256_unpacklo_epi64(vlo32_2, vlo32_3);
    __m256i vlo64_2 = _mm256_unpacklo_epi64(vhi32_0, vhi32_1);
    __m256i vlo64_3 = _mm256_unpacklo_epi64(vhi32_2, vhi32_3);
    __m256i vlo64_4 = _mm256_unpacklo_epi64(vlo32_4, vlo32_5);
    __m256i vlo64_5 = _mm256_unpacklo_epi64(vlo32_6, vlo32_7);
    __m256i vlo64_6 = _mm256_unpacklo_epi64(vhi32_4, vhi32_5);
    __m256i vlo64_7 = _mm256_unpacklo_epi64(vhi32_6, vhi32_7);
    __m256i vlo64_8 = _mm256_unpacklo_epi64(vlo32_8, vlo32_9);
    __m256i vlo64_9 = _mm256_unpacklo_epi64(vlo32_a, vlo32_b);
    __m256i vlo64_a = _mm256_unpacklo_epi64(vhi32_8, vhi32_9);
    __m256i vlo64_b = _mm256_unpacklo_epi64(vhi32_a, vhi32_b);
    __m256i vlo64_c = _mm256_unpacklo_epi64(vlo32_c, vlo32_d);
    __m256i vlo64_d = _mm256_unpacklo_epi64(vlo32_e, vlo32_f);
    __m256i vlo64_e = _mm256_unpacklo_epi64(vhi32_c, vhi32_d);
    __m256i vlo64_f = _mm256_unpacklo_epi64(vhi32_e, vhi32_f);

    __m256i vhi64_0 = _mm256_unpackhi_epi64(vlo32_0, vlo32_1);
    __m256i vhi64_1 = _mm256_unpackhi_epi64(vlo32_2, vlo32_3);
    __m256i vhi64_2 = _mm256_unpackhi_epi64(vhi32_0, vhi32_1);
    __m256i vhi64_3 = _mm256_unpackhi_epi64(vhi32_2, vhi32_3);
    __m256i vhi64_4 = _mm256_unpackhi_epi64(vlo32_4, vlo32_5);
    __m256i vhi64_5 = _mm256_unpackhi_epi64(vlo32_6, vlo32_7);
    __m256i vhi64_6 = _mm256_unpackhi_epi64(vhi32_4, vhi32_5);
    __m256i vhi64_7 = _mm256_unpackhi_epi64(vhi32_6, vhi32_7);
    __m256i vhi64_8 = _mm256_unpackhi_epi64(vlo32_8, vlo32_9);
    __m256i vhi64_9 = _mm256_unpackhi_epi64(vlo32_a, vlo32_b);
    __m256i vhi64_a = _mm256_unpackhi_epi64(vhi32_8, vhi32_9);
    __m256i vhi64_b = _mm256_unpackhi_epi64(vhi32_a, vhi32_b);
    __m256i vhi64_c = _mm256_unpackhi_epi64(vlo32_c, vlo32_d);
    __m256i vhi64_d = _mm256_unpackhi_epi64(vlo32_e, vlo32_f);
    __m256i vhi64_e = _mm256_unpackhi_epi64(vhi32_c, vhi32_d);
    __m256i vhi64_f = _mm256_unpackhi_epi64(vhi32_e, vhi32_f);

    __m256i vrow00 = _mm256_permute2x128_si256(vlo64_0, vlo64_1, 0x20);
    __m256i vrow01 = _mm256_permute2x128_si256(vhi64_0, vhi64_1, 0x20);
    __m256i vrow02 = _mm256_permute2x128_si256(vlo64_2, vlo64_3, 0x20);
    __m256i vrow03 = _mm256_permute2x128_si256(vhi64_2, vhi64_3, 0x20);
    __m256i vrow04 = _mm256_permute2x128_si256(vlo64_4, vlo64_5, 0x20);
    __m256i vrow05 = _mm256_permute2x128_si256(vhi64_4, vhi64_5, 0x20);
    __m256i vrow06 = _mm256_permute2x128_si256(vlo64_6, vlo64_7, 0x20);
    __m256i vrow07 = _mm256_permute2x128_si256(vhi64_6, vhi64_7, 0x20);

    __m256i vrow08 = _mm256_permute2x128_si256(vlo64_8, vlo64_9, 0x20);
    __m256i vrow09 = _mm256_permute2x128_si256(vhi64_8, vhi64_9, 0x20);
    __m256i vrow0a = _mm256_permute2x128_si256(vlo64_a, vlo64_b, 0x20);
    __m256i vrow0b = _mm256_permute2x128_si256(vhi64_a, vhi64_b, 0x20);
    __m256i vrow0c = _mm256_permute2x128_si256(vlo64_c, vlo64_d, 0x20);
    __m256i vrow0d = _mm256_permute2x128_si256(vhi64_c, vhi64_d, 0x20);
    __m256i vrow0e = _mm256_permute2x128_si256(vlo64_e, vlo64_f, 0x20);
    __m256i vrow0f = _mm256_permute2x128_si256(vhi64_e, vhi64_f, 0x20);

    __m256i vrow10 = _mm256_permute2x128_si256(vlo64_0, vlo64_1, 0x31);
    __m256i vrow11 = _mm256_permute2x128_si256(vhi64_0, vhi64_1, 0x31);
    __m256i vrow12 = _mm256_permute2x128_si256(vlo64_2, vlo64_3, 0x31);
    __m256i vrow13 = _mm256_permute2x128_si256(vhi64_2, vhi64_3, 0x31);
    __m256i vrow14 = _mm256_permute2x128_si256(vlo64_4, vlo64_5, 0x31);
    __m256i vrow15 = _mm256_permute2x128_si256(vhi64_4, vhi64_5, 0x31);
    __m256i vrow16 = _mm256_permute2x128_si256(vlo64_6, vlo64_7, 0x31);
    __m256i vrow17 = _mm256_permute2x128_si256(vhi64_6, vhi64_7, 0x31);

    __m256i vrow18 = _mm256_permute2x128_si256(vlo64_8, vlo64_9, 0x31);
    __m256i vrow19 = _mm256_permute2x128_si256(vhi64_8, vhi64_9, 0x31);
    __m256i vrow1a = _mm256_permute2x128_si256(vlo64_a, vlo64_b, 0x31);
    __m256i vrow1b = _mm256_permute2x128_si256(vhi64_a, vhi64_b, 0x31);
    __m256i vrow1c = _mm256_permute2x128_si256(vlo64_c, vlo64_d, 0x31);
    __m256i vrow1d = _mm256_permute2x128_si256(vhi64_c, vhi64_d, 0x31);
    __m256i vrow1e = _mm256_permute2x128_si256(vlo64_e, vlo64_f, 0x31);
    __m256i vrow1f = _mm256_permute2x128_si256(vhi64_e, vhi64_f, 0x31);
    
    _mm256_store_si256((__m256i*)(dst + 0),   vrow00);
    _mm256_store_si256((__m256i*)(dst + 32),  vrow01);
    _mm256_store_si256((__m256i*)(dst + 64),  vrow02);
    _mm256_store_si256((__m256i*)(dst + 96),  vrow03);
    _mm256_store_si256((__m256i*)(dst + 128), vrow04);
    _mm256_store_si256((__m256i*)(dst + 160), vrow05);
    _mm256_store_si256((__m256i*)(dst + 192), vrow06);
    _mm256_store_si256((__m256i*)(dst + 224), vrow07);
    _mm256_store_si256((__m256i*)(dst + 256), vrow08);
    _mm256_store_si256((__m256i*)(dst + 288), vrow09);
    _mm256_store_si256((__m256i*)(dst + 320), vrow0a);
    _mm256_store_si256((__m256i*)(dst + 352), vrow0b);
    _mm256_store_si256((__m256i*)(dst + 384), vrow0c);
    _mm256_store_si256((__m256i*)(dst + 416), vrow0d);
    _mm256_store_si256((__m256i*)(dst + 448), vrow0e);
    _mm256_store_si256((__m256i*)(dst + 480), vrow0f);
    _mm256_store_si256((__m256i*)(dst + 512), vrow10);
    _mm256_store_si256((__m256i*)(dst + 544), vrow11);
    _mm256_store_si256((__m256i*)(dst + 576), vrow12);
    _mm256_store_si256((__m256i*)(dst + 608), vrow13);
    _mm256_store_si256((__m256i*)(dst + 640), vrow14);
    _mm256_store_si256((__m256i*)(dst + 672), vrow15);
    _mm256_store_si256((__m256i*)(dst + 704), vrow16);
    _mm256_store_si256((__m256i*)(dst + 736), vrow17);
    _mm256_store_si256((__m256i*)(dst + 768), vrow18);
    _mm256_store_si256((__m256i*)(dst + 800), vrow19);
    _mm256_store_si256((__m256i*)(dst + 832), vrow1a);
    _mm256_store_si256((__m256i*)(dst + 864), vrow1b);
    _mm256_store_si256((__m256i*)(dst + 896), vrow1c);
    _mm256_store_si256((__m256i*)(dst + 928), vrow1d);
    _mm256_store_si256((__m256i*)(dst + 960), vrow1e);
    _mm256_store_si256((__m256i*)(dst + 992), vrow1f);

    dst += 1024;
  }
}



static void angular_pdpc_ver_w8_high_angle_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int height, const int scale, const int mode_disp)
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

    __m128i vdst = _mm_i64gather_epi64((const long long int*)(dst + y * width), vseq, 8);
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


// Mode 18

static void angular_pdpc_mode18_w4_avx2(uvg_pixel* dst, const uvg_pixel top_left, const uvg_pixel* ref_side, const int height, const int scale)
{
  const int width = 4;
  const int limit = MIN(3 << scale, height);

  __m256i v32s = _mm256_set1_epi16(32);

  const uint32_t ref4 = *(uint32_t*)&ref_side[1];

  __m128i vref = _mm_set1_epi32(ref4);
  __m256i vref16 = _mm256_cvtepu8_epi16(vref);

  __m256i vtopleft = _mm256_set1_epi16((uint16_t)top_left);

  // Weight table offset
  const int table_offset = scale * 64;

  for (int y = 0, o = 0; y < limit; y += 4, o += 16) {
    const int offset = table_offset + o;

    __m128i vpred = _mm_load_si128((__m128i*)(dst + y * width));

    __m256i vpred16 = _mm256_cvtepu8_epi16(vpred);
    __m256i vweight = _mm256_load_si256((const __m256i*) & intra_pdpc_w4_hor_weight[offset]);

    __m256i accu = _mm256_sub_epi16(vref16, vtopleft);
    accu = _mm256_mullo_epi16(vweight, accu);
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

  __m256i v32s = _mm256_set1_epi16(32);

  const uint64_t ref8 = *(uint64_t*)&ref_side[1];

  __m128i vref = _mm_set1_epi64x(ref8);
  __m256i vref16 = _mm256_cvtepu8_epi16(vref);

  __m256i vtopleft = _mm256_set1_epi16((uint16_t)top_left);

  // Weight table offset
  const int table_offset = scale * 128;

  for (int y = 0, o = table_offset; y < limit; y += 2, o += 16) {
    const __m256i vwT = _mm256_load_si256((const __m256i*) & intra_pdpc_w8_hor_weight[o]);

    __m128i vpred = _mm_load_si128((__m128i*)(dst + y * width));
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

  __m128i vref = _mm_loadu_si128((const __m128i*) & ref_side[1]);
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

  __m128i vrefa = _mm_loadu_si128((const __m128i*) & ref_side[1]);
  __m256i vref16a = _mm256_cvtepu8_epi16(vrefa);

  __m128i vrefb = _mm_loadu_si128((const __m128i*) & ref_side[17]);
  __m256i vref16b = _mm256_cvtepu8_epi16(vrefb);

  __m256i vtopleft = _mm256_set1_epi16((uint16_t)top_left);

  // Handle one line at a time. Skip line if vertical limit reached.
  for (int y = 0; y < limit; ++y) {
    const int16_t wT = 32 >> (2 * (y + 0) >> scale);
    __m256i vwT = _mm256_set1_epi16(wT);

    // Calculate first half
    __m128i vpred = _mm_load_si128((__m128i*)(dst + (y * width + 0)));
    __m256i vpred16 = _mm256_cvtepu8_epi16(vpred);

    __m256i accu0 = _mm256_sub_epi16(vref16a, vtopleft);
    accu0 = _mm256_mullo_epi16(vwT, accu0);
    accu0 = _mm256_add_epi16(accu0, v32s);
    accu0 = _mm256_srai_epi16(accu0, 6);
    accu0 = _mm256_add_epi16(vpred16, accu0);

    // Calculate second half
    vpred = _mm_load_si128((__m128i*)(dst + (y * width + 16)));
    vpred16 = _mm256_cvtepu8_epi16(vpred);

    __m256i accu1 = _mm256_sub_epi16(vref16b, vtopleft);
    accu1 = _mm256_mullo_epi16(vwT, accu1);
    accu1 = _mm256_add_epi16(accu1, v32s);
    accu1 = _mm256_srai_epi16(accu1, 6);
    accu1 = _mm256_add_epi16(vpred16, accu1);

    // Store results
    __m256i packed = _mm256_packus_epi16(accu0, accu1);
    packed = _mm256_permute4x64_epi64(packed, _MM_SHUFFLE(3, 1, 2, 0));

    _mm256_store_si256((__m256i*)(dst + (y * width)), packed);
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

    __m256i accu0 = _mm256_sub_epi16(vref16a, vtopleft);
    accu0 = _mm256_mullo_epi16(vwT, accu0);
    accu0 = _mm256_add_epi16(accu0, v32s);
    accu0 = _mm256_srai_epi16(accu0, 6);
    accu0 = _mm256_add_epi16(vpred16, accu0);

    // Calculate second quarter
    vpred = _mm_load_si128((__m128i*)(dst + (y * width + 16)));
    vpred16 = _mm256_cvtepu8_epi16(vpred);

    __m256i accu1 = _mm256_sub_epi16(vref16b, vtopleft);
    accu1 = _mm256_mullo_epi16(vwT, accu1);
    accu1 = _mm256_add_epi16(accu1, v32s);
    accu1 = _mm256_srai_epi16(accu1, 6);
    accu1 = _mm256_add_epi16(vpred16, accu1);

    // Calculate third quarter
    vpred = _mm_load_si128((__m128i*)(dst + (y * width + 32)));
    vpred16 = _mm256_cvtepu8_epi16(vpred);

    __m256i accu2 = _mm256_sub_epi16(vref16c, vtopleft);
    accu2 = _mm256_mullo_epi16(vwT, accu2);
    accu2 = _mm256_add_epi16(accu2, v32s);
    accu2 = _mm256_srai_epi16(accu2, 6);
    accu2 = _mm256_add_epi16(vpred16, accu2);

    // Calculate fourth quarter
    vpred = _mm_load_si128((__m128i*)(dst + (y * width + 48)));
    vpred16 = _mm256_cvtepu8_epi16(vpred);

    __m256i accu3 = _mm256_sub_epi16(vref16d, vtopleft);
    accu3 = _mm256_mullo_epi16(vwT, accu3);
    accu3 = _mm256_add_epi16(accu3, v32s);
    accu3 = _mm256_srai_epi16(accu3, 6);
    accu3 = _mm256_add_epi16(vpred16, accu3);

    __m256i packed0 = _mm256_packus_epi16(accu0, accu1);
    __m256i packed1 = _mm256_packus_epi16(accu2, accu3);
    packed0 = _mm256_permute4x64_epi64(packed0, _MM_SHUFFLE(3, 1, 2, 0));
    packed1 = _mm256_permute4x64_epi64(packed1, _MM_SHUFFLE(3, 1, 2, 0));

    _mm256_store_si256((__m256i*)(dst + (y * width + 0)),  packed0);
    _mm256_store_si256((__m256i*)(dst + (y * width + 32)), packed1);
  }
}


// Vertical modes

static void angular_pdpc_ver_w4_high_angle_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int height, const int scale, const int mode_disp)
{
  const int width = 4;
  //ALIGNED(32) uint8_t left[4][4];
  __m128i v32s = _mm_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int offset = scale * 16;
  const __m128i vweight = _mm_load_si128((const __m128i*) &intra_pdpc_w4_ver_improved_weight[offset]);

  const int inv_angle_offset = mode_disp * 64;
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  const __m128i vleftshuf = _mm_setr_epi8(
    0x00, 0x04, 0x08, 0x0c, 0x01, 0x05, 0x09, 0x0d, 
    0x02, 0x06, 0x0a, 0x0e, 0x03, 0x07, 0x0b, 0x0f
  );

  __m128i vidx = _mm_setr_epi32(shifted_inv_angle_sum[0], shifted_inv_angle_sum[1], 
                                shifted_inv_angle_sum[2], shifted_inv_angle_sum[3]);

  // For a 4 width block, height must be at least 4. Handle 4 lines at once.
  for (int y = 0; y < height; y += 4) {
    /*for (int xx = 0; xx < width; ++xx) {
      memcpy(left[xx], &ref_side[(y + 0) + shifted_inv_angle_sum[xx] + 1], 4 * sizeof(uint8_t));
    }*/

    __m128i vdst = _mm_loadu_si128((const __m128i*)(dst + y * width));
    //__m128i vleft = _mm_load_si128((__m128i*)left);
    __m128i vleft = _mm_i32gather_epi32((const int32_t*)&ref_side[y + 1], vidx, 1);
    vleft = _mm_shuffle_epi8(vleft, vleftshuf);

    __m128i vlo = _mm_unpacklo_epi8(vdst, vleft);
    __m128i vhi = _mm_unpackhi_epi8(vdst, vleft);

    __m128i vmaddlo = _mm_maddubs_epi16(vlo, vweight);
    __m128i vmaddhi = _mm_maddubs_epi16(vhi, vweight);

    vmaddlo = _mm_add_epi16(vmaddlo, v32s);
    vmaddhi = _mm_add_epi16(vmaddhi, v32s);

    vmaddlo = _mm_srai_epi16(vmaddlo, 6);
    vmaddhi = _mm_srai_epi16(vmaddhi, 6);

    __m128i packed = _mm_packus_epi16(vmaddlo, vmaddhi);

    _mm_store_si128((__m128i*)(dst + (y * width)), packed);
  }
}

static void angular_pdpc_ver_w4_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int height, const int scale, const int mode_disp)
{
  const int width = 4;
  __m128i v32s = _mm_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int offset = scale * 16;
  const int inv_angle_offset = mode_disp * 64;
  const int shuf_offset = mode_disp * 16;

  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];
  const __m128i vweight = _mm_load_si128((const __m128i*) &intra_pdpc_w4_ver_improved_weight[offset]);
  const __m128i vshuf = _mm_loadu_si128((__m128i*) &intra_pdpc_shuffle_vectors_w4_ver[shuf_offset]);

  // For a 4 width block, height must be at least 4. Handle 4 lines at once.
  for (int y = 0; y < height; y += 4) {
    __m128i vleft = _mm_loadu_si128((__m128i*) &ref_side[y + shifted_inv_angle_sum[0] + 1]);
    vleft = _mm_shuffle_epi8(vleft, vshuf);

    //__m128i vdst = _mm_i32gather_epi32((const int32_t*)(dst + y * width), vseq, 4);
    __m128i vdst = _mm_loadu_si128((const __m128i*)(dst + y * width));
    
    __m128i vlo = _mm_unpacklo_epi8(vdst, vleft);
    __m128i vhi = _mm_unpackhi_epi8(vdst, vleft);

    __m128i vmaddlo = _mm_maddubs_epi16(vlo, vweight);
    __m128i vmaddhi = _mm_maddubs_epi16(vhi, vweight);

    vmaddlo = _mm_add_epi16(vmaddlo, v32s);
    vmaddhi = _mm_add_epi16(vmaddhi, v32s);

    vmaddlo = _mm_srai_epi16(vmaddlo, 6);
    vmaddhi = _mm_srai_epi16(vmaddhi, 6);

    __m128i packed = _mm_packus_epi16(vmaddlo, vmaddhi);

    _mm_store_si128((__m128i*)(dst + (y * width)), packed);
  }
}


static void angular_pdpc_ver_4x4_scale0_high_angle_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  // This function is just the w4 function, retrofitted to work with any width when scale is 0. If width is 4, use a specialized function instead.
  // Since scale is 0, limit is 3 and therefore there is no meaningful work to be done when x > 3, so only the first column of 4x4 chunks is handled.
  const int scale = 0;
  const int log2_width = uvg_g_convert_to_log2[width];
  __m128i v32s = _mm_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int offset = scale * 16;
  const __m128i vweight = _mm_load_si128((const __m128i*) & intra_pdpc_w4_ver_improved_weight[offset]);

  const int inv_angle_offset = mode_disp * 64;
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  __m128i vseq = _mm_setr_epi32(0, 1, 2, 3);
  __m128i vidx = _mm_slli_epi32(vseq, log2_width);

  __m128i vidx_left = _mm_setr_epi32(shifted_inv_angle_sum[0], shifted_inv_angle_sum[1],
                                     shifted_inv_angle_sum[2], shifted_inv_angle_sum[3]);

  const __m128i vleftshuf = _mm_setr_epi8(
    0x00, 0x04, 0x08, 0x0c, 0x01, 0x05, 0x09, 0x0d,
    0x02, 0x06, 0x0a, 0x0e, 0x03, 0x07, 0x0b, 0x0f
  );

  // For a 4 width block, height must be at least 4. Handle 4 lines at once.
  for (int y = 0; y < height; y += 4) {
    __m128i vdst = _mm_i32gather_epi32((const int32_t*)(dst + y * width), vidx, 1);
    __m128i vleft = _mm_i32gather_epi32((const int32_t*)&ref_side[y + 1], vidx_left, 1);
    vleft = _mm_shuffle_epi8(vleft, vleftshuf);

    __m128i vlo = _mm_unpacklo_epi8(vdst, vleft);
    __m128i vhi = _mm_unpackhi_epi8(vdst, vleft);

    __m128i vmaddlo = _mm_maddubs_epi16(vlo, vweight);
    __m128i vmaddhi = _mm_maddubs_epi16(vhi, vweight);

    vmaddlo = _mm_add_epi16(vmaddlo, v32s);
    vmaddhi = _mm_add_epi16(vmaddhi, v32s);

    vmaddlo = _mm_srai_epi16(vmaddlo, 6);
    vmaddhi = _mm_srai_epi16(vmaddhi, 6);

    __m128i packed = _mm_packus_epi16(vmaddlo, vmaddhi);

    *(uint32_t*)(dst + (y + 0) * width) = _mm_extract_epi32(packed, 0);
    *(uint32_t*)(dst + (y + 1) * width) = _mm_extract_epi32(packed, 1);
    *(uint32_t*)(dst + (y + 2) * width) = _mm_extract_epi32(packed, 2);
    *(uint32_t*)(dst + (y + 3) * width) = _mm_extract_epi32(packed, 3);
  }
}

static void angular_pdpc_ver_4x4_scale0_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  // This function is just the w4 function, retrofitted to work with any width when scale is 0. If width is 4, use a specialized function instead.
  // Since scale is 0, limit is 3 and therefore there is no meaningful work to be done when x > 3, so only the first column of 4x4 chunks is handled.
  // This function handles cases where prediction angle is high. For PDPC, this means the needed reference samples are close together, enabling more effective loading.
  const int scale = 0;
  const int log2_width = uvg_g_convert_to_log2[width];

  const int limit = 3;

  __m128i vseq = _mm_setr_epi32(0, 1, 2, 3);
  __m128i vidx = _mm_slli_epi32(vseq, log2_width);
  __m128i v32s = _mm_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int offset = scale * 16;
  const int inv_angle_offset = mode_disp * 64;
  const int shuf_offset = mode_disp * 16;

  const __m128i vweight = _mm_load_si128((const __m128i*) &intra_pdpc_w4_ver_improved_weight[offset]);
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];
  const __m128i vshuf = _mm_loadu_si128((__m128i*) &intra_pdpc_shuffle_vectors_w4_ver[shuf_offset]);

  // For a 4 width block, height must be at least 4. Handle 4 lines at once.
  for (int y = 0; y < height; y += 4) {
    __m128i vleft = _mm_loadu_si128((__m128i*) & ref_side[y + shifted_inv_angle_sum[0] + 1]);
    vleft = _mm_shuffle_epi8(vleft, vshuf);
    __m128i vdst = _mm_i32gather_epi32((const int32_t*)(dst + y * width), vidx, 1);

    __m128i vlo = _mm_unpacklo_epi8(vdst, vleft);
    __m128i vhi = _mm_unpackhi_epi8(vdst, vleft);

    __m128i vmaddlo = _mm_maddubs_epi16(vlo, vweight);
    __m128i vmaddhi = _mm_maddubs_epi16(vhi, vweight);

    vmaddlo = _mm_add_epi16(vmaddlo, v32s);
    vmaddhi = _mm_add_epi16(vmaddhi, v32s);

    vmaddlo = _mm_srai_epi16(vmaddlo, 6);
    vmaddhi = _mm_srai_epi16(vmaddhi, 6);

    __m128i packed = _mm_packus_epi16(vmaddlo, vmaddhi);

    *(uint32_t*)(dst + (y + 0) * width) = _mm_extract_epi32(packed, 0);
    *(uint32_t*)(dst + (y + 1) * width) = _mm_extract_epi32(packed, 1);
    *(uint32_t*)(dst + (y + 2) * width) = _mm_extract_epi32(packed, 2);
    *(uint32_t*)(dst + (y + 3) * width) = _mm_extract_epi32(packed, 3);
  }
}


static void angular_pdpc_ver_8x4_scale1_high_angle_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  // NOTE: This function is just the w8 function, retrofitted to work with width 16 and up when scale is 1.
  // Since scale is 1, limit is 6 and therefore there is no meaningful work to be done when x > 6, so only the first column of 8x2 chunks is handled.
  const int scale = 1;
  const int log2_width = uvg_g_convert_to_log2[width];

  const int limit = 6;

  __m256i vseq = _mm256_setr_epi64x(0, 1, 2, 3);
  __m256i vidx = _mm256_slli_epi32(vseq, log2_width);
  __m256i v32s = _mm256_set1_epi16(32);

  const int offset = scale * 32;
  const int inv_angle_offset = mode_disp * 64;
  const int shuf_offset = mode_disp * 16;

  const __m256i vweight = _mm256_load_si256((const __m256i*) &intra_pdpc_w8_ver_improved_weight[offset]);
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];
  //const __m128i vshuf = _mm_loadu_si128((__m128i*) &intra_pdpc_shuffle_vectors_8x2_scale1_ver[shuf_offset]);

  __m256i vidxleft = _mm256_setr_epi32(shifted_inv_angle_sum[0], shifted_inv_angle_sum[1],
                                       shifted_inv_angle_sum[2], shifted_inv_angle_sum[3],
                                       shifted_inv_angle_sum[4], shifted_inv_angle_sum[5],
                                       shifted_inv_angle_sum[6], shifted_inv_angle_sum[7]); // These two are not needed.

  const __m256i vtranspose0 = _mm256_setr_epi8(
    0x00, 0x04, 0x08, 0x0c, 0x01, 0x05, 0x09, 0x0d,
    0x02, 0x06, 0x0a, 0x0e, 0x03, 0x07, 0x0b, 0x0f,
    0x00, 0x04, 0x08, 0x0c, 0x01, 0x05, 0x09, 0x0d,
    0x02, 0x06, 0x0a, 0x0e, 0x03, 0x07, 0x0b, 0x0f
  );

  const __m256i vtranspose1 = _mm256_setr_epi8(
    0x00, 0x01, 0x02, 0x03, 0x08, 0x09, 0x0a, 0x0b,
    0x04, 0x05, 0x06, 0x07, 0x0c, 0x0d, 0x0e, 0x0f,
    0x00, 0x01, 0x02, 0x03, 0x08, 0x09, 0x0a, 0x0b,
    0x04, 0x05, 0x06, 0x07, 0x0c, 0x0d, 0x0e, 0x0f
  );

  // For width 8, height must be at least 4 as PDPC is not done when height < 4. Handle 4 lines at once, this enables us to use gather for ref pixels.
  for (int y = 0; y < height; y += 4) {
    __m256i vdst =  _mm256_i64gather_epi64((const long long int*)(dst + y * width), vidx, 1);
    __m256i vleft = _mm256_i32gather_epi32((const int32_t*)&ref_side[y + 1], vidxleft, 1);
    
    // Transpose vleft
    vleft = _mm256_shuffle_epi8(vleft, vtranspose0);
    vleft = _mm256_permute4x64_epi64(vleft, _MM_SHUFFLE(3, 1, 2, 0));
    vleft = _mm256_shuffle_epi8(vleft, vtranspose1);
    
    __m256i vlo = _mm256_unpacklo_epi8(vdst, vleft);
    __m256i vhi = _mm256_unpackhi_epi8(vdst, vleft);

    __m256i vmaddlo = _mm256_maddubs_epi16(vlo, vweight);
    __m256i vmaddhi = _mm256_maddubs_epi16(vhi, vweight);

    vmaddlo = _mm256_add_epi16(vmaddlo, v32s);
    vmaddhi = _mm256_add_epi16(vmaddhi, v32s);

    vmaddlo = _mm256_srai_epi16(vmaddlo, 6);
    vmaddhi = _mm256_srai_epi16(vmaddhi, 6);

    __m256i packed = _mm256_packus_epi16(vmaddlo, vmaddhi);

    // TODO: if this if branch is deemed to cause slow down, make another version of this, where this check is not needed.
    // If this does not slow down significantly, make this same check in other functions to reduce the function call switch case complexity
    if (width == 8) {
      _mm256_store_si256((__m256i*)(dst + (y * width)), packed);
    }
    else {
      *(uint64_t*)(dst + (y + 0) * width) = _mm256_extract_epi64(packed, 0);
      *(uint64_t*)(dst + (y + 1) * width) = _mm256_extract_epi64(packed, 1);
      *(uint64_t*)(dst + (y + 2) * width) = _mm256_extract_epi64(packed, 2);
      *(uint64_t*)(dst + (y + 3) * width) = _mm256_extract_epi64(packed, 3);
    }
  }
}

static void angular_pdpc_ver_8x4_scale1_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  // NOTE: This function is just the w8 function, retrofitted to work with width 16 and up when scale is 1.
  // Since scale is 1, limit is 6 and therefore there is no meaningful work to be done when x > 6, so only the first column of 8x2 chunks is handled.
  // This function handles cases where prediction angle is high. For PDPC, this means the needed reference samples are close together, enabling more effective loading.
  const int scale = 1;
  const int log2_width = uvg_g_convert_to_log2[width];

  const int limit = 6;

  __m256i vseq = _mm256_setr_epi64x(0, 1, 2, 3);
  __m256i vidx = _mm256_slli_epi32(vseq, log2_width);
  __m256i v32s = _mm256_set1_epi16(32);

  const int offset = scale * 32;
  const int inv_angle_offset = mode_disp * 64;
  const int shuf_offset = mode_disp * 16;

  const __m256i vweight = _mm256_load_si256((const __m256i*) &intra_pdpc_w8_ver_improved_weight[offset]);
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];
  const __m128i vshuf = _mm_loadu_si128((__m128i*) &intra_pdpc_shuffle_vectors_8x2_scale1_ver[shuf_offset]);

  // For width 8, height must be at least 4 as PDPC is not done when height < 4. Handle 4 lines at once.
  for (int y = 0; y < height; y += 4) {
    __m128i vleft0 = _mm_loadu_si128((__m128i*) &ref_side[(y + 0) + shifted_inv_angle_sum[0] + 1]);
    __m128i vleft1 = _mm_loadu_si128((__m128i*) &ref_side[(y + 2) + shifted_inv_angle_sum[0] + 1]);
    vleft0 = _mm_shuffle_epi8(vleft0, vshuf);
    vleft1 = _mm_shuffle_epi8(vleft1, vshuf);

    __m256i vleft = _mm256_inserti128_si256(_mm256_castsi128_si256(vleft0), vleft1, 1);
    __m256i vdst = _mm256_i64gather_epi64((const long long int*)(dst + y * width), vidx, 1);

    __m256i vlo = _mm256_unpacklo_epi8(vdst, vleft);
    __m256i vhi = _mm256_unpackhi_epi8(vdst, vleft);

    __m256i vmaddlo = _mm256_maddubs_epi16(vlo, vweight);
    __m256i vmaddhi = _mm256_maddubs_epi16(vhi, vweight);

    vmaddlo = _mm256_add_epi16(vmaddlo, v32s);
    vmaddhi = _mm256_add_epi16(vmaddhi, v32s);

    vmaddlo = _mm256_srai_epi16(vmaddlo, 6);
    vmaddhi = _mm256_srai_epi16(vmaddhi, 6);

    __m256i packed = _mm256_packus_epi16(vmaddlo, vmaddhi);

    // TODO: if this if branch is deemed to cause slow down, make another version of this, where this check is not needed.
    // If this does not slow down significantly, make this same check in other functions to reduce the function call switch case complexity
    if (width == 8) {
      _mm256_store_si256((__m256i*)(dst + (y * width)), packed);
    }
    else {
      *(uint64_t*)(dst + (y + 0) * width) = _mm256_extract_epi64(packed, 0);
      *(uint64_t*)(dst + (y + 1) * width) = _mm256_extract_epi64(packed, 1);
      *(uint64_t*)(dst + (y + 2) * width) = _mm256_extract_epi64(packed, 2);
      *(uint64_t*)(dst + (y + 3) * width) = _mm256_extract_epi64(packed, 3);
    }
  }
}


static void angular_pdpc_ver_8x2_scale2_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  // NOTE: This function is just the w8 function, retrofitted to work with width 16 and up when scale is 1.
  // Since scale is 1, limit is 6 and therefore there is no meaningful work to be done when x > 6, so only the first column of 8x2 chunks is handled.
  // This function handles cases where prediction angle is high. For PDPC, this means the needed reference samples are close together, enabling more effective loading.
  const int scale = 2;
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

    __m128i vdst = _mm_i64gather_epi64((const long long int*)(dst + y * width), vidx, 1);
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


static void angular_pdpc_ver_w16_high_angle_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  __m256i v32s = _mm256_set1_epi16(32);
  const int scale = 2; // Other functions handle scales 0 and 1
  int limit = 12; // With scale 2, limit is always 12.

  const int offset = scale * 32;
  const __m256i vweight = _mm256_load_si256((const __m256i*) &intra_pdpc_w16_ver_improved_weight[offset]);

  const int inv_angle_offset = mode_disp * 64;
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  const __m256i vidx0 = _mm256_setr_epi32(shifted_inv_angle_sum[0], shifted_inv_angle_sum[1],
                                          shifted_inv_angle_sum[2], shifted_inv_angle_sum[3],
                                          shifted_inv_angle_sum[4], shifted_inv_angle_sum[5],
                                          shifted_inv_angle_sum[6], shifted_inv_angle_sum[7]);
  const __m256i vidx1 = _mm256_setr_epi32(shifted_inv_angle_sum[8], shifted_inv_angle_sum[9],
                                          shifted_inv_angle_sum[10], shifted_inv_angle_sum[11],
                                          shifted_inv_angle_sum[12], shifted_inv_angle_sum[13],  // These are not used.
                                          shifted_inv_angle_sum[14], shifted_inv_angle_sum[15]); // These are not used.

  const __m256i transpose = _mm256_setr_epi8(
    0x00, 0x04, 0x08, 0x0c, 0x01, 0x05, 0x09, 0x0d,
    0x02, 0x06, 0x0a, 0x0e, 0x03, 0x07, 0x0b, 0x0f,
    0x00, 0x04, 0x08, 0x0c, 0x01, 0x05, 0x09, 0x0d,
    0x02, 0x06, 0x0a, 0x0e, 0x03, 0x07, 0x0b, 0x0f
  );

  // 0xff are 'don't care' values, they will be zeroed out by coefficients
  const __m256i vpermute = _mm256_setr_epi32(
    0x00, 0x04, 0x02, 0xff, 0x01, 0x05, 0x03, 0xff
  );

  // Handle 4 rows at once to enable gather for ref pixels.
  for (int y = 0; y < height; y += 4) {
    __m128i vdstraw0 = _mm_load_si128((const __m128i*)(dst + ((y + 0) * width)));
    __m128i vdstraw1 = _mm_load_si128((const __m128i*)(dst + ((y + 1) * width)));
    __m128i vdstraw2 = _mm_load_si128((const __m128i*)(dst + ((y + 2) * width)));
    __m128i vdstraw3 = _mm_load_si128((const __m128i*)(dst + ((y + 3) * width)));

    __m256i vdst0 = _mm256_inserti128_si256(_mm256_castsi128_si256(vdstraw0), vdstraw1, 1);
    __m256i vdst1 = _mm256_inserti128_si256(_mm256_castsi128_si256(vdstraw2), vdstraw3, 1);

    __m256i vleft0 = _mm256_i32gather_epi32((const int32_t*)&ref_side[y + 1], vidx0, 1);
    __m256i vleft1 = _mm256_i32gather_epi32((const int32_t*)&ref_side[y + 1], vidx1, 1);
    vleft0 = _mm256_shuffle_epi8(vleft0, transpose);
    vleft1 = _mm256_shuffle_epi8(vleft1, transpose);

    __m256i vtmplo = _mm256_unpacklo_epi64(vleft0, vleft1);
    __m256i vtmphi = _mm256_unpackhi_epi64(vleft0, vleft1);

    vleft0 = _mm256_permutevar8x32_epi32(vtmplo, vpermute);
    vleft1 = _mm256_permutevar8x32_epi32(vtmphi, vpermute);

    __m256i vlo0 = _mm256_unpacklo_epi8(vdst0, vleft0);
    __m256i vhi0 = _mm256_unpackhi_epi8(vdst0, vleft0);
    __m256i vlo1 = _mm256_unpacklo_epi8(vdst1, vleft1);
    __m256i vhi1 = _mm256_unpackhi_epi8(vdst1, vleft1);

    __m256i v0 = _mm256_permute2x128_si256(vlo0, vhi0, 0x20);
    __m256i v1 = _mm256_permute2x128_si256(vlo0, vhi0, 0x31);
    __m256i v2 = _mm256_permute2x128_si256(vlo1, vhi1, 0x20);
    __m256i v3 = _mm256_permute2x128_si256(vlo1, vhi1, 0x31);

    __m256i vmadd0 = _mm256_maddubs_epi16(v0, vweight);
    __m256i vmadd1 = _mm256_maddubs_epi16(v1, vweight);
    __m256i vmadd2 = _mm256_maddubs_epi16(v2, vweight);
    __m256i vmadd3 = _mm256_maddubs_epi16(v3, vweight);

    vmadd0 = _mm256_add_epi16(vmadd0, v32s);
    vmadd1 = _mm256_add_epi16(vmadd1, v32s);
    vmadd2 = _mm256_add_epi16(vmadd2, v32s);
    vmadd3 = _mm256_add_epi16(vmadd3, v32s);

    vmadd0 = _mm256_srai_epi16(vmadd0, 6);
    vmadd1 = _mm256_srai_epi16(vmadd1, 6);
    vmadd2 = _mm256_srai_epi16(vmadd2, 6);
    vmadd3 = _mm256_srai_epi16(vmadd3, 6);

    __m256i packed0 = _mm256_packus_epi16(vmadd0, vmadd1);
    __m256i packed1 = _mm256_packus_epi16(vmadd2, vmadd3);
    packed0 = _mm256_permute4x64_epi64(packed0, _MM_SHUFFLE(3, 1, 2, 0));
    packed1 = _mm256_permute4x64_epi64(packed1, _MM_SHUFFLE(3, 1, 2, 0));

    if (width == 16) {
      _mm256_store_si256((__m256i*)(dst + ((y + 0) * width)), packed0);
      _mm256_store_si256((__m256i*)(dst + ((y + 2) * width)), packed1);
    }
    else {
      _mm_store_si128((__m128i*)(dst + ((y + 0) * width)), _mm256_extracti128_si256(packed0, 0));
      _mm_store_si128((__m128i*)(dst + ((y + 1) * width)), _mm256_extracti128_si256(packed0, 1));
      _mm_store_si128((__m128i*)(dst + ((y + 2) * width)), _mm256_extracti128_si256(packed1, 0));
      _mm_store_si128((__m128i*)(dst + ((y + 3) * width)), _mm256_extracti128_si256(packed1, 1));
    }
  }
}

static void angular_pdpc_ver_w16_scale2_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int mode_disp)
{
  __m128i v32s = _mm_set1_epi16(32);
  const int scale = 2; // Other functions handle scales 0 and 1
  int limit = 12; // With scale 2, limit is always 12.

  const int offset = scale * 32;
  const int inv_angle_offset = mode_disp * 64;
  const int shuf_offset = mode_disp * 16;

  const __m128i vweightlo = _mm_load_si128((const __m128i*) &intra_pdpc_w16_ver_improved_weight[offset + 0]);
  const __m128i vweighthi = _mm_load_si128((const __m128i*) &intra_pdpc_w16_ver_improved_weight[offset + 16]);
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];
  const __m128i vshuf = _mm_load_si128((const __m128i*) & intra_pdpc_shuffle_vectors_w16_scale2_ver[shuf_offset]);

  // Handle 2 rows at once.
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < limit; x += 16) {
      __m128i vleft = _mm_loadu_si128((__m128i*) &ref_side[(y + 0) + shifted_inv_angle_sum[0] + 1]);
      vleft = _mm_shuffle_epi8(vleft, vshuf);
      
      __m128i vdst = _mm_load_si128((const __m128i*)(dst + ((y + 0) * width + x)));

      __m128i vlo = _mm_unpacklo_epi8(vdst, vleft);
      __m128i vhi = _mm_unpackhi_epi8(vdst, vleft);

      __m128i vmaddlo = _mm_maddubs_epi16(vlo, vweightlo);
      __m128i vmaddhi = _mm_maddubs_epi16(vhi, vweighthi);

      vmaddlo = _mm_add_epi16(vmaddlo, v32s);
      vmaddhi = _mm_add_epi16(vmaddhi, v32s);

      vmaddlo = _mm_srai_epi16(vmaddlo, 6);
      vmaddhi = _mm_srai_epi16(vmaddhi, 6);

      __m128i packed = _mm_packus_epi16(vmaddlo, vmaddhi);

      _mm_store_si128((__m128i*)(dst + (y * width + x)), packed);
    }
  }
}


// Horizontal modes

static void angular_pdpc_hor_w4_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int height, const int scale, const int mode_disp)
{
  const int width = 4;

  int limit = MIN(3 << scale, height);
  
  __m128i v32s = _mm_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int table_offset = scale * 128;
  const int shuf_offset = mode_disp * 256;
  const int inv_angle_offset = mode_disp * 64;
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  for (int y = 0, so = 0, wo = 0; y < limit; y += 4, so += 16, wo += 32) {
    const __m128i vshuf = _mm_loadu_si128((__m128i*) & intra_pdpc_shuffle_vectors_w4_hor[shuf_offset + so]);

    __m128i vtop = _mm_loadu_si128((__m128i*) & ref_side[shifted_inv_angle_sum[y] + 1]);
    vtop = _mm_shuffle_epi8(vtop, vshuf);

    const int offset = table_offset + wo;

    __m128i vdst = _mm_load_si128((const __m128i*)(dst + y * width));
    __m128i vweightlo = _mm_load_si128((const __m128i*) &intra_pdpc_w4_hor_improved_weight[offset + 0]);
    __m128i vweighthi = _mm_load_si128((const __m128i*) &intra_pdpc_w4_hor_improved_weight[offset + 16]);

    __m128i vlo = _mm_unpacklo_epi8(vdst, vtop);
    __m128i vhi = _mm_unpackhi_epi8(vdst, vtop);

    __m128i vmaddlo = _mm_maddubs_epi16(vlo, vweightlo);
    __m128i vmaddhi = _mm_maddubs_epi16(vhi, vweighthi);

    vmaddlo = _mm_add_epi16(vmaddlo, v32s);
    vmaddhi = _mm_add_epi16(vmaddhi, v32s);

    vmaddlo = _mm_srai_epi16(vmaddlo, 6);
    vmaddhi = _mm_srai_epi16(vmaddhi, 6);

    __m128i packed = _mm_packus_epi16(vmaddlo, vmaddhi);

    _mm_storeu_si128((__m128i*)(dst + (y * width)), packed);
  }
}

static void angular_pdpc_hor_w8_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int height, const int scale, const int mode_disp)
{
  const int width = 8;

  int limit = MIN(3 << scale, height);

  __m256i v32s = _mm256_set1_epi16(32);

  // Scale can be 0, 1 or 2
  const int table_offset = scale * 256;
  const int inv_angle_offset = mode_disp * 64;
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  // Handle 4 lines at once since PDPC is not done on 8x2 blocks.
  for (int y = 0, o = table_offset; y < limit; y += 4, o += 64) {
    const __m256i vweight01 = _mm256_load_si256((const __m256i*) &intra_pdpc_w8_hor_improved_weight[o + 0]);
    const __m256i vweight23 = _mm256_load_si256((const __m256i*) &intra_pdpc_w8_hor_improved_weight[o + 32]);

    const __m256i vidx = _mm256_set_epi64x(shifted_inv_angle_sum[y + 3], shifted_inv_angle_sum[y + 2],
                                           shifted_inv_angle_sum[y + 1], shifted_inv_angle_sum[y + 0]);

    __m256i vdst = _mm256_load_si256((const __m256i*)(dst + y * width));
    __m256i vtop = _mm256_i64gather_epi64((const long long int*)&ref_side[1], vidx, 1);
    
    __m256i vlo = _mm256_unpacklo_epi8(vdst, vtop);
    __m256i vhi = _mm256_unpackhi_epi8(vdst, vtop);

    __m256i v01 = _mm256_permute2x128_si256(vlo, vhi, 0x20);
    __m256i v23 = _mm256_permute2x128_si256(vlo, vhi, 0x31);

    __m256i vmadd01 = _mm256_maddubs_epi16(v01, vweight01);
    __m256i vmadd23 = _mm256_maddubs_epi16(v23, vweight23);

    vmadd01 = _mm256_add_epi16(vmadd01, v32s);
    vmadd23 = _mm256_add_epi16(vmadd23, v32s);

    vmadd01 = _mm256_srai_epi16(vmadd01, 6);
    vmadd23 = _mm256_srai_epi16(vmadd23, 6);

    __m256i packed = _mm256_packus_epi16(vmadd01, vmadd23);
    packed = _mm256_permute4x64_epi64(packed, _MM_SHUFFLE(3, 1, 2, 0));

    _mm256_storeu_si256((__m256i*)(dst + (y * width)), packed);
  }
}

static void angular_pdpc_hor_w16_avx2(uvg_pixel* dst, const uvg_pixel* ref_side, const int width, const int height, const int scale, const int mode_disp)
{
  int limit = MIN(3 << scale, height);
  __m128i v32s = _mm_set1_epi16(32);

  const int inv_angle_offset = mode_disp * 64;
  const int16_t* shifted_inv_angle_sum = &intra_pdpc_shifted_inv_angle_sum[inv_angle_offset];

  // Handle one line at a time. Skip line if vertical limit reached.
  for (int y = 0; y < limit; ++y) {
    const uint8_t weight1 = 32 >> (2 * y >> scale);
    const uint8_t weight0 = 64 - weight1;
    ALIGNED(2) const uint8_t tmp[2] = { weight0, weight1 };

    __m128i vweight = _mm_set1_epi16(*(uint16_t*)tmp);

    for (int x = 0; x < width; x += 16) {
      __m128i vdst = _mm_load_si128((__m128i*)(dst + (y * width + x)));
      __m128i vtop = _mm_loadu_si128((__m128i*) &ref_side[x + shifted_inv_angle_sum[y] + 1]);
      
      __m128i vlo = _mm_unpacklo_epi8(vdst, vtop);
      __m128i vhi = _mm_unpackhi_epi8(vdst, vtop);

      __m128i vmaddlo = _mm_maddubs_epi16(vlo, vweight);
      __m128i vmaddhi = _mm_maddubs_epi16(vhi, vweight);

      vmaddlo = _mm_add_epi16(vmaddlo, v32s);
      vmaddhi = _mm_add_epi16(vmaddhi, v32s);

      vmaddlo = _mm_srai_epi16(vmaddlo, 6);
      vmaddhi = _mm_srai_epi16(vmaddhi, 6);

      __m128i packed = _mm_packus_epi16(vmaddlo, vmaddhi);

      _mm_storeu_si128((__m128i*)(dst + (y * width + x)), packed);
    }
  }
}


// Prediction mode 50 versions of PDPC functions.
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
    __m128i vdst = _mm_i64gather_epi64((const long long int*)(dst + y * width), vidx, 1);
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

// The main angular prediction entry point for AVX2.
/**
 * \brief AVX2 version of angular intra prediction.
 * \param cu_loc        CU location and size data.
 * \param intra_mode    Intra prediction mode.
 * \param channel_type  Color channel.
 * \param in_ref_above  Pointer to -1 index of above reference.
 * \param in_ref_left   Pointer to -1 index of left reference.
 * \param dst           Buffer of size MAX_PRED_WIDTH * MAX_PRED_WIDTH.
 * \param multi_ref_idx Multi reference index.
 * \param isp_mode      Intra sub-partition mode.
 * \param cu_dim        CU dimension, used along ISP mode.
 */
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
  uvg_pixel temp_main[2 * 128 + 3 + 33 * MAX_REF_LINE_IDX];
  uvg_pixel temp_side[2 * 128 + 3 + 33 * MAX_REF_LINE_IDX];

  int32_t pred_mode = intra_mode; // ToDo: handle WAIP

  // Whether to swap references to always project on the left reference row.
  const bool vertical_mode = intra_mode >= 34;
  
  // Modes distance to horizontal or vertical mode. Possible values: [-16, 16]
  // For pure vertical or horizontal modes, this is 0. For pure diagonal modes, this is either -16 or 16.
  const int_fast8_t mode_disp = vertical_mode ? pred_mode - 50 : -(pred_mode - 18);
  const int_fast8_t abs_mode_disp = abs(mode_disp);
  const bool wide_angle_mode = mode_disp > 16;

  // Sample displacement per column in fractions of 32.
  const int_fast16_t sample_disp = (mode_disp < 0 ? -1 : 1) * modedisp2sampledisp[abs_mode_disp];

  const int side_size = vertical_mode ? log2_height : log2_width;
  int scale = MIN(2, side_size - pre_scale[abs_mode_disp]);

  // Pointer for the reference we are interpolating from.
  uvg_pixel* ref_main;
  // Pointer for the other reference.
  const uvg_pixel* ref_side;

  const int top_ref_length = isp_mode == ISP_MODE_VER ? width + cu_dim : width << 1;
  const int left_ref_length = isp_mode == ISP_MODE_HOR ? height + cu_dim : height << 1;

  // Set ref_main and ref_side such that, when indexed with 0, they point to
  // index 0 in block coordinates.
  if (sample_disp < 0) {
    // In cases where sample_disp is negative, references are needed from both sides.
    // This step combines the main and side reference.
    if (vertical_mode) {
      memcpy(&temp_main[height], in_ref_above, (width + 2 + multi_ref_index) * sizeof(uvg_pixel));
    }
    else {
      memcpy(&temp_main[width], in_ref_left, (height + 2 + multi_ref_index) * sizeof(uvg_pixel));
    }
    //memcpy(&temp_main[height], &in_ref_above[0], (width + 2 + multi_ref_index) * sizeof(uvg_pixel));
    //memcpy(&temp_side[width], &in_ref_left[0], (height + 2 + multi_ref_index) * sizeof(uvg_pixel));

    ref_main = vertical_mode ? &temp_main[height] : &temp_main[width];
    ref_side = vertical_mode ? in_ref_left : in_ref_above;

    int size_side = vertical_mode ? height : width;
    switch (size_side) {
      case 4:
      {
        int shuf_offset = abs_mode_disp * 16;
        __m128i vshuf = _mm_load_si128((__m128i*) &intra_refbuild_shuffle_vectors_sidesize_4[shuf_offset]);
        __m128i vref = _mm_loadu_si128((const __m128i*) &ref_side[0]);
        vref = _mm_shuffle_epi8(vref, vshuf);
        /*uint32_t tmp = _mm_extract_epi32(vref, 0);
        memcpy(&temp_main[0], &tmp, sizeof(uint32_t));*/
        _mm_maskstore_epi32((int32_t*)&temp_main[0], _mm_setr_epi32(0xffffffff, 0, 0, 0), vref);
        break;
      }
      case 8:
      {
        int shuf_offset = abs_mode_disp * 16;
        __m128i vshuf = _mm_load_si128((__m128i*) &intra_refbuild_shuffle_vectors_sidesize_8[shuf_offset]);
        __m128i vref = _mm_loadu_si128((const __m128i*) &ref_side[0]);
        vref = _mm_shuffle_epi8(vref, vshuf);
        /*uint64_t tmp = _mm_extract_epi64(vref, 0);
        memcpy(&temp_main[0], &tmp, sizeof(uint64_t));*/
        _mm_maskstore_epi32((int32_t*)&temp_main[0], _mm_setr_epi32(0xffffffff, 0xffffffff, 0, 0), vref);
        break;
      }
      case 16:
      {
        int shuf_offset = abs_mode_disp * 16;
        __m128i vshuf = _mm_load_si128((__m128i*) &intra_refbuild_shuffle_vectors_sidesize_16[shuf_offset]);
        __m128i vref = _mm_loadu_si128((const __m128i*) &ref_side[1]); // Offset ref by one to fit all necessary 16 refs. Offset accounted for in shuffle vectors.
        vref = _mm_shuffle_epi8(vref, vshuf);
        _mm_store_si128((__m128i*) &temp_main[0], vref);
        break;
      }
      case 32:
      {
        int shuf_offset = abs_mode_disp * 32;
        __m128i vshufhi = _mm_load_si128((__m128i*) &intra_refbuild_shuffle_vectors_sidesize_32[shuf_offset + 0]);
        __m128i vshuflo = _mm_load_si128((__m128i*) &intra_refbuild_shuffle_vectors_sidesize_32[shuf_offset + 16]);
        __m128i vblend = _mm_cmpgt_epi8(vshuflo, _mm_set1_epi8(15));
        
        __m128i vreflo = _mm_loadu_si128((const __m128i*) & ref_side[1]); // Offset ref by one to fit all necessary 16 refs. Offset accounted for in shuffle vectors.
        __m128i vrefhi = _mm_loadu_si128((const __m128i*) & ref_side[17]);
        
        // Second half of references requires samples from both sides
        __m128i vreftmphi = _mm_shuffle_epi8(vrefhi, vshuflo);
        __m128i vreftmplo = _mm_shuffle_epi8(vreflo, vshuflo);
        vreflo = _mm_blendv_epi8(vreftmplo, vreftmphi, vblend);
        
        // First half of references use references from the hi side only
        vrefhi = _mm_shuffle_epi8(vrefhi, vshufhi);
        
        _mm_store_si128((__m128i*) &temp_main[0], vrefhi);
        _mm_store_si128((__m128i*) &temp_main[16], vreflo);
        break;
      }
      case 64:
      {
        int shuf_offset = abs_mode_disp * 64;
        __m128i vshuf0 = _mm_load_si128((__m128i*) &intra_refbuild_shuffle_vectors_sidesize_64[shuf_offset + 0]);
        __m128i vshuf1 = _mm_load_si128((__m128i*) &intra_refbuild_shuffle_vectors_sidesize_64[shuf_offset + 16]);
        __m128i vshuf2 = _mm_load_si128((__m128i*) &intra_refbuild_shuffle_vectors_sidesize_64[shuf_offset + 32]);
        __m128i vshuf3 = _mm_load_si128((__m128i*) &intra_refbuild_shuffle_vectors_sidesize_64[shuf_offset + 48]);

        __m128i vref0 = _mm_loadu_si128((const __m128i*) &ref_side[ 0 + 1]); // Offset ref by one to fit all necessary 16 refs. Offset accounted for in shuffle vectors.
        __m128i vref1 = _mm_loadu_si128((const __m128i*) &ref_side[16 + 1]);
        __m128i vref2 = _mm_loadu_si128((const __m128i*) &ref_side[32 + 1]);
        __m128i vref3 = _mm_loadu_si128((const __m128i*) &ref_side[48 + 1]);

        // First quarter of references use references from vref3 only
        __m128i vrefout0 = _mm_shuffle_epi8(vref3, vshuf0);

        // Second quarter can require samples from vref3 and vref2
        __m128i vreftmp0 = _mm_shuffle_epi8(vref3, vshuf1);
        __m128i vreftmp1 = _mm_shuffle_epi8(vref2, vshuf1);
        __m128i vblend0 = _mm_cmpgt_epi8(vshuf1, _mm_set1_epi8(47));
        __m128i vrefout1 = _mm_blendv_epi8(vreftmp1, vreftmp0, vblend0);

        // Third quarter can require samples from vref3, vref2 and vref1
        vreftmp0 = _mm_shuffle_epi8(vref3, vshuf2);
        vreftmp1 = _mm_shuffle_epi8(vref2, vshuf2);
        __m128i vreftmp2 = _mm_shuffle_epi8(vref1, vshuf2);
        vblend0 = _mm_cmpgt_epi8(vshuf2, _mm_set1_epi8(47));
        __m128i vblend1 = _mm_cmpgt_epi8(vshuf2, _mm_set1_epi8(31));

        vreftmp0 = _mm_blendv_epi8(vreftmp1, vreftmp0, vblend0);
        __m128i vrefout2 = _mm_blendv_epi8(vreftmp2, vreftmp0, vblend1);

        // Fourth quarter can require samples from vref3, vref2, vref1 and vref0
        vreftmp0 = _mm_shuffle_epi8(vref3, vshuf3);
        vreftmp1 = _mm_shuffle_epi8(vref2, vshuf3);
        vreftmp2 = _mm_shuffle_epi8(vref1, vshuf3);
        __m128i vreftmp3 = _mm_shuffle_epi8(vref0, vshuf3);

        vblend0 = _mm_cmpgt_epi8(vshuf3, _mm_set1_epi8(47));
        vblend1 = _mm_cmpgt_epi8(vshuf3, _mm_set1_epi8(31));
        __m128i vblend2 = _mm_cmpgt_epi8(vshuf3, _mm_set1_epi8(15));

        vreftmp0 = _mm_blendv_epi8(vreftmp1, vreftmp0, vblend0);
        vreftmp0 = _mm_blendv_epi8(vreftmp2, vreftmp0, vblend1);
        __m128i vrefout3 = _mm_blendv_epi8(vreftmp3, vreftmp0, vblend2);

        _mm_store_si128((__m128i*) &temp_main[0],  vrefout0);
        _mm_store_si128((__m128i*) &temp_main[16], vrefout1);
        _mm_store_si128((__m128i*) &temp_main[32], vrefout2);
        _mm_store_si128((__m128i*) &temp_main[48], vrefout3);
        break;
      }
      default:
        // This should work in the case everything else fails.
        for (int i = 0; i < size_side; ++i) {
          const int modedisp2invsampledisp_abs = modedisp2invsampledisp[abs_mode_disp];
          ref_main[i] = ref_side[MIN((-i * modedisp2invsampledisp_abs + 256) >> 9, size_side)];
        }
    }
  }
  else {
    ref_main = (uvg_pixel*)(vertical_mode ? in_ref_above : in_ref_left);
    ref_side = vertical_mode ? in_ref_left : in_ref_above;
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

  //const int8_t* pfilter = use_cubic ? &cubic_filter_8bit_c[0][0] : &cubic_filter_8bit_g[0][0];
  const int8_t (*pfilter)[4] = use_cubic ? cubic_filter_8bit_c : cubic_filter_8bit_g;



  if (sample_disp != 0) {
    // The mode is not horizontal or vertical, we have to do interpolation.

    // Set delta table pointers
    const int table_offset = wide_angle_mode ? (pred_mode < 2 ? (pred_mode + 13) * DELTA_TABLE_ROW_LENGTH : (81 - pred_mode) * DELTA_TABLE_ROW_LENGTH) : (pred_mode <= 34 ? (pred_mode - 2) * DELTA_TABLE_ROW_LENGTH : (66 - pred_mode) * DELTA_TABLE_ROW_LENGTH);
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
            case  4: angular_pred_w4_ver_avx2(dst, ref_main, delta_int, delta_fract, height, pfilter); break;
            case  8: 
              if (height < 4)
                angular_pred_w8_h2_ver_avx2(dst, ref_main, delta_int, delta_fract, height, pfilter);
              else
                angular_pred_w8_ver_avx2(dst, ref_main, delta_int, delta_fract, height, pfilter);
              break;
            case 16: // Use w16 function for all widths 16 and up
            case 32: 
            case 64: angular_pred_w16_ver_avx2(dst, ref_main, delta_int, delta_fract, width, height, pfilter); break;
            default:
              assert(false && "Intra angular predicion: illegal width.\n");
              break;
          }
        }
        else {
          switch (width) {
            case  4: 
              if (pred_mode < -7 || (multi_ref_index == 2 && pred_mode == -7)) // High angles need special handling
                angular_pred_w4_hor_high_angle_avx2(dst, ref_main, delta_int, delta_fract, height, pfilter); 
              else
                angular_pred_w4_hor_avx2(dst, ref_main, pred_mode, multi_ref_index, delta_int, delta_fract, height, pfilter);
              
              break;
            case  8: 
              if (pred_mode < -2)
                angular_pred_w8_hor_high_angle_avx2(dst, ref_main, delta_int, delta_fract, height, pfilter);
              else
                angular_pred_w8_hor_avx2(dst, ref_main, pred_mode, multi_ref_index, delta_int, delta_fract, height, pfilter);

              break;
            case 16: 
              if (pred_mode < 5 || pred_mode == 33)
                angular_pred_w16_hor_high_angle_avx2(dst, ref_main, delta_int, delta_fract, width, height, pfilter);
              else
                angular_pred_w16_hor_avx2(dst, ref_main, pred_mode, multi_ref_index, delta_int, delta_fract, height, pfilter);
              
              break;
            case 32: 
              if (pred_mode < 5 || pred_mode == 33)
                angular_pred_w16_hor_high_angle_avx2(dst, ref_main, delta_int, delta_fract, width, height, pfilter);
              else
                angular_pred_w32_hor_avx2(dst, ref_main, pred_mode, multi_ref_index, delta_int, delta_fract, width, height, pfilter);

              break;
            case 64: 
              if (pred_mode < 5 || pred_mode == 33)
                angular_pred_w16_hor_high_angle_avx2(dst, ref_main, delta_int, delta_fract, width, height, pfilter);
              else
                angular_pred_w32_hor_avx2(dst, ref_main, pred_mode, multi_ref_index, delta_int, delta_fract, width, height, pfilter);

              break;
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
        if (pred_mode == 2) {
          switch (width) {
            // Note: these functions do not need the delta int table as the mode is known
            case 4: angular_pred_non_fractional_angle_pxl_copy_w4_mode2_hor_avx2(dst, ref_main, height, multi_ref_index); break;
            case 8: angular_pred_non_fractional_angle_pxl_copy_w8_mode2_hor_avx2(dst, ref_main, height, multi_ref_index); break;
            // Cases 16 onward can be solved with a simple memcpy
            case 16:
              for (int y = 0; y < height; ++y) {
                // Offset indices by one since index 0 is top left and plus one since delta_int[0] for mode 2 is 1.
                memcpy(&dst[y * 16], &ref_main[2 + y + multi_ref_index], 16 * sizeof(uvg_pixel));
              }
              break;
            case 32:
              for (int y = 0; y < height; ++y) {
                memcpy(&dst[y * 32], &ref_main[2 + y + multi_ref_index], 32 * sizeof(uvg_pixel));
              }
              break;
            case 64:
              for (int y = 0; y < height; ++y) {
                memcpy(&dst[y * 64], &ref_main[2 + y + multi_ref_index], 64 * sizeof(uvg_pixel));
              }
              break;
            default:
              assert(false && "Intra angular predicion: illegal width.\n");
              break;
          }

        }
        else {
          // Wide angle modes -12, -10, -8 and -4
          switch (width) {
            case 4: angular_pred_non_fractional_angle_pxl_copy_w4_wide_angle_hor_avx2(dst, ref_main, height, delta_int); break;
            case 8: angular_pred_non_fractional_angle_pxl_copy_w8_wide_angle_hor_avx2(dst, ref_main, height, delta_int); break;
            case 16: angular_pred_non_fractional_angle_pxl_copy_w16_wide_angle_hor_avx2(dst, ref_main, height, delta_int); break;
            case 32: angular_pred_non_fractional_angle_pxl_copy_w32_wide_angle_hor_avx2(dst, ref_main, height, delta_int); break;
            // Width 64 never goes into this branch. Leave an assert here to catch future problems.
            case 64: 
              //angular_pred_non_fractional_angle_pxl_copy_hor_avx2(dst, ref_main, width, height, delta_int); break;
              assert(false && "Intra angular predicion: Non fractional angle pixel copy with width 64. This should never happen.\n");
              break;
            default:
              assert(false && "Intra angular predicion: illegal width.\n");
              break;
          }
        }
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

  
  bool PDPC_filter = (width >= TR_MIN_WIDTH && height >= TR_MIN_WIDTH && multi_ref_index == 0);
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
          // Low mode disp -> high angle. For pdpc, this causes the needed references to be extremely sparse making loads without using gathers impossible.
          // Handle low angles with more tight reference spacing with separate functions with more optimized loads.
          if (mode_disp < 6)
            angular_pdpc_ver_w4_high_angle_avx2(dst, ref_side, height, scale, mode_disp);
          else
            angular_pdpc_ver_w4_avx2(dst, ref_side, height, scale, mode_disp);
          break;
        case 8:
          if (scale == 0) {
            if (mode_disp < 6)
              angular_pdpc_ver_4x4_scale0_high_angle_avx2(dst, ref_side, width, height, mode_disp);
            else
              angular_pdpc_ver_4x4_scale0_avx2(dst, ref_side, width, height, mode_disp);
          }
          else if (scale == 1) {
            if (mode_disp < 8)
              angular_pdpc_ver_8x4_scale1_high_angle_avx2(dst, ref_side, width, height, mode_disp);
            else
              angular_pdpc_ver_8x4_scale1_avx2(dst, ref_side, width, height, mode_disp);
          }
          else {
            if (mode_disp < 10)
              angular_pdpc_ver_w8_high_angle_avx2(dst, ref_side, height, 2, mode_disp);
            else
              angular_pdpc_ver_8x2_scale2_avx2(dst, ref_side, width, height, mode_disp);
          }
          break;
        case 16: // 16 width and higher done with the same functions
        case 32:
        case 64:
          switch (scale) {
          case 0:
            if (mode_disp < 6)
              angular_pdpc_ver_4x4_scale0_high_angle_avx2(dst, ref_side, width, height, mode_disp);
            else
              angular_pdpc_ver_4x4_scale0_avx2(dst, ref_side, width, height, mode_disp);
            break;
          case 1:
            if (mode_disp < 8)
              angular_pdpc_ver_8x4_scale1_high_angle_avx2(dst, ref_side, width, height, mode_disp);
            else
              angular_pdpc_ver_8x4_scale1_avx2(dst, ref_side, width, height, mode_disp);
            break;
          case 2:
            if (mode_disp < 14)
              angular_pdpc_ver_w16_high_angle_avx2(dst, ref_side, width, height, mode_disp);
            else
              angular_pdpc_ver_w16_scale2_avx2(dst, ref_side, width, height, mode_disp);
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
          // Low mode disp -> high angle. For pdpc, this causes the needed references to be extremely sparse making loads without using gathers impossible.
          // Handle low angles with more tight reference spacing with separate functions with more optimized loads.
          /*if (mode_disp < 6)
            angular_pdpc_hor_w4_high_angle_improved_avx2(dst, ref_side, height, scale, mode_disp);
          else*/
          // The above code was not accessed ever. There is no case where width == 4 and and mode disp < 6 for horizontal modes where PDPC is enabled.
          angular_pdpc_hor_w4_avx2(dst, ref_side, height, scale, mode_disp);
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


/**
* \brief Position Dependent Prediction Combination for Planar and DC modes.
* \param mode          Intra mode, 0 for planar, 1 for DC.
* \param cu_loc        Pointer to the CU location information.
* \param color         Color component.
* \param used_ref      Pointer to the used reference pixels.
* \param dst           Buffer of size MAX_PRED_WIDTH * MAX_PRED_WIDTH.
*/
// TODO: allegedly does not work with blocks with height 1 and 2. Test this.
// TODO: or just rework the whole thing. We might be able to optimize this further.
static void uvg_pdpc_planar_dc_avx2(
  const int mode,
  const cu_loc_t* const cu_loc,
  const color_t color,
  const uvg_intra_ref *const used_ref,
  uvg_pixel *const dst)
{
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

  const int scale = ((log2_width - 2 + log2_height - 2 + 2) >> 2);

  // Same weights regardless of axis, compute once
  int16_t w[LCU_WIDTH];
  for (int i = 0; i < MAX(width, height); i += 4) {
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

  ALIGNED(16) uint32_t ref[4];
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
    weight += input_size * 4;
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
    __m256i vbehind = _mm256_loadu_si256((__m256i*)src_ptr);

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

  _mm256_storeu_si256((__m256i*)(dst + 0), vres0);
  _mm256_storeu_si256((__m256i*)(dst + 32), vres1);
  
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
  __m128i vbefore1 = _mm_loadu_si128((__m128i*)(src + 8));
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

  vbefore0 = _mm_loadu_si128((__m128i*)(src + 24));
  vbehind0 = _mm_load_si128((__m128i*)(src + 32));
  vbefore1 = _mm_loadu_si128((__m128i*)(src + 40));
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


static void mip_upsampling_w32_ups2_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
{
  __m256i vbefore = _mm256_loadu_si256((__m256i*)ref);

  for (int i = 0; i < 8; ++i) {
    __m256i vbehind = _mm256_loadu_si256((__m256i*)(src + (i * 64)));
    __m256i vavg = _mm256_avg_epu8(vbefore, vbehind);

    _mm256_storeu_si256((__m256i*)(dst + (i * 64)), vavg);

    vbefore = vbehind;
  }
}

static void mip_upsampling_w32_ups4_ver_avx2(uvg_pixel* const dst, const uvg_pixel* const src, const uvg_pixel* const ref)
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
          mip_upsampling_w8_ups2_hor_avx2(hor_dst, reduced_pred, ref_samples_left, ver_src_step, ups_ver_factor);
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
            mip_upsampling_w32_ups4_hor_avx2(hor_dst, reduced_pred, ref_samples_left, ver_src_step, ups_ver_factor); // Works for height 8, 16, 32 and 64. Upsamples 1 to 4.
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
            mip_upsampling_w32_ups4_ver_avx2(result, ver_src, ref_samples_top);
          }
          else {
            mip_upsampling_w32_ups8_ver_avx2(result, ver_src, ref_samples_top);
          }
          break;
          
        case 64: 
          if (ups_ver_factor == 2) {
            mip_upsampling_w64_ups2_ver_avx2(result, ver_src, ref_samples_top);
          }
          else if (ups_ver_factor == 4) {
            mip_upsampling_w64_ups4_ver_avx2(result, ver_src, ref_samples_top);
          }
          else {
            mip_upsampling_w64_ups8_ver_avx2(result, ver_src, ref_samples_top);
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
    success &= uvg_strategyselector_register(opaque, "pdpc_planar_dc", "avx2", 40, &uvg_pdpc_planar_dc_avx2);
    success &= uvg_strategyselector_register(opaque, "mip_predict", "avx2", 40, &mip_predict_avx2);
  }
#endif //UVG_BIT_DEPTH == 8
#endif //COMPILE_INTEL_AVX2 && defined X86_64
  return success;
}

