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

/*
* \file
*/

#include "strategies/avx2/dct-avx2.h"


#include "strategyselector.h"
#include "tables.h"

extern const int16_t uvg_g_dst_4[4][4];
extern const int16_t uvg_g_dct_4[4][4];
extern const int16_t uvg_g_dct_8[8][8];
extern const int16_t uvg_g_dct_16[16][16];
extern const int16_t uvg_g_dct_32[32][32];

extern const int16_t uvg_g_dst_4_t[4][4];
extern const int16_t uvg_g_dct_4_t[4][4];
extern const int16_t uvg_g_dct_8_t[8][8];
extern const int16_t uvg_g_dct_16_t[16][16];
extern const int16_t uvg_g_dct_32_t[32][32];

#if COMPILE_INTEL_AVX2 
#include "uvg266.h"
#if UVG_BIT_DEPTH == 8
#include <immintrin.h>
#include "strategies/avx2/dct_avx2_tables.h"
#define MAX_LOG2_TR_DYNAMIC_RANGE 15
#define TRANSFORM_MATRIX_SHIFT    6
#define INVERSE_SHIFT_1ST (TRANSFORM_MATRIX_SHIFT + 1)
#define INVERSE_SHIFT_2ND (TRANSFORM_MATRIX_SHIFT + MAX_LOG2_TR_DYNAMIC_RANGE - 1 - UVG_BIT_DEPTH)

/*
* \file
* \brief AVX2 transformations.
*/

static INLINE __m256i swap_lanes(__m256i v)
{
  return _mm256_permute4x64_epi64(v, _MM_SHUFFLE(1, 0, 3, 2));
}

static INLINE __m256i truncate_avx2(__m256i v, __m256i debias, int32_t shift)
{
  __m256i truncable = _mm256_add_epi32 (v,         debias);
  return              _mm256_srai_epi32(truncable, shift);
}


// TODO: find avx2 solution for transpose
// TODO: attempt to make a generic transpose for avx2. Needs some extra logic for different widths and heights.
// TODO: make a few solutions for exact sizes and see if some pattern emerges...
static void transpose_matrix(const int16_t* src, int16_t* dst, const int width, const int height) {
  const int sample_num = width * height;
  const int vectors = sample_num / 16;

  int16_t* d_ptr = dst;
  if (vectors == 0) {
    return;
  }
  else if (vectors == 1) {

  }
  else {
    // Reserve enough storage for max transform size 32x32
    __m256i v_16b_result[64];
    __m256i v_32b_result[64];
    __m256i v_64b_result[64];
    __m256i v_128b_result[64];

    // Handle two source vectors at a time
    for (int i = 0; i < vectors; i += 2) {
      __m256i v_src_0 = _mm256_load_si256((const __m256i*)src);
      __m256i v_src_1 = _mm256_load_si256((const __m256i*)(src + 16));

      v_16b_result[i] = _mm256_unpacklo_epi16(v_src_0, v_src_1);
      v_16b_result[i + 1] = _mm256_unpackhi_epi16(v_src_0, v_src_1);

      src += 32;
    }

    // 32 bit shuffle pass
    int loop_idx = 0;
    for (int i = 0; i < vectors; i += 2) {
      const int idx_a = loop_idx;
      const int idx_b = loop_idx + 2;

      v_32b_result[i] = _mm256_unpacklo_epi32(v_16b_result[idx_a], v_16b_result[idx_b]);
      v_32b_result[i + 1] = _mm256_unpackhi_epi32(v_16b_result[idx_a], v_16b_result[idx_b]);
      loop_idx++;
    }

    // 64 bit shuffle pass
    loop_idx = 0;
    for (int i = 0; i < vectors; i += 2) {
      const int idx_a = loop_idx;
      const int idx_b = loop_idx + 4;

      v_64b_result[i] = _mm256_unpacklo_epi32(v_32b_result[idx_a], v_32b_result[idx_b]);
      v_64b_result[i + 1] = _mm256_unpackhi_epi32(v_32b_result[idx_a], v_32b_result[idx_b]);
      loop_idx++;
    }

    // Final 128 bit shuffle pass
    for (int i = 0; i < vectors; i += 2) {
      const int idx_a = 0;
      const int idx_b = 0;

      v_128b_result[i] = _mm256_unpacklo_epi32(v_64b_result[idx_a], v_64b_result[idx_b]);
      v_128b_result[i + 1] = _mm256_unpackhi_epi32(v_64b_result[idx_a], v_64b_result[idx_b]);
    }

    // Store loop
    for (int i = 0; i < vectors; ++i) {
      _mm256_store_si256((__m256i*)dst, v_128b_result[i]);
      dst += 16;
    }
  }
}

static void transpose_generic(const int16_t* src, int16_t* dst, const int width, const int height)
{
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      dst[x * height + y] = src[y * width + x];
    }
  }
}


typedef void (transpose_func)(const __m256i* src, __m256i* dst);


static void transpose_2x2_avx2(const __m256i* src, __m256i* dst){}
static void transpose_2x4_avx2(const __m256i* src, __m256i* dst){}
static void transpose_2x8_avx2(const __m256i* src, __m256i* dst){}
static void transpose_2x16_avx2(const __m256i* src, __m256i* dst){}
static void transpose_2x32_avx2(const __m256i* src, __m256i* dst)
{
  const __m256i v_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0246);
  __m256i v_tmp[4];
  v_tmp[0] = _mm256_shuffle_epi8(src[0], v_shuffle);
  v_tmp[1] = _mm256_shuffle_epi8(src[1], v_shuffle);
  v_tmp[2] = _mm256_shuffle_epi8(src[2], v_shuffle);
  v_tmp[3] = _mm256_shuffle_epi8(src[3], v_shuffle);

  v_tmp[0] = _mm256_permute4x64_epi64(v_tmp[0], _MM_SHUFFLE(3, 1, 2, 0));
  v_tmp[1] = _mm256_permute4x64_epi64(v_tmp[1], _MM_SHUFFLE(3, 1, 2, 0));
  v_tmp[2] = _mm256_permute4x64_epi64(v_tmp[2], _MM_SHUFFLE(3, 1, 2, 0));
  v_tmp[3] = _mm256_permute4x64_epi64(v_tmp[3], _MM_SHUFFLE(3, 1, 2, 0));

  dst[0] = _mm256_permute2x128_si256(v_tmp[0], v_tmp[1], 0x20);
  dst[1] = _mm256_permute2x128_si256(v_tmp[2], v_tmp[3], 0x20);
  dst[2] = _mm256_permute2x128_si256(v_tmp[0], v_tmp[1], 0x31);
  dst[3] = _mm256_permute2x128_si256(v_tmp[2], v_tmp[3], 0x31);
}
static void transpose_2x64_avx2(const __m256i* src, __m256i* dst){}
static void transpose_4x2_avx2(const __m256i* src, __m256i* dst){}
static void transpose_4x4_avx2(const __m256i* src, __m256i* dst){}
static void transpose_4x8_avx2(const __m256i* src, __m256i* dst){}
static void transpose_4x16_avx2(const __m256i* src, __m256i* dst)
{
  const __m256i v_shuffle = _mm256_set_epi8(15, 14,  7,  6, 13, 12,  5,  4, 11, 10,  3,  2,  9,  8,  1,  0,
                                            31, 30, 23, 22, 29, 28, 21, 20, 27, 26, 19, 18, 25, 24, 17, 16);

  // const __m256i v_shuffle = _mm256_set_epi8( 0, 1, 8, 9, 2, 3, 10, 11, 4, 5, 12, 13, 6, 7, 14, 15,
  //                                           16, 17, 24, 25, 18, 19, 26, 27, 20, 21, 28, 29, 22, 23, 30, 31);

  __m256i v_src_tmp[4];
  v_src_tmp[0] = _mm256_shuffle_epi8(src[0], v_shuffle);
  v_src_tmp[1] = _mm256_shuffle_epi8(src[1], v_shuffle);
  v_src_tmp[2] = _mm256_shuffle_epi8(src[2], v_shuffle);
  v_src_tmp[3] = _mm256_shuffle_epi8(src[3], v_shuffle);

  __m256i v_tmp[4];
  v_tmp[0] = _mm256_permute2x128_si256(v_src_tmp[0], v_src_tmp[1], 0x20);
  v_tmp[1] = _mm256_permute2x128_si256(v_src_tmp[0], v_src_tmp[1], 0x31);
  v_tmp[2] = _mm256_permute2x128_si256(v_src_tmp[2], v_src_tmp[3], 0x20);
  v_tmp[3] = _mm256_permute2x128_si256(v_src_tmp[2], v_src_tmp[3], 0x31);

  __m256i v_tmp16_lo[2];
  __m256i v_tmp16_hi[2];
  v_tmp16_lo[0] = _mm256_unpacklo_epi32(v_tmp[0], v_tmp[1]);
  v_tmp16_lo[1] = _mm256_unpacklo_epi32(v_tmp[2], v_tmp[3]);
  v_tmp16_hi[0] = _mm256_unpackhi_epi32(v_tmp[0], v_tmp[1]);
  v_tmp16_hi[1] = _mm256_unpackhi_epi32(v_tmp[2], v_tmp[3]);

  v_tmp[0] = _mm256_permute4x64_epi64(v_tmp16_lo[0], _MM_SHUFFLE(3, 1, 2, 0));
  v_tmp[1] = _mm256_permute4x64_epi64(v_tmp16_lo[1], _MM_SHUFFLE(3, 1, 2, 0));
  v_tmp[2] = _mm256_permute4x64_epi64(v_tmp16_hi[0], _MM_SHUFFLE(3, 1, 2, 0));
  v_tmp[3] = _mm256_permute4x64_epi64(v_tmp16_hi[1], _MM_SHUFFLE(3, 1, 2, 0));

  dst[0] = _mm256_permute2x128_si256(v_tmp[0], v_tmp[1], 0x20);
  dst[1] = _mm256_permute2x128_si256(v_tmp[0], v_tmp[1], 0x31);
  dst[2] = _mm256_permute2x128_si256(v_tmp[2], v_tmp[3], 0x20);
  dst[3] = _mm256_permute2x128_si256(v_tmp[2], v_tmp[3], 0x31);
}
static void transpose_4x32_avx2(const __m256i* src, __m256i* dst)
{
  __m256i v_tmp[8];
  const __m256i v_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0415);
  for (int i = 0; i < 8; ++i) {
    v_tmp[i] = _mm256_shuffle_epi8(src[i], v_shuffle);
    v_tmp[i] = _mm256_permute4x64_epi64(v_tmp[i], _MM_SHUFFLE(3, 1, 2, 0));
    v_tmp[i] = _mm256_shuffle_epi32(v_tmp[i], _MM_SHUFFLE(3, 1, 2, 0));
  }

  __m256i v_tmp64_lo[4];
  __m256i v_tmp64_hi[4];
  v_tmp64_lo[0] = _mm256_unpacklo_epi64(v_tmp[0], v_tmp[1]);
  v_tmp64_lo[1] = _mm256_unpacklo_epi64(v_tmp[2], v_tmp[3]);
  v_tmp64_lo[2] = _mm256_unpacklo_epi64(v_tmp[4], v_tmp[5]);
  v_tmp64_lo[3] = _mm256_unpacklo_epi64(v_tmp[6], v_tmp[7]);

  v_tmp64_hi[0] = _mm256_unpackhi_epi64(v_tmp[0], v_tmp[1]);
  v_tmp64_hi[1] = _mm256_unpackhi_epi64(v_tmp[2], v_tmp[3]);
  v_tmp64_hi[2] = _mm256_unpackhi_epi64(v_tmp[4], v_tmp[5]);
  v_tmp64_hi[3] = _mm256_unpackhi_epi64(v_tmp[6], v_tmp[7]);

  dst[0] = _mm256_permute2x128_si256(v_tmp64_lo[0], v_tmp64_lo[1], 0x20);
  dst[1] = _mm256_permute2x128_si256(v_tmp64_lo[2], v_tmp64_lo[3], 0x20);
  dst[2] = _mm256_permute2x128_si256(v_tmp64_hi[0], v_tmp64_hi[1], 0x20);
  dst[3] = _mm256_permute2x128_si256(v_tmp64_hi[2], v_tmp64_hi[3], 0x20);

  dst[4] = _mm256_permute2x128_si256(v_tmp64_lo[0], v_tmp64_lo[1], 0x31);
  dst[5] = _mm256_permute2x128_si256(v_tmp64_lo[2], v_tmp64_lo[3], 0x31);
  dst[6] = _mm256_permute2x128_si256(v_tmp64_hi[0], v_tmp64_hi[1], 0x31);
  dst[7] = _mm256_permute2x128_si256(v_tmp64_hi[2], v_tmp64_hi[3], 0x31);
}
static void transpose_4x64_avx2(const __m256i* src, __m256i* dst){}
static void transpose_8x2_avx2(const __m256i* src, __m256i* dst){}
static void transpose_8x4_avx2(const __m256i* src, __m256i* dst){}
static void transpose_8x8_avx2(const __m256i* src, __m256i* dst){}
static void transpose_8x16_avx2(const __m256i* src, __m256i* dst)
{
  __m256i v_tmp16_lo[4];
  __m256i v_tmp16_hi[4];
  __m256i v_tmp32_lo[4];
  __m256i v_tmp32_hi[4];
  __m256i v_tmp64_lo[4];
  __m256i v_tmp64_hi[4];
  __m256i v_tmp128[8];

  v_tmp128[0] = _mm256_permute2x128_si256(src[0], src[4], 0x20);
  v_tmp128[1] = _mm256_permute2x128_si256(src[0], src[4], 0x31);
  v_tmp128[2] = _mm256_permute2x128_si256(src[1], src[5], 0x20);
  v_tmp128[3] = _mm256_permute2x128_si256(src[1], src[5], 0x31);
  v_tmp128[4] = _mm256_permute2x128_si256(src[2], src[6], 0x20);
  v_tmp128[5] = _mm256_permute2x128_si256(src[2], src[6], 0x31);
  v_tmp128[6] = _mm256_permute2x128_si256(src[3], src[7], 0x20);
  v_tmp128[7] = _mm256_permute2x128_si256(src[3], src[7], 0x31);

  v_tmp16_lo[0] = _mm256_unpacklo_epi16(v_tmp128[0], v_tmp128[1]);
  v_tmp16_lo[1] = _mm256_unpacklo_epi16(v_tmp128[2], v_tmp128[3]);
  v_tmp16_lo[2] = _mm256_unpacklo_epi16(v_tmp128[4], v_tmp128[5]);
  v_tmp16_lo[3] = _mm256_unpacklo_epi16(v_tmp128[6], v_tmp128[7]);
  v_tmp16_hi[0] = _mm256_unpackhi_epi16(v_tmp128[0], v_tmp128[1]);
  v_tmp16_hi[1] = _mm256_unpackhi_epi16(v_tmp128[2], v_tmp128[3]);
  v_tmp16_hi[2] = _mm256_unpackhi_epi16(v_tmp128[4], v_tmp128[5]);
  v_tmp16_hi[3] = _mm256_unpackhi_epi16(v_tmp128[6], v_tmp128[7]);

  v_tmp32_lo[0] = _mm256_unpacklo_epi32(v_tmp16_lo[0], v_tmp16_lo[1]);
  v_tmp32_lo[1] = _mm256_unpacklo_epi32(v_tmp16_lo[2], v_tmp16_lo[3]);
  v_tmp32_lo[2] = _mm256_unpacklo_epi32(v_tmp16_hi[0], v_tmp16_hi[1]);
  v_tmp32_lo[3] = _mm256_unpacklo_epi32(v_tmp16_hi[2], v_tmp16_hi[3]);
  v_tmp32_hi[0] = _mm256_unpackhi_epi32(v_tmp16_lo[0], v_tmp16_lo[1]);
  v_tmp32_hi[1] = _mm256_unpackhi_epi32(v_tmp16_lo[2], v_tmp16_lo[3]);
  v_tmp32_hi[2] = _mm256_unpackhi_epi32(v_tmp16_hi[0], v_tmp16_hi[1]);
  v_tmp32_hi[3] = _mm256_unpackhi_epi32(v_tmp16_hi[2], v_tmp16_hi[3]);

  dst[0] = _mm256_unpacklo_epi64(v_tmp32_lo[0], v_tmp32_lo[1]);
  dst[1] = _mm256_unpackhi_epi64(v_tmp32_lo[0], v_tmp32_lo[1]);
  dst[2] = _mm256_unpacklo_epi64(v_tmp32_hi[0], v_tmp32_hi[1]);
  dst[3] = _mm256_unpackhi_epi64(v_tmp32_hi[0], v_tmp32_hi[1]);
  dst[4] = _mm256_unpacklo_epi64(v_tmp32_lo[2], v_tmp32_lo[3]);
  dst[5] = _mm256_unpackhi_epi64(v_tmp32_lo[2], v_tmp32_lo[3]);
  dst[6] = _mm256_unpacklo_epi64(v_tmp32_hi[2], v_tmp32_hi[3]);
  dst[7] = _mm256_unpackhi_epi64(v_tmp32_hi[2], v_tmp32_hi[3]);
}
static void transpose_8x32_avx2(const __m256i* src, __m256i* dst)
{
  __m256i v_tmp16_lo[8];
  __m256i v_tmp16_hi[8];
  __m256i v_tmp32_lo[8];
  __m256i v_tmp32_hi[8];
  __m256i v_tmp64_lo[8];
  __m256i v_tmp64_hi[8];

  const __m256i v_shuffle = _mm256_setr_epi8(0, 1, 8, 9, 2, 3, 10, 11, 4, 5, 12, 13, 6, 7, 14, 15,
    16, 17, 24, 25, 18, 19, 26, 27, 20, 21, 28, 29, 22, 23, 30, 31);
  for (int i = 0; i < 8; ++i) {
    const int offset = i * 2;
    v_tmp16_lo[i] = _mm256_unpacklo_epi16(src[offset], src[offset + 1]);
    v_tmp16_hi[i] = _mm256_unpackhi_epi16(src[offset], src[offset + 1]);
  }

  for (int i = 0; i < 8; i += 4) {
    v_tmp32_lo[i + 0] = _mm256_unpacklo_epi32(v_tmp16_lo[i + 0], v_tmp16_lo[i + 1]);
    v_tmp32_lo[i + 1] = _mm256_unpacklo_epi32(v_tmp16_lo[i + 2], v_tmp16_lo[i + 3]);
    v_tmp32_lo[i + 2] = _mm256_unpacklo_epi32(v_tmp16_hi[i + 0], v_tmp16_hi[i + 1]);
    v_tmp32_lo[i + 3] = _mm256_unpacklo_epi32(v_tmp16_hi[i + 2], v_tmp16_hi[i + 3]);

    v_tmp32_hi[i + 0] = _mm256_unpackhi_epi32(v_tmp16_lo[i + 0], v_tmp16_lo[i + 1]);
    v_tmp32_hi[i + 1] = _mm256_unpackhi_epi32(v_tmp16_lo[i + 2], v_tmp16_lo[i + 3]);
    v_tmp32_hi[i + 2] = _mm256_unpackhi_epi32(v_tmp16_hi[i + 0], v_tmp16_hi[i + 1]);
    v_tmp32_hi[i + 3] = _mm256_unpackhi_epi32(v_tmp16_hi[i + 2], v_tmp16_hi[i + 3]);
  }

  for (int i = 0; i < 8; i += 4) {
    v_tmp64_lo[i + 0] = _mm256_unpacklo_epi64(v_tmp32_lo[i + 0], v_tmp32_lo[i + 1]);
    v_tmp64_lo[i + 1] = _mm256_unpacklo_epi64(v_tmp32_lo[i + 2], v_tmp32_lo[i + 3]);
    v_tmp64_lo[i + 2] = _mm256_unpacklo_epi64(v_tmp32_hi[i + 0], v_tmp32_hi[i + 1]);
    v_tmp64_lo[i + 3] = _mm256_unpacklo_epi64(v_tmp32_hi[i + 2], v_tmp32_hi[i + 3]);

    v_tmp64_hi[i + 0] = _mm256_unpackhi_epi64(v_tmp32_lo[i + 0], v_tmp32_lo[i + 1]);
    v_tmp64_hi[i + 1] = _mm256_unpackhi_epi64(v_tmp32_lo[i + 2], v_tmp32_lo[i + 3]);
    v_tmp64_hi[i + 2] = _mm256_unpackhi_epi64(v_tmp32_hi[i + 0], v_tmp32_hi[i + 1]);
    v_tmp64_hi[i + 3] = _mm256_unpackhi_epi64(v_tmp32_hi[i + 2], v_tmp32_hi[i + 3]);
  }

  for (int i = 0; i < 8; ++i) {
    v_tmp64_lo[i] = _mm256_permute4x64_epi64(v_tmp64_lo[i], _MM_SHUFFLE(3, 1, 2, 0));
    v_tmp64_hi[i] = _mm256_permute4x64_epi64(v_tmp64_hi[i], _MM_SHUFFLE(3, 1, 2, 0));
  }

  dst[0] = _mm256_shuffle_epi8(v_tmp64_lo[0], v_shuffle);
  dst[1] = _mm256_shuffle_epi8(v_tmp64_lo[4], v_shuffle);
  dst[2] = _mm256_shuffle_epi8(v_tmp64_hi[0], v_shuffle);
  dst[3] = _mm256_shuffle_epi8(v_tmp64_hi[4], v_shuffle);

  dst[4] = _mm256_shuffle_epi8(v_tmp64_lo[2], v_shuffle);
  dst[5] = _mm256_shuffle_epi8(v_tmp64_lo[6], v_shuffle);
  dst[6] = _mm256_shuffle_epi8(v_tmp64_hi[2], v_shuffle);
  dst[7] = _mm256_shuffle_epi8(v_tmp64_hi[6], v_shuffle);

  dst[8] = _mm256_shuffle_epi8(v_tmp64_lo[1], v_shuffle);
  dst[9] = _mm256_shuffle_epi8(v_tmp64_lo[5], v_shuffle);
  dst[10] = _mm256_shuffle_epi8(v_tmp64_hi[1], v_shuffle);
  dst[11] = _mm256_shuffle_epi8(v_tmp64_hi[5], v_shuffle);

  dst[12] = _mm256_shuffle_epi8(v_tmp64_lo[3], v_shuffle);
  dst[13] = _mm256_shuffle_epi8(v_tmp64_lo[7], v_shuffle);
  dst[14] = _mm256_shuffle_epi8(v_tmp64_hi[3], v_shuffle);
  dst[15] = _mm256_shuffle_epi8(v_tmp64_hi[7], v_shuffle);
}
static void transpose_8x64_avx2(const __m256i* src, __m256i* dst){}
static void transpose_16x2_avx2(const __m256i* src, __m256i* dst){}
static void transpose_16x4_avx2(const __m256i* src, __m256i* dst)
{
  __m256i v_tmp16_lo[2];
  __m256i v_tmp16_hi[2];
  __m256i v_tmp32_lo[2];
  __m256i v_tmp32_hi[2];

  v_tmp16_lo[0] = _mm256_unpacklo_epi16(src[0], src[1]);
  v_tmp16_lo[1] = _mm256_unpacklo_epi16(src[2], src[3]);
  v_tmp16_hi[0] = _mm256_unpackhi_epi16(src[0], src[1]);
  v_tmp16_hi[1] = _mm256_unpackhi_epi16(src[2], src[3]);

  v_tmp32_lo[0] = _mm256_unpacklo_epi32(v_tmp16_lo[0], v_tmp16_lo[1]);
  v_tmp32_lo[1] = _mm256_unpacklo_epi32(v_tmp16_hi[0], v_tmp16_hi[1]);

  v_tmp32_hi[0] = _mm256_unpackhi_epi32(v_tmp16_lo[0], v_tmp16_lo[1]);
  v_tmp32_hi[1] = _mm256_unpackhi_epi32(v_tmp16_hi[0], v_tmp16_hi[1]);

  dst[0] = _mm256_permute2x128_si256(v_tmp32_lo[0], v_tmp32_hi[0], 0x20);
  dst[1] = _mm256_permute2x128_si256(v_tmp32_lo[1], v_tmp32_hi[1], 0x20);
  dst[2] = _mm256_permute2x128_si256(v_tmp32_lo[0], v_tmp32_hi[0], 0x31);
  dst[3] = _mm256_permute2x128_si256(v_tmp32_lo[1], v_tmp32_hi[1], 0x31);
}
static void transpose_16x8_avx2(const __m256i* src, __m256i* dst)
{
  __m256i v_tmp16_lo[4];
  __m256i v_tmp16_hi[4];
  __m256i v_tmp32_lo[4];
  __m256i v_tmp32_hi[4];
  __m256i v_tmp64_lo[4];
  __m256i v_tmp64_hi[4];
  v_tmp16_lo[0] = _mm256_unpacklo_epi16(src[0], src[1]);
  v_tmp16_lo[1] = _mm256_unpacklo_epi16(src[2], src[3]);
  v_tmp16_lo[2] = _mm256_unpacklo_epi16(src[4], src[5]);
  v_tmp16_lo[3] = _mm256_unpacklo_epi16(src[6], src[7]);
  v_tmp16_hi[0] = _mm256_unpackhi_epi16(src[0], src[1]);
  v_tmp16_hi[1] = _mm256_unpackhi_epi16(src[2], src[3]);
  v_tmp16_hi[2] = _mm256_unpackhi_epi16(src[4], src[5]);
  v_tmp16_hi[3] = _mm256_unpackhi_epi16(src[6], src[7]);

  v_tmp32_lo[0] = _mm256_unpacklo_epi32(v_tmp16_lo[0], v_tmp16_lo[1]);
  v_tmp32_lo[1] = _mm256_unpacklo_epi32(v_tmp16_lo[2], v_tmp16_lo[3]);
  v_tmp32_lo[2] = _mm256_unpacklo_epi32(v_tmp16_hi[0], v_tmp16_hi[1]);
  v_tmp32_lo[3] = _mm256_unpacklo_epi32(v_tmp16_hi[2], v_tmp16_hi[3]);
  v_tmp32_hi[0] = _mm256_unpackhi_epi32(v_tmp16_lo[0], v_tmp16_lo[1]);
  v_tmp32_hi[1] = _mm256_unpackhi_epi32(v_tmp16_lo[2], v_tmp16_lo[3]);
  v_tmp32_hi[2] = _mm256_unpackhi_epi32(v_tmp16_hi[0], v_tmp16_hi[1]);
  v_tmp32_hi[3] = _mm256_unpackhi_epi32(v_tmp16_hi[2], v_tmp16_hi[3]);

  v_tmp64_lo[0] = _mm256_unpacklo_epi64(v_tmp32_lo[0], v_tmp32_lo[1]);
  v_tmp64_lo[1] = _mm256_unpacklo_epi64(v_tmp32_lo[2], v_tmp32_lo[3]);
  v_tmp64_lo[2] = _mm256_unpacklo_epi64(v_tmp32_hi[0], v_tmp32_hi[1]);
  v_tmp64_lo[3] = _mm256_unpacklo_epi64(v_tmp32_hi[2], v_tmp32_hi[3]);
  v_tmp64_hi[0] = _mm256_unpackhi_epi64(v_tmp32_lo[0], v_tmp32_lo[1]);
  v_tmp64_hi[1] = _mm256_unpackhi_epi64(v_tmp32_lo[2], v_tmp32_lo[3]);
  v_tmp64_hi[2] = _mm256_unpackhi_epi64(v_tmp32_hi[0], v_tmp32_hi[1]);
  v_tmp64_hi[3] = _mm256_unpackhi_epi64(v_tmp32_hi[2], v_tmp32_hi[3]);

  dst[0] = _mm256_permute2x128_si256(v_tmp64_lo[0], v_tmp64_hi[0], 0x20);
  dst[1] = _mm256_permute2x128_si256(v_tmp64_lo[2], v_tmp64_hi[2], 0x20);
  dst[2] = _mm256_permute2x128_si256(v_tmp64_lo[1], v_tmp64_hi[1], 0x20);
  dst[3] = _mm256_permute2x128_si256(v_tmp64_lo[3], v_tmp64_hi[3], 0x20);
  dst[4] = _mm256_permute2x128_si256(v_tmp64_lo[0], v_tmp64_hi[0], 0x31);
  dst[5] = _mm256_permute2x128_si256(v_tmp64_lo[2], v_tmp64_hi[2], 0x31);
  dst[6] = _mm256_permute2x128_si256(v_tmp64_lo[1], v_tmp64_hi[1], 0x31);
  dst[7] = _mm256_permute2x128_si256(v_tmp64_lo[3], v_tmp64_hi[3], 0x31);
}

static void transpose_16x16_avx2_stride(const int16_t* src, int16_t* dst, const int src_stride, const int dst_stride) {
  __m256i v_tmp16_lo[8];
  __m256i v_tmp16_hi[8];
  for (int d = 0, s = 0; d < 8; ++d, s += 2) {
    v_tmp16_lo[d] = _mm256_unpacklo_epi16(*(__m256i*)(src + s * src_stride), *(__m256i*)(src + (s + 1) * src_stride));
    v_tmp16_hi[d] = _mm256_unpackhi_epi16(*(__m256i*)(src + s * src_stride), *(__m256i*)(src + (s + 1) * src_stride));
  }

  __m256i v_tmp32_lo[8];
  __m256i v_tmp32_hi[8];
  for (int d = 0, s = 0; d < 8; d += 2, s += 2) {
    v_tmp32_lo[d + 0] = _mm256_unpacklo_epi32(v_tmp16_lo[s + 0], v_tmp16_lo[s + 1]);
    v_tmp32_lo[d + 1] = _mm256_unpacklo_epi32(v_tmp16_hi[s + 0], v_tmp16_hi[s + 1]);
    v_tmp32_hi[d + 0] = _mm256_unpackhi_epi32(v_tmp16_lo[s + 0], v_tmp16_lo[s + 1]);
    v_tmp32_hi[d + 1] = _mm256_unpackhi_epi32(v_tmp16_hi[s + 0], v_tmp16_hi[s + 1]);
  }

  __m256i v_tmp64_lo[8];
  __m256i v_tmp64_hi[8];
  for (int d = 0, s = 0; d < 8; d += 4, s += 4) {
    v_tmp64_lo[d + 0] = _mm256_unpacklo_epi64(v_tmp32_lo[s + 0], v_tmp32_lo[s + 2]);
    v_tmp64_lo[d + 1] = _mm256_unpacklo_epi64(v_tmp32_lo[s + 1], v_tmp32_lo[s + 3]);
    v_tmp64_hi[d + 0] = _mm256_unpackhi_epi64(v_tmp32_lo[s + 0], v_tmp32_lo[s + 2]);
    v_tmp64_hi[d + 1] = _mm256_unpackhi_epi64(v_tmp32_lo[s + 1], v_tmp32_lo[s + 3]);

    v_tmp64_lo[d + 2] = _mm256_unpacklo_epi64(v_tmp32_hi[s + 0], v_tmp32_hi[s + 2]);
    v_tmp64_lo[d + 3] = _mm256_unpacklo_epi64(v_tmp32_hi[s + 1], v_tmp32_hi[s + 3]);
    v_tmp64_hi[d + 2] = _mm256_unpackhi_epi64(v_tmp32_hi[s + 0], v_tmp32_hi[s + 2]);
    v_tmp64_hi[d + 3] = _mm256_unpackhi_epi64(v_tmp32_hi[s + 1], v_tmp32_hi[s + 3]);
  }

  _mm256_storeu_si256((__m256i*)(dst + 0 * dst_stride), _mm256_permute2x128_si256(v_tmp64_lo[0], v_tmp64_lo[4], 0x20));
  _mm256_storeu_si256((__m256i*)(dst + 1 * dst_stride), _mm256_permute2x128_si256(v_tmp64_hi[0], v_tmp64_hi[4], 0x20));
  _mm256_storeu_si256((__m256i*)(dst + 2 * dst_stride), _mm256_permute2x128_si256(v_tmp64_lo[2], v_tmp64_lo[6], 0x20));
  _mm256_storeu_si256((__m256i*)(dst + 3 * dst_stride), _mm256_permute2x128_si256(v_tmp64_hi[2], v_tmp64_hi[6], 0x20));
  _mm256_storeu_si256((__m256i*)(dst + 4 * dst_stride), _mm256_permute2x128_si256(v_tmp64_lo[1], v_tmp64_lo[5], 0x20));
  _mm256_storeu_si256((__m256i*)(dst + 5 * dst_stride), _mm256_permute2x128_si256(v_tmp64_hi[1], v_tmp64_hi[5], 0x20));
  _mm256_storeu_si256((__m256i*)(dst + 6 * dst_stride), _mm256_permute2x128_si256(v_tmp64_lo[3], v_tmp64_lo[7], 0x20));
  _mm256_storeu_si256((__m256i*)(dst + 7 * dst_stride), _mm256_permute2x128_si256(v_tmp64_hi[3], v_tmp64_hi[7], 0x20));

  _mm256_storeu_si256((__m256i*)(dst + 8 * dst_stride), _mm256_permute2x128_si256(v_tmp64_lo[0], v_tmp64_lo[4], 0x31));
  _mm256_storeu_si256((__m256i*)(dst + 9 * dst_stride), _mm256_permute2x128_si256(v_tmp64_hi[0], v_tmp64_hi[4], 0x31));
  _mm256_storeu_si256((__m256i*)(dst + 10 * dst_stride), _mm256_permute2x128_si256(v_tmp64_lo[2], v_tmp64_lo[6], 0x31));
  _mm256_storeu_si256((__m256i*)(dst + 11 * dst_stride), _mm256_permute2x128_si256(v_tmp64_hi[2], v_tmp64_hi[6], 0x31));
  _mm256_storeu_si256((__m256i*)(dst + 12 * dst_stride), _mm256_permute2x128_si256(v_tmp64_lo[1], v_tmp64_lo[5], 0x31));
  _mm256_storeu_si256((__m256i*)(dst + 13 * dst_stride), _mm256_permute2x128_si256(v_tmp64_hi[1], v_tmp64_hi[5], 0x31));
  _mm256_storeu_si256((__m256i*)(dst + 14 * dst_stride), _mm256_permute2x128_si256(v_tmp64_lo[3], v_tmp64_lo[7], 0x31));
  _mm256_storeu_si256((__m256i*)(dst + 15 * dst_stride), _mm256_permute2x128_si256(v_tmp64_hi[3], v_tmp64_hi[7], 0x31));
}

static void transpose_16x16_avx2(const __m256i* src, __m256i* dst) {
  transpose_16x16_avx2_stride((int16_t const *)src, (int16_t*)dst, 16, 16);
}

static void transpose_16x32_avx2(const __m256i* src, __m256i* dst) {
  transpose_16x16_avx2_stride((int16_t const *)src, (int16_t*)dst, 16, 32);
  transpose_16x16_avx2_stride((int16_t const *)src + 16 * 16, (int16_t*)dst + 16, 16, 32);

}
static void transpose_16x64_avx2(const __m256i* src, __m256i* dst){}
static void transpose_32x2_avx2(const __m256i* src, __m256i* dst)
{
  __m256i v_tmp16_lo0 = _mm256_unpacklo_epi16(src[0], src[2]);
  __m256i v_tmp16_lo1 = _mm256_unpacklo_epi16(src[1], src[3]);
  __m256i v_tmp16_hi0 = _mm256_unpackhi_epi16(src[0], src[2]);
  __m256i v_tmp16_hi1 = _mm256_unpackhi_epi16(src[1], src[3]);

  dst[0] = _mm256_permute2x128_si256(v_tmp16_lo0, v_tmp16_hi0, 0x20);
  dst[1] = _mm256_permute2x128_si256(v_tmp16_lo0, v_tmp16_hi0, 0x31);
  dst[2] = _mm256_permute2x128_si256(v_tmp16_lo1, v_tmp16_hi1, 0x20);
  dst[3] = _mm256_permute2x128_si256(v_tmp16_lo1, v_tmp16_hi1, 0x31);
}
static void transpose_32x4_avx2(const __m256i* src, __m256i* dst)
{
  __m256i v_tmp16_lo[4];
  __m256i v_tmp16_hi[4];
  v_tmp16_lo[0] = _mm256_unpacklo_epi16(src[0], src[2]);
  v_tmp16_lo[1] = _mm256_unpacklo_epi16(src[1], src[3]);
  v_tmp16_lo[2] = _mm256_unpacklo_epi16(src[4], src[6]);
  v_tmp16_lo[3] = _mm256_unpacklo_epi16(src[5], src[7]);

  v_tmp16_hi[0] = _mm256_unpackhi_epi16(src[0], src[2]);
  v_tmp16_hi[1] = _mm256_unpackhi_epi16(src[1], src[3]);
  v_tmp16_hi[2] = _mm256_unpackhi_epi16(src[4], src[6]);
  v_tmp16_hi[3] = _mm256_unpackhi_epi16(src[5], src[7]);

  __m256i v_tmp32_lo[4];
  __m256i v_tmp32_hi[4];
  v_tmp32_lo[0] = _mm256_unpacklo_epi32(v_tmp16_lo[0], v_tmp16_lo[2]);
  v_tmp32_lo[1] = _mm256_unpacklo_epi32(v_tmp16_lo[1], v_tmp16_lo[3]);
  v_tmp32_lo[2] = _mm256_unpacklo_epi32(v_tmp16_hi[0], v_tmp16_hi[2]);
  v_tmp32_lo[3] = _mm256_unpacklo_epi32(v_tmp16_hi[1], v_tmp16_hi[3]);

  v_tmp32_hi[0] = _mm256_unpackhi_epi32(v_tmp16_lo[0], v_tmp16_lo[2]);
  v_tmp32_hi[1] = _mm256_unpackhi_epi32(v_tmp16_lo[1], v_tmp16_lo[3]);
  v_tmp32_hi[2] = _mm256_unpackhi_epi32(v_tmp16_hi[0], v_tmp16_hi[2]);
  v_tmp32_hi[3] = _mm256_unpackhi_epi32(v_tmp16_hi[1], v_tmp16_hi[3]);

  dst[0] = _mm256_permute2x128_si256(v_tmp32_lo[0], v_tmp32_hi[0], 0x20);
  dst[1] = _mm256_permute2x128_si256(v_tmp32_lo[2], v_tmp32_hi[2], 0x20);
  dst[2] = _mm256_permute2x128_si256(v_tmp32_lo[0], v_tmp32_hi[0], 0x31);
  dst[3] = _mm256_permute2x128_si256(v_tmp32_lo[2], v_tmp32_hi[2], 0x31);

  dst[4] = _mm256_permute2x128_si256(v_tmp32_lo[1], v_tmp32_hi[1], 0x20);
  dst[5] = _mm256_permute2x128_si256(v_tmp32_lo[3], v_tmp32_hi[3], 0x20);
  dst[6] = _mm256_permute2x128_si256(v_tmp32_lo[1], v_tmp32_hi[1], 0x31);
  dst[7] = _mm256_permute2x128_si256(v_tmp32_lo[3], v_tmp32_hi[3], 0x31);
}
static void transpose_32x8_avx2(const __m256i* src, __m256i* dst)
{
  __m256i v_tmp16_lo[8];
  __m256i v_tmp16_hi[8];
  for (int d = 0, s = 0; d < 8; d += 2, s += 4) {
    v_tmp16_lo[d + 0] = _mm256_unpacklo_epi16(src[s + 0], src[s + 2]);
    v_tmp16_lo[d + 1] = _mm256_unpacklo_epi16(src[s + 1], src[s + 3]);

    v_tmp16_hi[d + 0] = _mm256_unpackhi_epi16(src[s + 0], src[s + 2]);
    v_tmp16_hi[d + 1] = _mm256_unpackhi_epi16(src[s + 1], src[s + 3]);
  }

  __m256i v_tmp32_lo[8];
  __m256i v_tmp32_hi[8];
  for (int d = 0, s = 0; d < 4; d += 2, s += 4) {
    v_tmp32_lo[d + 0] = _mm256_unpacklo_epi32(v_tmp16_lo[s + 0], v_tmp16_lo[s + 2]);
    v_tmp32_lo[d + 1] = _mm256_unpacklo_epi32(v_tmp16_lo[s + 1], v_tmp16_lo[s + 3]);
    v_tmp32_lo[d + 4] = _mm256_unpacklo_epi32(v_tmp16_hi[s + 0], v_tmp16_hi[s + 2]);
    v_tmp32_lo[d + 5] = _mm256_unpacklo_epi32(v_tmp16_hi[s + 1], v_tmp16_hi[s + 3]);

    v_tmp32_hi[d + 0] = _mm256_unpackhi_epi32(v_tmp16_lo[s + 0], v_tmp16_lo[s + 2]);
    v_tmp32_hi[d + 1] = _mm256_unpackhi_epi32(v_tmp16_lo[s + 1], v_tmp16_lo[s + 3]);
    v_tmp32_hi[d + 4] = _mm256_unpackhi_epi32(v_tmp16_hi[s + 0], v_tmp16_hi[s + 2]);
    v_tmp32_hi[d + 5] = _mm256_unpackhi_epi32(v_tmp16_hi[s + 1], v_tmp16_hi[s + 3]);
  }

  __m256i v_tmp64_lo[8];
  __m256i v_tmp64_hi[8];
  for (int d = 0, s = 0; d < 4; d += 2, s += 4) {
    v_tmp64_lo[d + 0] = _mm256_unpacklo_epi64(v_tmp32_lo[s + 0], v_tmp32_lo[s + 2]);
    v_tmp64_lo[d + 1] = _mm256_unpacklo_epi64(v_tmp32_lo[s + 1], v_tmp32_lo[s + 3]);
    v_tmp64_lo[d + 4] = _mm256_unpacklo_epi64(v_tmp32_hi[s + 0], v_tmp32_hi[s + 2]);
    v_tmp64_lo[d + 5] = _mm256_unpacklo_epi64(v_tmp32_hi[s + 1], v_tmp32_hi[s + 3]);

    v_tmp64_hi[d + 0] = _mm256_unpackhi_epi64(v_tmp32_lo[s + 0], v_tmp32_lo[s + 2]);
    v_tmp64_hi[d + 1] = _mm256_unpackhi_epi64(v_tmp32_lo[s + 1], v_tmp32_lo[s + 3]);
    v_tmp64_hi[d + 4] = _mm256_unpackhi_epi64(v_tmp32_hi[s + 0], v_tmp32_hi[s + 2]);
    v_tmp64_hi[d + 5] = _mm256_unpackhi_epi64(v_tmp32_hi[s + 1], v_tmp32_hi[s + 3]);
  }

  dst[0] = _mm256_permute2x128_si256(v_tmp64_lo[0], v_tmp64_hi[0], 0x20);
  dst[1] = _mm256_permute2x128_si256(v_tmp64_lo[4], v_tmp64_hi[4], 0x20);
  dst[2] = _mm256_permute2x128_si256(v_tmp64_lo[2], v_tmp64_hi[2], 0x20);
  dst[3] = _mm256_permute2x128_si256(v_tmp64_lo[6], v_tmp64_hi[6], 0x20);

  dst[4] = _mm256_permute2x128_si256(v_tmp64_lo[0], v_tmp64_hi[0], 0x31);
  dst[5] = _mm256_permute2x128_si256(v_tmp64_lo[4], v_tmp64_hi[4], 0x31);
  dst[6] = _mm256_permute2x128_si256(v_tmp64_lo[2], v_tmp64_hi[2], 0x31);
  dst[7] = _mm256_permute2x128_si256(v_tmp64_lo[6], v_tmp64_hi[6], 0x31);

  dst[8]  = _mm256_permute2x128_si256(v_tmp64_lo[1], v_tmp64_hi[1], 0x20);
  dst[9]  = _mm256_permute2x128_si256(v_tmp64_lo[5], v_tmp64_hi[5], 0x20);
  dst[10] = _mm256_permute2x128_si256(v_tmp64_lo[3], v_tmp64_hi[3], 0x20);
  dst[11] = _mm256_permute2x128_si256(v_tmp64_lo[7], v_tmp64_hi[7], 0x20);

  dst[12] = _mm256_permute2x128_si256(v_tmp64_lo[1], v_tmp64_hi[1], 0x31);
  dst[13] = _mm256_permute2x128_si256(v_tmp64_lo[5], v_tmp64_hi[5], 0x31);
  dst[14] = _mm256_permute2x128_si256(v_tmp64_lo[3], v_tmp64_hi[3], 0x31);
  dst[15] = _mm256_permute2x128_si256(v_tmp64_lo[7], v_tmp64_hi[7], 0x31);
}
static void transpose_32x16_avx2(const __m256i* src, __m256i* dst) {
  transpose_16x16_avx2_stride((int16_t const *)src,                        (int16_t *)dst, 32, 16);
  transpose_16x16_avx2_stride((int16_t const *)src + 16,           (int16_t *)dst + 16 * 16, 32, 16);
}
static void transpose_32x32_avx2(const __m256i* src, __m256i* dst) {
  transpose_16x16_avx2_stride((int16_t const *)src,                        (int16_t *)dst, 32, 32);
  transpose_16x16_avx2_stride((int16_t const *)src + 16,           (int16_t *)dst + 16 * 32, 32, 32);
  transpose_16x16_avx2_stride((int16_t const *)src + 16 * 32,      (int16_t *)dst + 16, 32, 32);
  transpose_16x16_avx2_stride((int16_t const *)src + 16 * 32 + 16, (int16_t *)dst + 16 * 32 + 16, 32, 32);
}
static void transpose_32x64_avx2(const __m256i* src, __m256i* dst){}
static void transpose_64x2_avx2(const __m256i* src, __m256i* dst){}
static void transpose_64x4_avx2(const __m256i* src, __m256i* dst){}
static void transpose_64x8_avx2(const __m256i* src, __m256i* dst){}
static void transpose_64x16_avx2(const __m256i* src, __m256i* dst){}
static void transpose_64x32_avx2(const __m256i* src, __m256i* dst){}
static void transpose_64x64_avx2(const __m256i* src, __m256i* dst){}



static transpose_func* transpose_func_table[6][6] = {
  { transpose_2x2_avx2,  transpose_4x2_avx2,  transpose_8x2_avx2,  transpose_16x2_avx2,  transpose_32x2_avx2,  transpose_64x2_avx2},
  { transpose_2x4_avx2,  transpose_4x4_avx2,  transpose_8x4_avx2,  transpose_16x4_avx2,  transpose_32x4_avx2,  transpose_64x4_avx2},
  { transpose_2x8_avx2,  transpose_4x8_avx2,  transpose_8x8_avx2,  transpose_16x8_avx2,  transpose_32x8_avx2,  transpose_64x8_avx2},
  {transpose_2x16_avx2, transpose_4x16_avx2, transpose_8x16_avx2, transpose_16x16_avx2, transpose_32x16_avx2, transpose_64x16_avx2},
  {transpose_2x32_avx2, transpose_4x32_avx2, transpose_8x32_avx2, transpose_16x32_avx2, transpose_32x32_avx2, transpose_64x32_avx2},
  {transpose_2x64_avx2, transpose_4x64_avx2, transpose_8x64_avx2, transpose_16x64_avx2, transpose_32x64_avx2, transpose_64x64_avx2},
};


// Dispatcher function for avx2 transposes. This calls the proper subfunction
static void transpose_avx2(const __m256i* src, __m256i* dst, const int width, const int height)
{
  // No need to transpose something of width or height 1
  const int w_log2_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int h_log2_minus1 = uvg_g_convert_to_log2[height] - 1;
  
  transpose_func* func = transpose_func_table[h_log2_minus1][w_log2_minus1];
  func(src, dst);
}


// 4x4 matrix multiplication with value clipping.
// Parameters: Two 4x4 matrices containing 16-bit values in consecutive addresses,
//             destination for the result and the shift value for clipping.
static __m256i mul_clip_matrix_4x4_avx2(const __m256i left, const __m256i right, int shift)
{
  const int32_t add    = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  __m256i right_los = _mm256_permute4x64_epi64(right, _MM_SHUFFLE(2, 0, 2, 0));
  __m256i right_his = _mm256_permute4x64_epi64(right, _MM_SHUFFLE(3, 1, 3, 1));

  __m256i right_cols_up = _mm256_unpacklo_epi16(right_los, right_his);
  __m256i right_cols_dn = _mm256_unpackhi_epi16(right_los, right_his);

  __m256i left_slice1 = _mm256_shuffle_epi32(left, _MM_SHUFFLE(0, 0, 0, 0));
  __m256i left_slice2 = _mm256_shuffle_epi32(left, _MM_SHUFFLE(1, 1, 1, 1));
  __m256i left_slice3 = _mm256_shuffle_epi32(left, _MM_SHUFFLE(2, 2, 2, 2));
  __m256i left_slice4 = _mm256_shuffle_epi32(left, _MM_SHUFFLE(3, 3, 3, 3));

  __m256i prod1 = _mm256_madd_epi16(left_slice1, right_cols_up);
  __m256i prod2 = _mm256_madd_epi16(left_slice2, right_cols_dn);
  __m256i prod3 = _mm256_madd_epi16(left_slice3, right_cols_up);
  __m256i prod4 = _mm256_madd_epi16(left_slice4, right_cols_dn);

  __m256i rows_up = _mm256_add_epi32(prod1, prod2);
  __m256i rows_dn = _mm256_add_epi32(prod3, prod4);

  __m256i rows_up_tr = truncate_avx2(rows_up, debias, shift);
  __m256i rows_dn_tr = truncate_avx2(rows_dn, debias, shift);

  __m256i result = _mm256_packs_epi32(rows_up_tr, rows_dn_tr);
  return result;
}

static void matrix_dst_4x4_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = uvg_g_convert_to_bit[4] + 1 + (bitdepth - 8);
  int32_t shift_2nd = uvg_g_convert_to_bit[4] + 8;
  const int16_t *tdst = &uvg_g_dst_4_t[0][0];
  const int16_t *dst  = &uvg_g_dst_4  [0][0];

  __m256i tdst_v = _mm256_load_si256((const __m256i *) tdst);
  __m256i  dst_v = _mm256_load_si256((const __m256i *)  dst);
  __m256i   in_v = _mm256_load_si256((const __m256i *)input);

  __m256i tmp    = mul_clip_matrix_4x4_avx2(in_v,  tdst_v, shift_1st);
  __m256i result = mul_clip_matrix_4x4_avx2(dst_v, tmp,    shift_2nd);

  _mm256_store_si256((__m256i *)output, result);
}

static void matrix_idst_4x4_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);

  const int16_t *tdst = &uvg_g_dst_4_t[0][0];
  const int16_t *dst  = &uvg_g_dst_4  [0][0];

  __m256i tdst_v = _mm256_load_si256((const __m256i *)tdst);
  __m256i  dst_v = _mm256_load_si256((const __m256i *) dst);
  __m256i   in_v = _mm256_load_si256((const __m256i *)input);

  __m256i tmp    = mul_clip_matrix_4x4_avx2(tdst_v, in_v,  shift_1st);
  __m256i result = mul_clip_matrix_4x4_avx2(tmp,    dst_v, shift_2nd);

  _mm256_store_si256((__m256i *)output, result);
}

static void matrix_dct_4x4_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = uvg_g_convert_to_bit[4] + 1 + (bitdepth - 8);
  int32_t shift_2nd = uvg_g_convert_to_bit[4] + 8;
  const int16_t *tdct = &uvg_g_dct_4_t[0][0];
  const int16_t *dct  = &uvg_g_dct_4  [0][0];

  __m256i tdct_v = _mm256_load_si256((const __m256i *) tdct);
  __m256i  dct_v = _mm256_load_si256((const __m256i *)  dct);
  __m256i   in_v = _mm256_load_si256((const __m256i *)input);

  __m256i tmp    = mul_clip_matrix_4x4_avx2(in_v,  tdct_v, shift_1st);
  __m256i result = mul_clip_matrix_4x4_avx2(dct_v, tmp,    shift_2nd);

  _mm256_store_si256((__m256i *)output, result);
}

static void matrix_idct_4x4_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);

  const int16_t *tdct = &uvg_g_dct_4_t[0][0];
  const int16_t *dct  = &uvg_g_dct_4  [0][0];

  __m256i tdct_v = _mm256_load_si256((const __m256i *)tdct);
  __m256i  dct_v = _mm256_load_si256((const __m256i *) dct);
  __m256i   in_v = _mm256_load_si256((const __m256i *)input);

  __m256i tmp    = mul_clip_matrix_4x4_avx2(tdct_v, in_v,  shift_1st);
  __m256i result = mul_clip_matrix_4x4_avx2(tmp,    dct_v, shift_2nd);

  _mm256_store_si256((__m256i *)output, result);
}

static void mul_clip_matrix_8x8_avx2(const int16_t *left, const int16_t *right, int16_t *dst, const int32_t shift)
{
  const __m256i transp_mask = _mm256_broadcastsi128_si256(_mm_setr_epi8(0, 1, 8, 9, 2, 3, 10, 11, 4, 5, 12, 13, 6, 7, 14, 15));

  const int32_t add    = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  __m256i left_dr[4] = {
    _mm256_load_si256((const __m256i *)left + 0),
    _mm256_load_si256((const __m256i *)left + 1),
    _mm256_load_si256((const __m256i *)left + 2),
    _mm256_load_si256((const __m256i *)left + 3),
  };
  __m256i right_dr[4] = {
    _mm256_load_si256((const __m256i *)right + 0),
    _mm256_load_si256((const __m256i *)right + 1),
    _mm256_load_si256((const __m256i *)right + 2),
    _mm256_load_si256((const __m256i *)right + 3),
  };

  __m256i rdrs_rearr[8];

  // Rearrange right matrix
  for (int32_t dry = 0; dry < 4; dry++) {
    __m256i rdr = right_dr[dry];
    __m256i rdr_los = _mm256_permute4x64_epi64(rdr, _MM_SHUFFLE(2, 0, 2, 0));
    __m256i rdr_his = _mm256_permute4x64_epi64(rdr, _MM_SHUFFLE(3, 1, 3, 1));

    __m256i rdr_lo_rearr = _mm256_shuffle_epi8(rdr_los, transp_mask);
    __m256i rdr_hi_rearr = _mm256_shuffle_epi8(rdr_his, transp_mask);

    rdrs_rearr[dry * 2 + 0] = rdr_lo_rearr;
    rdrs_rearr[dry * 2 + 1] = rdr_hi_rearr;
  }

  // Double-Row Y for destination matrix
  for (int32_t dry = 0; dry < 4; dry++) {
    __m256i ldr = left_dr[dry];

    __m256i ldr_slice12 = _mm256_shuffle_epi32(ldr, _MM_SHUFFLE(0, 0, 0, 0));
    __m256i ldr_slice34 = _mm256_shuffle_epi32(ldr, _MM_SHUFFLE(1, 1, 1, 1));
    __m256i ldr_slice56 = _mm256_shuffle_epi32(ldr, _MM_SHUFFLE(2, 2, 2, 2));
    __m256i ldr_slice78 = _mm256_shuffle_epi32(ldr, _MM_SHUFFLE(3, 3, 3, 3));

    __m256i prod1 = _mm256_madd_epi16(ldr_slice12, rdrs_rearr[0]);
    __m256i prod2 = _mm256_madd_epi16(ldr_slice12, rdrs_rearr[1]);
    __m256i prod3 = _mm256_madd_epi16(ldr_slice34, rdrs_rearr[2]);
    __m256i prod4 = _mm256_madd_epi16(ldr_slice34, rdrs_rearr[3]);
    __m256i prod5 = _mm256_madd_epi16(ldr_slice56, rdrs_rearr[4]);
    __m256i prod6 = _mm256_madd_epi16(ldr_slice56, rdrs_rearr[5]);
    __m256i prod7 = _mm256_madd_epi16(ldr_slice78, rdrs_rearr[6]);
    __m256i prod8 = _mm256_madd_epi16(ldr_slice78, rdrs_rearr[7]);

    __m256i lo_1 = _mm256_add_epi32(prod1, prod3);
    __m256i hi_1 = _mm256_add_epi32(prod2, prod4);
    __m256i lo_2 = _mm256_add_epi32(prod5, prod7);
    __m256i hi_2 = _mm256_add_epi32(prod6, prod8);

    __m256i lo   = _mm256_add_epi32(lo_1,  lo_2);
    __m256i hi   = _mm256_add_epi32(hi_1,  hi_2);

    __m256i lo_tr = truncate_avx2(lo, debias, shift);
    __m256i hi_tr = truncate_avx2(hi, debias, shift);

    __m256i final_dr = _mm256_packs_epi32(lo_tr, hi_tr);

    _mm256_store_si256((__m256i *)dst + dry, final_dr);
  }
}

// Multiplies A by B_T's transpose and stores result's transpose in output,
// which should be an array of 4 __m256i's
static void matmul_8x8_a_bt_t(const int16_t *a, const int16_t *b_t,
    __m256i *output, const int8_t shift)
{
  const int32_t add    = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  // Keep upper row intact and swap neighboring 16-bit words in lower row
  const __m256i shuf_lorow_mask =
      _mm256_setr_epi8(0,  1,  2,  3,  4,  5,  6,  7,
                       8,  9,  10, 11, 12, 13, 14, 15,
                       18, 19, 16, 17, 22, 23, 20, 21,
                       26, 27, 24, 25, 30, 31, 28, 29);

  const __m256i *b_t_256 = (const __m256i *)b_t;

  // Dual Rows, because two 8x16b words fit in one YMM
  __m256i a_dr_0      = _mm256_load_si256((__m256i *)a + 0);
  __m256i a_dr_1      = _mm256_load_si256((__m256i *)a + 1);
  __m256i a_dr_2      = _mm256_load_si256((__m256i *)a + 2);
  __m256i a_dr_3      = _mm256_load_si256((__m256i *)a + 3);

  __m256i a_dr_0_swp  = swap_lanes(a_dr_0);
  __m256i a_dr_1_swp  = swap_lanes(a_dr_1);
  __m256i a_dr_2_swp  = swap_lanes(a_dr_2);
  __m256i a_dr_3_swp  = swap_lanes(a_dr_3);

  for (int dry = 0; dry < 4; dry++) {

    // Read dual columns of B matrix by reading rows of its transpose
    __m256i b_dc        = _mm256_load_si256(b_t_256 + dry);

    __m256i prod0       = _mm256_madd_epi16(b_dc,     a_dr_0);
    __m256i prod0_swp   = _mm256_madd_epi16(b_dc,     a_dr_0_swp);
    __m256i prod1       = _mm256_madd_epi16(b_dc,     a_dr_1);
    __m256i prod1_swp   = _mm256_madd_epi16(b_dc,     a_dr_1_swp);
    __m256i prod2       = _mm256_madd_epi16(b_dc,     a_dr_2);
    __m256i prod2_swp   = _mm256_madd_epi16(b_dc,     a_dr_2_swp);
    __m256i prod3       = _mm256_madd_epi16(b_dc,     a_dr_3);
    __m256i prod3_swp   = _mm256_madd_epi16(b_dc,     a_dr_3_swp);

    __m256i hsum0       = _mm256_hadd_epi32(prod0,    prod0_swp);
    __m256i hsum1       = _mm256_hadd_epi32(prod1,    prod1_swp);
    __m256i hsum2       = _mm256_hadd_epi32(prod2,    prod2_swp);
    __m256i hsum3       = _mm256_hadd_epi32(prod3,    prod3_swp);

    __m256i hsum2c_0    = _mm256_hadd_epi32(hsum0,    hsum1);
    __m256i hsum2c_1    = _mm256_hadd_epi32(hsum2,    hsum3);

    __m256i hsum2c_0_tr = truncate_avx2(hsum2c_0, debias, shift);
    __m256i hsum2c_1_tr = truncate_avx2(hsum2c_1, debias, shift);

    __m256i tmp_dc      = _mm256_packs_epi32(hsum2c_0_tr, hsum2c_1_tr);

    output[dry]         = _mm256_shuffle_epi8(tmp_dc, shuf_lorow_mask);
  }
}

// Multiplies A by B_T's transpose and stores result in output
// which should be an array of 4 __m256i's
static void matmul_8x8_a_bt(const int16_t *a, const __m256i *b_t,
    int16_t *output, const int8_t shift)
{
  const int32_t add    = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i shuf_lorow_mask =
      _mm256_setr_epi8(0,  1,  2,  3,  4,  5,  6,  7,
                       8,  9,  10, 11, 12, 13, 14, 15,
                       18, 19, 16, 17, 22, 23, 20, 21,
                       26, 27, 24, 25, 30, 31, 28, 29);

  const __m256i *a_256 = (const __m256i *)a;

  __m256i b_dc_0      = b_t[0];
  __m256i b_dc_1      = b_t[1];
  __m256i b_dc_2      = b_t[2];
  __m256i b_dc_3      = b_t[3];

  __m256i b_dc_0_swp  = swap_lanes(b_dc_0);
  __m256i b_dc_1_swp  = swap_lanes(b_dc_1);
  __m256i b_dc_2_swp  = swap_lanes(b_dc_2);
  __m256i b_dc_3_swp  = swap_lanes(b_dc_3);

  for (int dry = 0; dry < 4; dry++) {
    __m256i a_dr        = _mm256_load_si256(a_256 + dry);

    __m256i prod0       = _mm256_madd_epi16(a_dr,     b_dc_0);
    __m256i prod0_swp   = _mm256_madd_epi16(a_dr,     b_dc_0_swp);
    __m256i prod1       = _mm256_madd_epi16(a_dr,     b_dc_1);
    __m256i prod1_swp   = _mm256_madd_epi16(a_dr,     b_dc_1_swp);
    __m256i prod2       = _mm256_madd_epi16(a_dr,     b_dc_2);
    __m256i prod2_swp   = _mm256_madd_epi16(a_dr,     b_dc_2_swp);
    __m256i prod3       = _mm256_madd_epi16(a_dr,     b_dc_3);
    __m256i prod3_swp   = _mm256_madd_epi16(a_dr,     b_dc_3_swp);

    __m256i hsum0       = _mm256_hadd_epi32(prod0,    prod0_swp);
    __m256i hsum1       = _mm256_hadd_epi32(prod1,    prod1_swp);
    __m256i hsum2       = _mm256_hadd_epi32(prod2,    prod2_swp);
    __m256i hsum3       = _mm256_hadd_epi32(prod3,    prod3_swp);

    __m256i hsum2c_0    = _mm256_hadd_epi32(hsum0,    hsum1);
    __m256i hsum2c_1    = _mm256_hadd_epi32(hsum2,    hsum3);

    __m256i hsum2c_0_tr = truncate_avx2(hsum2c_0, debias, shift);
    __m256i hsum2c_1_tr = truncate_avx2(hsum2c_1, debias, shift);

    __m256i tmp_dr      = _mm256_packs_epi32(hsum2c_0_tr, hsum2c_1_tr);

    __m256i final_dr    = _mm256_shuffle_epi8(tmp_dr, shuf_lorow_mask);

    _mm256_store_si256((__m256i *)output + dry, final_dr);
  }
}

static void matrix_dct_8x8_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = uvg_g_convert_to_bit[8] + 1 + (bitdepth - 8);
  int32_t shift_2nd = uvg_g_convert_to_bit[8] + 8;

  const int16_t *dct  = &uvg_g_dct_8[0][0];

  /*
   * Multiply input by the tranpose of DCT matrix into tmpres, and DCT matrix
   * by tmpres - this is then our output matrix
   *
   * It's easier to implement an AVX2 matrix multiplication if you can multiply
   * the left term with the transpose of the right term. Here things are stored
   * row-wise, not column-wise, so we can effectively read DCT_T column-wise
   * into YMM registers by reading DCT row-wise. Also because of this, the
   * first multiplication is hacked to produce the transpose of the result
   * instead, since it will be used in similar fashion as the right operand
   * in the second multiplication.
   */

  __m256i tmpres[4];

  matmul_8x8_a_bt_t(input,  dct, tmpres, shift_1st);
  matmul_8x8_a_bt  (dct, tmpres, output, shift_2nd);
}

static void matrix_idct_8x8_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);
  ALIGNED(64) int16_t tmp[8 * 8];

  const int16_t *tdct = &uvg_g_dct_8_t[0][0];
  const int16_t *dct  = &uvg_g_dct_8  [0][0];

  mul_clip_matrix_8x8_avx2(tdct, input, tmp,    shift_1st);
  mul_clip_matrix_8x8_avx2(tmp,  dct,   output, shift_2nd);

  /*
   * Because:
   * out = tdct * input * dct = tdct * (input * dct) = tdct * (input * transpose(tdct))
   * This could almost be done this way:
   *
   * matmul_8x8_a_bt_t(input, tdct, debias1, shift_1st, tmp);
   * matmul_8x8_a_bt  (tdct,  tmp,  debias2, shift_2nd, output);
   *
   * But not really, since it will fall victim to some very occasional
   * rounding errors. Sadly.
   */
}

static void matmul_16x16_a_bt(const __m256i *a,
                              const __m256i *b_t,
                                    __m256i *output,
                              const int32_t  shift)
{
  const int32_t add    = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  for (int32_t y = 0; y < 16; y++) {
    __m256i a_r = a[y];
    __m256i results_32[2];

    for (int32_t fco = 0; fco < 2; fco++) {
      // Read first cols 0, 1, 2, 3, 8, 9, 10, 11, and then next 4
      __m256i bt_c0  = b_t[fco * 4 + 0];
      __m256i bt_c1  = b_t[fco * 4 + 1];
      __m256i bt_c2  = b_t[fco * 4 + 2];
      __m256i bt_c3  = b_t[fco * 4 + 3];
      __m256i bt_c8  = b_t[fco * 4 + 8];
      __m256i bt_c9  = b_t[fco * 4 + 9];
      __m256i bt_c10 = b_t[fco * 4 + 10];
      __m256i bt_c11 = b_t[fco * 4 + 11];

      __m256i p0  = _mm256_madd_epi16(a_r, bt_c0);
      __m256i p1  = _mm256_madd_epi16(a_r, bt_c1);
      __m256i p2  = _mm256_madd_epi16(a_r, bt_c2);
      __m256i p3  = _mm256_madd_epi16(a_r, bt_c3);
      __m256i p8  = _mm256_madd_epi16(a_r, bt_c8);
      __m256i p9  = _mm256_madd_epi16(a_r, bt_c9);
      __m256i p10 = _mm256_madd_epi16(a_r, bt_c10);
      __m256i p11 = _mm256_madd_epi16(a_r, bt_c11);

      // Combine low lanes from P0 and P8, high lanes from them, and the same
      // with P1:P9 and so on
      __m256i p0l = _mm256_permute2x128_si256(p0, p8,  0x20);
      __m256i p0h = _mm256_permute2x128_si256(p0, p8,  0x31);
      __m256i p1l = _mm256_permute2x128_si256(p1, p9,  0x20);
      __m256i p1h = _mm256_permute2x128_si256(p1, p9,  0x31);
      __m256i p2l = _mm256_permute2x128_si256(p2, p10, 0x20);
      __m256i p2h = _mm256_permute2x128_si256(p2, p10, 0x31);
      __m256i p3l = _mm256_permute2x128_si256(p3, p11, 0x20);
      __m256i p3h = _mm256_permute2x128_si256(p3, p11, 0x31);

      __m256i s0  = _mm256_add_epi32(p0l, p0h);
      __m256i s1  = _mm256_add_epi32(p1l, p1h);
      __m256i s2  = _mm256_add_epi32(p2l, p2h);
      __m256i s3  = _mm256_add_epi32(p3l, p3h);

      __m256i s4  = _mm256_unpacklo_epi64(s0, s1);
      __m256i s5  = _mm256_unpackhi_epi64(s0, s1);
      __m256i s6  = _mm256_unpacklo_epi64(s2, s3);
      __m256i s7  = _mm256_unpackhi_epi64(s2, s3);

      __m256i s8  = _mm256_add_epi32(s4, s5);
      __m256i s9  = _mm256_add_epi32(s6, s7);

      __m256i res = _mm256_hadd_epi32(s8, s9);
      results_32[fco] = truncate_avx2(res, debias, shift);
    }
    output[y] = _mm256_packs_epi32(results_32[0], results_32[1]);
  }
}

// NOTE: The strides measured by s_stride_log2 and d_stride_log2 are in units
// of 16 coeffs, not 1!
static void transpose_16x16_stride(const int16_t *src,
                                         int16_t *dst,
                                         uint8_t  s_stride_log2,
                                         uint8_t  d_stride_log2)
{
  __m256i tmp_128[16];
  for (uint32_t i = 0; i < 16; i += 8) {

    // After every n-bit unpack, 2n-bit units in the vectors will be in
    // correct order. Pair words first, then dwords, then qwords. After that,
    // whole lanes will be correct.
    __m256i tmp_32[8];
    __m256i tmp_64[8];

    __m256i m[8] = {
      _mm256_load_si256((const __m256i *)src + ((i + 0) << s_stride_log2)),
      _mm256_load_si256((const __m256i *)src + ((i + 1) << s_stride_log2)),
      _mm256_load_si256((const __m256i *)src + ((i + 2) << s_stride_log2)),
      _mm256_load_si256((const __m256i *)src + ((i + 3) << s_stride_log2)),
      _mm256_load_si256((const __m256i *)src + ((i + 4) << s_stride_log2)),
      _mm256_load_si256((const __m256i *)src + ((i + 5) << s_stride_log2)),
      _mm256_load_si256((const __m256i *)src + ((i + 6) << s_stride_log2)),
      _mm256_load_si256((const __m256i *)src + ((i + 7) << s_stride_log2)),
    };

    tmp_32[0]      = _mm256_unpacklo_epi16(     m[0],      m[1]);
    tmp_32[1]      = _mm256_unpacklo_epi16(     m[2],      m[3]);
    tmp_32[2]      = _mm256_unpackhi_epi16(     m[0],      m[1]);
    tmp_32[3]      = _mm256_unpackhi_epi16(     m[2],      m[3]);

    tmp_32[4]      = _mm256_unpacklo_epi16(     m[4],      m[5]);
    tmp_32[5]      = _mm256_unpacklo_epi16(     m[6],      m[7]);
    tmp_32[6]      = _mm256_unpackhi_epi16(     m[4],      m[5]);
    tmp_32[7]      = _mm256_unpackhi_epi16(     m[6],      m[7]);


    tmp_64[0]      = _mm256_unpacklo_epi32(tmp_32[0], tmp_32[1]);
    tmp_64[1]      = _mm256_unpacklo_epi32(tmp_32[2], tmp_32[3]);
    tmp_64[2]      = _mm256_unpackhi_epi32(tmp_32[0], tmp_32[1]);
    tmp_64[3]      = _mm256_unpackhi_epi32(tmp_32[2], tmp_32[3]);

    tmp_64[4]      = _mm256_unpacklo_epi32(tmp_32[4], tmp_32[5]);
    tmp_64[5]      = _mm256_unpacklo_epi32(tmp_32[6], tmp_32[7]);
    tmp_64[6]      = _mm256_unpackhi_epi32(tmp_32[4], tmp_32[5]);
    tmp_64[7]      = _mm256_unpackhi_epi32(tmp_32[6], tmp_32[7]);


    tmp_128[i + 0] = _mm256_unpacklo_epi64(tmp_64[0], tmp_64[4]);
    tmp_128[i + 1] = _mm256_unpackhi_epi64(tmp_64[0], tmp_64[4]);
    tmp_128[i + 2] = _mm256_unpacklo_epi64(tmp_64[2], tmp_64[6]);
    tmp_128[i + 3] = _mm256_unpackhi_epi64(tmp_64[2], tmp_64[6]);

    tmp_128[i + 4] = _mm256_unpacklo_epi64(tmp_64[1], tmp_64[5]);
    tmp_128[i + 5] = _mm256_unpackhi_epi64(tmp_64[1], tmp_64[5]);
    tmp_128[i + 6] = _mm256_unpacklo_epi64(tmp_64[3], tmp_64[7]);
    tmp_128[i + 7] = _mm256_unpackhi_epi64(tmp_64[3], tmp_64[7]);
  }

  for (uint32_t i = 0; i < 8; i++) {
    uint32_t loid     = i + 0;
    uint32_t hiid     = i + 8;

    uint32_t dst_loid = loid << d_stride_log2;
    uint32_t dst_hiid = hiid << d_stride_log2;

    __m256i lo       = tmp_128[loid];
    __m256i hi       = tmp_128[hiid];
    __m256i final_lo = _mm256_permute2x128_si256(lo, hi, 0x20);
    __m256i final_hi = _mm256_permute2x128_si256(lo, hi, 0x31);

    _mm256_store_si256((__m256i *)dst + dst_loid, final_lo);
    _mm256_store_si256((__m256i *)dst + dst_hiid, final_hi);
  }
}

static void transpose_16x16(const int16_t *src, int16_t *dst)
{
  transpose_16x16_stride(src, dst, 0, 0);
}

static __m256i truncate_inv(__m256i v, int32_t shift)
{
  int32_t add = 1 << (shift - 1);

  __m256i debias  = _mm256_set1_epi32(add);
  __m256i v2      = _mm256_add_epi32 (v,  debias);
  __m256i trunced = _mm256_srai_epi32(v2, shift);
  return  trunced;
}

static __m256i extract_odds(__m256i v)
{
  // 0 1 2 3 4 5 6 7 | 8 9 a b c d e f => 1 3 5 7 1 3 5 7 | 9 b d f 9 b d f
  const __m256i oddmask = _mm256_setr_epi8( 2,  3,  6,  7, 10, 11, 14, 15,
                                            2,  3,  6,  7, 10, 11, 14, 15,
                                            2,  3,  6,  7, 10, 11, 14, 15,
                                            2,  3,  6,  7, 10, 11, 14, 15);

  __m256i tmp = _mm256_shuffle_epi8 (v,   oddmask);
  return _mm256_permute4x64_epi64   (tmp, _MM_SHUFFLE(3, 1, 2, 0));
}

static __m256i extract_combine_odds(__m256i v0, __m256i v1)
{
  // 0 1 2 3 4 5 6 7 | 8 9 a b c d e f => 1 3 5 7 1 3 5 7 | 9 b d f 9 b d f
  const __m256i oddmask = _mm256_setr_epi8( 2,  3,  6,  7, 10, 11, 14, 15,
                                            2,  3,  6,  7, 10, 11, 14, 15,
                                            2,  3,  6,  7, 10, 11, 14, 15,
                                            2,  3,  6,  7, 10, 11, 14, 15);

  __m256i tmp0 = _mm256_shuffle_epi8(v0,   oddmask);
  __m256i tmp1 = _mm256_shuffle_epi8(v1,   oddmask);

  __m256i tmp2 = _mm256_blend_epi32 (tmp0, tmp1, 0xcc); // 1100 1100

  return _mm256_permute4x64_epi64   (tmp2, _MM_SHUFFLE(3, 1, 2, 0));
}

// Extract items 2, 6, A and E from first four columns of DCT, order them as
// follows:
// D0,2 D0,6 D1,2 D1,6 D1,a D1,e D0,a D0,e | D2,2 D2,6 D3,2 D3,6 D3,a D3,e D2,a D2,e
static __m256i extract_26ae(const __m256i *tdct)
{
  // 02 03 22 23 06 07 26 27 | 0a 0b 2a 2b 02 0f 2e 2f
  // =>
  // 02 06 22 26 02 06 22 26 | 2a 2e 0a 0e 2a 2e 0a 0e
  const __m256i evens_mask = _mm256_setr_epi8( 0,  1,  8,  9,  4,  5, 12, 13,
                                               0,  1,  8,  9,  4,  5, 12, 13,
                                               4,  5, 12, 13,  0,  1,  8,  9,
                                               4,  5, 12, 13,  0,  1,  8,  9);

  __m256i shufd_0 = _mm256_shuffle_epi32(tdct[0], _MM_SHUFFLE(2, 3, 0, 1));
  __m256i shufd_2 = _mm256_shuffle_epi32(tdct[2], _MM_SHUFFLE(2, 3, 0, 1));

  __m256i cmbd_01 = _mm256_blend_epi32(shufd_0, tdct[1], 0xaa); // 1010 1010
  __m256i cmbd_23 = _mm256_blend_epi32(shufd_2, tdct[3], 0xaa); // 1010 1010

  __m256i evens_01 = _mm256_shuffle_epi8(cmbd_01, evens_mask);
  __m256i evens_23 = _mm256_shuffle_epi8(cmbd_23, evens_mask);

  __m256i evens_0123 = _mm256_unpacklo_epi64(evens_01, evens_23);

  return _mm256_permute4x64_epi64(evens_0123, _MM_SHUFFLE(3, 1, 2, 0));
}

// 2 6 2 6 a e a e | 2 6 2 6 a e a e
static __m256i extract_26ae_vec(__m256i col)
{
  const __m256i mask_26ae = _mm256_set1_epi32(0x0d0c0504);

  // 2 6 2 6 2 6 2 6 | a e a e a e a e
  __m256i reord = _mm256_shuffle_epi8     (col,   mask_26ae);
  __m256i final = _mm256_permute4x64_epi64(reord, _MM_SHUFFLE(3, 1, 2, 0));
  return  final;
}

// D00 D80 D01 D81 D41 Dc1 D40 Dc0 | D40 Dc0 D41 Dc1 D01 D81 D00 D80
static __m256i extract_d048c(const __m256i *tdct)
{
  const __m256i final_shuf = _mm256_setr_epi8( 0,  1,  8,  9,  2,  3, 10, 11,
                                               6,  7, 14, 15,  4,  5, 12, 13,
                                               4,  5, 12, 13,  6,  7, 14, 15,
                                               2,  3, 10, 11,  0,  1,  8,  9);
  __m256i c0 = tdct[0];
  __m256i c1 = tdct[1];

  __m256i c1_2  = _mm256_slli_epi32       (c1,    16);
  __m256i cmbd  = _mm256_blend_epi16      (c0,    c1_2, 0x22); // 0010 0010
  __m256i cmbd2 = _mm256_shuffle_epi32    (cmbd,  _MM_SHUFFLE(2, 0, 2, 0));
  __m256i cmbd3 = _mm256_permute4x64_epi64(cmbd2, _MM_SHUFFLE(3, 1, 2, 0));
  __m256i final = _mm256_shuffle_epi8     (cmbd3, final_shuf);

  return final;
}

// 0 8 0 8 4 c 4 c | 4 c 4 c 0 8 0 8
static __m256i extract_d048c_vec(__m256i col)
{
  const __m256i shufmask = _mm256_setr_epi8( 0,  1,  0,  1,  8,  9,  8,  9,
                                             8,  9,  8,  9,  0,  1,  0,  1,
                                             0,  1,  0,  1,  8,  9,  8,  9,
                                             8,  9,  8,  9,  0,  1,  0,  1);

  __m256i col_db4s = _mm256_shuffle_epi8     (col, shufmask);
  __m256i col_los  = _mm256_permute4x64_epi64(col_db4s, _MM_SHUFFLE(1, 1, 0, 0));
  __m256i col_his  = _mm256_permute4x64_epi64(col_db4s, _MM_SHUFFLE(3, 3, 2, 2));

  __m256i final    = _mm256_unpacklo_epi16   (col_los,  col_his);
  return final;
}

static void partial_butterfly_inverse_16_avx2(const int16_t *src, int16_t *dst, int32_t shift)
{
  __m256i tsrc[16];

  const uint32_t width = 16;

  const int16_t *tdct = &uvg_g_dct_16_t[0][0];

  const __m256i  eo_signmask = _mm256_setr_epi32( 1,  1,  1,  1, -1, -1, -1, -1);
  const __m256i eeo_signmask = _mm256_setr_epi32( 1,  1, -1, -1, -1, -1,  1,  1);
  const __m256i   o_signmask = _mm256_set1_epi32(-1);

  const __m256i final_shufmask = _mm256_setr_epi8( 0,  1,  2,  3,  4,  5,  6,  7,
                                                   8,  9, 10, 11, 12, 13, 14, 15,
                                                   6,  7,  4,  5,  2,  3,  0,  1,
                                                  14, 15, 12, 13, 10, 11,  8,  9);
  transpose_16x16(src, (int16_t *)tsrc);

  const __m256i dct_cols[8] = {
    _mm256_load_si256((const __m256i *)tdct + 0),
    _mm256_load_si256((const __m256i *)tdct + 1),
    _mm256_load_si256((const __m256i *)tdct + 2),
    _mm256_load_si256((const __m256i *)tdct + 3),
    _mm256_load_si256((const __m256i *)tdct + 4),
    _mm256_load_si256((const __m256i *)tdct + 5),
    _mm256_load_si256((const __m256i *)tdct + 6),
    _mm256_load_si256((const __m256i *)tdct + 7),
  };

  // These contain: D1,0 D3,0 D5,0 D7,0 D9,0 Db,0 Dd,0 Df,0 | D1,4 D3,4 D5,4 D7,4 D9,4 Db,4 Dd,4 Df,4
  //                D1,1 D3,1 D5,1 D7,1 D9,1 Db,1 Dd,1 Df,1 | D1,5 D3,5 D5,5 D7,5 D9,5 Db,5 Dd,5 Df,5
  //                D1,2 D3,2 D5,2 D7,2 D9,2 Db,2 Dd,2 Df,2 | D1,6 D3,6 D5,6 D7,6 D9,6 Db,6 Dd,6 Df,6
  //                D1,3 D3,3 D5,3 D7,3 D9,3 Db,3 Dd,3 Df,3 | D1,7 D3,7 D5,7 D7,7 D9,7 Db,7 Dd,7 Df,7
  __m256i dct_col_odds[4];
  for (uint32_t j = 0; j < 4; j++) {
    dct_col_odds[j] = extract_combine_odds(dct_cols[j + 0], dct_cols[j + 4]);
  }
  for (uint32_t j = 0; j < width; j++) {
    __m256i col = tsrc[j];
    __m256i odds = extract_odds(col);

    __m256i o04   = _mm256_madd_epi16           (odds,     dct_col_odds[0]);
    __m256i o15   = _mm256_madd_epi16           (odds,     dct_col_odds[1]);
    __m256i o26   = _mm256_madd_epi16           (odds,     dct_col_odds[2]);
    __m256i o37   = _mm256_madd_epi16           (odds,     dct_col_odds[3]);

    __m256i o0145 = _mm256_hadd_epi32           (o04,      o15);
    __m256i o2367 = _mm256_hadd_epi32           (o26,      o37);

    __m256i o     = _mm256_hadd_epi32           (o0145,    o2367);

    // D0,2 D0,6 D1,2 D1,6 D1,a D1,e D0,a D0,e | D2,2 D2,6 D3,2 D3,6 D3,a D3,e D2,a D2,e
    __m256i d_db2 = extract_26ae(dct_cols);

    // 2 6 2 6 a e a e | 2 6 2 6 a e a e
    __m256i t_db2 = extract_26ae_vec            (col);

    __m256i eo_parts  = _mm256_madd_epi16       (d_db2,    t_db2);
    __m256i eo_parts2 = _mm256_shuffle_epi32    (eo_parts, _MM_SHUFFLE(0, 1, 2, 3));

    // EO0 EO1 EO1 EO0 | EO2 EO3 EO3 EO2
    __m256i eo        = _mm256_add_epi32        (eo_parts, eo_parts2);
    __m256i eo2       = _mm256_permute4x64_epi64(eo,       _MM_SHUFFLE(1, 3, 2, 0));
    __m256i eo3       = _mm256_sign_epi32       (eo2,      eo_signmask);

    __m256i d_db4     = extract_d048c           (dct_cols);
    __m256i t_db4     = extract_d048c_vec       (col);
    __m256i eee_eeo   = _mm256_madd_epi16       (d_db4,   t_db4);

    __m256i eee_eee   = _mm256_permute4x64_epi64(eee_eeo,  _MM_SHUFFLE(3, 0, 3, 0));
    __m256i eeo_eeo1  = _mm256_permute4x64_epi64(eee_eeo,  _MM_SHUFFLE(1, 2, 1, 2));

    __m256i eeo_eeo2  = _mm256_sign_epi32       (eeo_eeo1, eeo_signmask);

    // EE0 EE1 EE2 EE3 | EE3 EE2 EE1 EE0
    __m256i ee        = _mm256_add_epi32        (eee_eee,  eeo_eeo2);
    __m256i e         = _mm256_add_epi32        (ee,       eo3);

    __m256i o_neg     = _mm256_sign_epi32       (o,        o_signmask);
    __m256i o_lo      = _mm256_blend_epi32      (o,        o_neg, 0xf0); // 1111 0000
    __m256i o_hi      = _mm256_blend_epi32      (o,        o_neg, 0x0f); // 0000 1111

    __m256i res_lo    = _mm256_add_epi32        (e,        o_lo);
    __m256i res_hi    = _mm256_add_epi32        (e,        o_hi);
    __m256i res_hi2   = _mm256_permute4x64_epi64(res_hi,   _MM_SHUFFLE(1, 0, 3, 2));

    __m256i res_lo_t  = truncate_inv(res_lo,  shift);
    __m256i res_hi_t  = truncate_inv(res_hi2, shift);

    __m256i res_16_1  = _mm256_packs_epi32      (res_lo_t, res_hi_t);
    __m256i final     = _mm256_shuffle_epi8     (res_16_1, final_shufmask);

    _mm256_store_si256((__m256i *)dst + j, final);
  }
}

static void matrix_idct_16x16_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);
  ALIGNED(64) int16_t tmp[16 * 16];

  partial_butterfly_inverse_16_avx2(input, tmp,    shift_1st);
  partial_butterfly_inverse_16_avx2(tmp,   output, shift_2nd);
}

static void matrix_dct_16x16_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = uvg_g_convert_to_bit[16] + 1 + (bitdepth - 8);
  int32_t shift_2nd = uvg_g_convert_to_bit[16] + 8;

  const int16_t *dct  = &uvg_g_dct_16[0][0];

  /*
   * Multiply input by the tranpose of DCT matrix into tmpres, and DCT matrix
   * by tmpres - this is then our output matrix
   *
   * It's easier to implement an AVX2 matrix multiplication if you can multiply
   * the left term with the transpose of the right term. Here things are stored
   * row-wise, not column-wise, so we can effectively read DCT_T column-wise
   * into YMM registers by reading DCT row-wise. Also because of this, the
   * first multiplication is hacked to produce the transpose of the result
   * instead, since it will be used in similar fashion as the right operand
   * in the second multiplication.
   */

  const __m256i *d_v = (const __m256i *)dct;
  const __m256i *i_v = (const __m256i *)input;
        __m256i *o_v = (      __m256i *)output;
  __m256i tmp[16];

  // Hack! (A * B^T)^T = B * A^T, so we can dispatch the transpose-produciong
  // multiply completely
  matmul_16x16_a_bt(d_v, i_v, tmp, shift_1st);
  matmul_16x16_a_bt(d_v, tmp, o_v, shift_2nd);
}

// 32x32 matrix multiplication with value clipping.
// Parameters: Two 32x32 matrices containing 16-bit values in consecutive addresses,
//             destination for the result and the shift value for clipping.
static void mul_clip_matrix_32x32_avx2(const int16_t *left,
                                       const int16_t *right,
                                             int16_t *dst,
                                       const int32_t  shift)
{
  const int32_t add    = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const uint32_t *l_32  = (const uint32_t *)left;
  const __m256i  *r_v   = (const __m256i *)right;
        __m256i  *dst_v = (      __m256i *)dst;

  __m256i accu[128] = {_mm256_setzero_si256()};
  size_t i, j;

  for (j = 0; j < 64; j += 4) {
    const __m256i r0 = r_v[j + 0];
    const __m256i r1 = r_v[j + 1];
    const __m256i r2 = r_v[j + 2];
    const __m256i r3 = r_v[j + 3];

    __m256i r02l   = _mm256_unpacklo_epi16(r0, r2);
    __m256i r02h   = _mm256_unpackhi_epi16(r0, r2);
    __m256i r13l   = _mm256_unpacklo_epi16(r1, r3);
    __m256i r13h   = _mm256_unpackhi_epi16(r1, r3);

    __m256i r02_07 = _mm256_permute2x128_si256(r02l, r02h, 0x20);
    __m256i r02_8f = _mm256_permute2x128_si256(r02l, r02h, 0x31);

    __m256i r13_07 = _mm256_permute2x128_si256(r13l, r13h, 0x20);
    __m256i r13_8f = _mm256_permute2x128_si256(r13l, r13h, 0x31);

    for (i = 0; i < 32; i += 2) {
      size_t acc_base = i << 2;

      uint32_t curr_e    = l_32[(i + 0) * (32 / 2) + (j >> 2)];
      uint32_t curr_o    = l_32[(i + 1) * (32 / 2) + (j >> 2)];

      __m256i even       = _mm256_set1_epi32(curr_e);
      __m256i odd        = _mm256_set1_epi32(curr_o);

      __m256i p_e0       = _mm256_madd_epi16(even, r02_07);
      __m256i p_e1       = _mm256_madd_epi16(even, r02_8f);
      __m256i p_e2       = _mm256_madd_epi16(even, r13_07);
      __m256i p_e3       = _mm256_madd_epi16(even, r13_8f);

      __m256i p_o0       = _mm256_madd_epi16(odd,  r02_07);
      __m256i p_o1       = _mm256_madd_epi16(odd,  r02_8f);
      __m256i p_o2       = _mm256_madd_epi16(odd,  r13_07);
      __m256i p_o3       = _mm256_madd_epi16(odd,  r13_8f);

      accu[acc_base + 0] = _mm256_add_epi32 (p_e0, accu[acc_base + 0]);
      accu[acc_base + 1] = _mm256_add_epi32 (p_e1, accu[acc_base + 1]);
      accu[acc_base + 2] = _mm256_add_epi32 (p_e2, accu[acc_base + 2]);
      accu[acc_base + 3] = _mm256_add_epi32 (p_e3, accu[acc_base + 3]);

      accu[acc_base + 4] = _mm256_add_epi32 (p_o0, accu[acc_base + 4]);
      accu[acc_base + 5] = _mm256_add_epi32 (p_o1, accu[acc_base + 5]);
      accu[acc_base + 6] = _mm256_add_epi32 (p_o2, accu[acc_base + 6]);
      accu[acc_base + 7] = _mm256_add_epi32 (p_o3, accu[acc_base + 7]);
    }
  }

  for (i = 0; i < 32; i++) {
    size_t acc_base = i << 2;
    size_t dst_base = i << 1;

    __m256i q0  = truncate_avx2(accu[acc_base + 0], debias, shift);
    __m256i q1  = truncate_avx2(accu[acc_base + 1], debias, shift);
    __m256i q2  = truncate_avx2(accu[acc_base + 2], debias, shift);
    __m256i q3  = truncate_avx2(accu[acc_base + 3], debias, shift);

    __m256i h01 = _mm256_packs_epi32(q0, q1);
    __m256i h23 = _mm256_packs_epi32(q2, q3);

            h01 = _mm256_permute4x64_epi64(h01, _MM_SHUFFLE(3, 1, 2, 0));
            h23 = _mm256_permute4x64_epi64(h23, _MM_SHUFFLE(3, 1, 2, 0));

    _mm256_store_si256(dst_v + dst_base + 0, h01);
    _mm256_store_si256(dst_v + dst_base + 1, h23);
  }
}

// Macro that generates 2D transform functions with clipping values.
// Sets correct shift values and matrices according to transform type and
// block size. Performs matrix multiplication horizontally and vertically.
#define TRANSFORM(type, n) static void matrix_ ## type ## _ ## n ## x ## n ## _avx2(int8_t bitdepth, const int16_t *input, int16_t *output)\
{\
  int32_t shift_1st = uvg_g_convert_to_bit[n] + 1 + (bitdepth - 8); \
  int32_t shift_2nd = uvg_g_convert_to_bit[n] + 8; \
  ALIGNED(64) int16_t tmp[n * n];\
  const int16_t *tdct = &uvg_g_ ## type ## _ ## n ## _t[0][0];\
  const int16_t *dct = &uvg_g_ ## type ## _ ## n [0][0];\
\
  mul_clip_matrix_ ## n ## x ## n ## _avx2(input, tdct, tmp, shift_1st);\
  mul_clip_matrix_ ## n ## x ## n ## _avx2(dct, tmp, output, shift_2nd);\
}\

// Macro that generates 2D inverse transform functions with clipping values.
// Sets correct shift values and matrices according to transform type and
// block size. Performs matrix multiplication horizontally and vertically.
#define ITRANSFORM(type, n) \
static void matrix_i ## type ## _## n ## x ## n ## _avx2(int8_t bitdepth, const int16_t *input, int16_t *output)\
{\
  int32_t shift_1st = 7; \
  int32_t shift_2nd = 12 - (bitdepth - 8); \
  ALIGNED(64) int16_t tmp[n * n];\
  const int16_t *tdct = &uvg_g_ ## type ## _ ## n ## _t[0][0];\
  const int16_t *dct = &uvg_g_ ## type ## _ ## n [0][0];\
\
  mul_clip_matrix_ ## n ## x ## n ## _avx2(tdct, input, tmp, shift_1st);\
  mul_clip_matrix_ ## n ## x ## n ## _avx2(tmp, dct, output, shift_2nd);\
}\

// Ha, we've got a tailored implementation for these
// TRANSFORM(dst, 4);
// ITRANSFORM(dst, 4);
// TRANSFORM(dct, 4);
// ITRANSFORM(dct, 4);
// TRANSFORM(dct, 8);
// ITRANSFORM(dct, 8);
// TRANSFORM(dct, 16);
// ITRANSFORM(dct, 16);

// Generate all the transform functions

TRANSFORM(dct, 32);
ITRANSFORM(dct, 32);


/*****************************************************/
/********************** M T S ************************/
/*****************************************************/

// DST-7

#define DEFINE_DST7_P4_MATRIX_T(a,b,c,d) { \
    { a,  c,  d,  b},\
    { b,  c, -a, -d},\
    { c,  0, -c,  c},\
    { d, -c,  b, -a},\
}


#define DEFINE_DST7_P8_MATRIX_T(a,b,c,d,e,f,g,h) \
{\
   { a,  c,  e,  g,  h,  f,  d,  b,},\
   { b,  f,  g,  c, -a, -e, -h, -d,},\
   { c,  h,  b, -d, -g, -a,  e,  f,},\
   { d,  e, -c, -f,  b,  g, -a, -h,},\
   { e,  b, -h,  a,  f, -d, -c,  g,},\
   { f, -a, -d,  h, -c, -b,  g, -e,},\
   { g, -d,  a,  b, -e,  h, -f,  c,},\
   { h, -g,  f, -e,  d, -c,  b, -a,},\
}\


#define DEFINE_DST7_P16_MATRIX_T(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) \
{ \
    {a,  c,  e,  g,  i,  k,  m,  o,  p,  n,  l,  j,  h,  f,  d,  b,},\
    {b,  f,  j,  n,  o,  k,  g,  c, -a, -e, -i, -m, -p, -l, -h, -d,},\
    {c,  i,  o,  l,  f,  0, -f, -l, -o, -i, -c,  c,  i,  o,  l,  f,},\
    {d,  l,  m,  e, -c, -k, -n, -f,  b,  j,  o,  g, -a, -i, -p, -h,},\
    {e,  o,  h, -b, -l, -k, -a,  i,  n,  d, -f, -p, -g,  c,  m,  j,},\
    {f,  o,  c, -i, -l,  0,  l,  i, -c, -o, -f,  f,  o,  c, -i, -l,},\
    {g,  l, -b, -p, -c,  k,  h, -f, -m,  a,  o,  d, -j, -i,  e,  n,},\
    {h,  i, -g, -j,  f,  k, -e, -l,  d,  m, -c, -n,  b,  o, -a, -p,},\
    {i,  f, -l, -c,  o,  0, -o,  c,  l, -f, -i,  i,  f, -l, -c,  o,},\
    {j,  c, -p,  d,  i, -k, -b,  o, -e, -h,  l,  a, -n,  f,  g, -m,},\
    {k,  0, -k,  k,  0, -k,  k,  0, -k,  k,  0, -k,  k,  0, -k,  k,},\
    {l, -c, -f,  o, -i,  0,  i, -o,  f,  c, -l,  l, -c, -f,  o, -i,},\
    {m, -f, -a,  h, -o,  k, -d, -c,  j, -p,  i, -b, -e,  l, -n,  g,},\
    {n, -i,  d,  a, -f,  k, -p,  l, -g,  b,  c, -h,  m, -o,  j, -e,},\
    {o, -l,  i, -f,  c,  0, -c,  f, -i,  l, -o,  o, -l,  i, -f,  c,},\
    {p, -o,  n, -m,  l, -k,  j, -i,  h, -g,  f, -e,  d, -c,  b, -a,},\
}


#define DEFINE_DST7_P32_MATRIX_T(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,A,B,C,D,E,F) \
{ \
    {a,  c,  e,  g,  i,  k,  m,  o,  q,  s,  u,  w,  y,  A,  C,  E,  F,  D,  B,  z,  x,  v,  t,  r,  p,  n,  l,  j,  h,  f,  d,  b,},\
    {b,  f,  j,  n,  r,  v,  z,  D,  E,  A,  w,  s,  o,  k,  g,  c, -a, -e, -i, -m, -q, -u, -y, -C, -F, -B, -x, -t, -p, -l, -h, -d,},\
    {c,  i,  o,  u,  A,  F,  z,  t,  n,  h,  b, -d, -j, -p, -v, -B, -E, -y, -s, -m, -g, -a,  e,  k,  q,  w,  C,  D,  x,  r,  l,  f,},\
    {d,  l,  t,  B,  C,  u,  m,  e, -c, -k, -s, -A, -D, -v, -n, -f,  b,  j,  r,  z,  E,  w,  o,  g, -a, -i, -q, -y, -F, -x, -p, -h,},\
    {e,  o,  y,  D,  t,  j,  0, -j, -t, -D, -y, -o, -e,  e,  o,  y,  D,  t,  j,  0, -j, -t, -D, -y, -o, -e,  e,  o,  y,  D,  t,  j,},\
    {f,  r,  D,  w,  k, -a, -m, -y, -B, -p, -d,  h,  t,  F,  u,  i, -c, -o, -A, -z, -n, -b,  j,  v,  E,  s,  g, -e, -q, -C, -x, -l,},\
    {g,  u,  D,  p,  b, -l, -z, -y, -k,  c,  q,  E,  t,  f, -h, -v, -C, -o, -a,  m,  A,  x,  j, -d, -r, -F, -s, -e,  i,  w,  B,  n,},\
    {h,  x,  y,  i, -g, -w, -z, -j,  f,  v,  A,  k, -e, -u, -B, -l,  d,  t,  C,  m, -c, -s, -D, -n,  b,  r,  E,  o, -a, -q, -F, -p,},\
    {i,  A,  t,  b, -p, -E, -m,  e,  w,  x,  f, -l, -D, -q,  a,  s,  B,  j, -h, -z, -u, -c,  o,  F,  n, -d, -v, -y, -g,  k,  C,  r,},\
    {j,  D,  o, -e, -y, -t,  0,  t,  y,  e, -o, -D, -j,  j,  D,  o, -e, -y, -t,  0,  t,  y,  e, -o, -D, -j,  j,  D,  o, -e, -y, -t,},\
    {k,  F,  j, -l, -E, -i,  m,  D,  h, -n, -C, -g,  o,  B,  f, -p, -A, -e,  q,  z,  d, -r, -y, -c,  s,  x,  b, -t, -w, -a,  u,  v,},\
    {l,  C,  e, -s, -v,  b,  z,  o, -i, -F, -h,  p,  y,  a, -w, -r,  f,  D,  k, -m, -B, -d,  t,  u, -c, -A, -n,  j,  E,  g, -q, -x,},\
    {m,  z,  0, -z, -m,  m,  z,  0, -z, -m,  m,  z,  0, -z, -m,  m,  z,  0, -z, -m,  m,  z,  0, -z, -m,  m,  z,  0, -z, -m,  m,  z,},\
    {n,  w, -e, -F, -d,  x,  m, -o, -v,  f,  E,  c, -y, -l,  p,  u, -g, -D, -b,  z,  k, -q, -t,  h,  C,  a, -A, -j,  r,  s, -i, -B,},\
    {o,  t, -j, -y,  e,  D,  0, -D, -e,  y,  j, -t, -o,  o,  t, -j, -y,  e,  D,  0, -D, -e,  y,  j, -t, -o,  o,  t, -j, -y,  e,  D,},\
    {p,  q, -o, -r,  n,  s, -m, -t,  l,  u, -k, -v,  j,  w, -i, -x,  h,  y, -g, -z,  f,  A, -e, -B,  d,  C, -c, -D,  b,  E, -a, -F,},\
    {q,  n, -t, -k,  w,  h, -z, -e,  C,  b, -F,  a,  D, -d, -A,  g,  x, -j, -u,  m,  r, -p, -o,  s,  l, -v, -i,  y,  f, -B, -c,  E,},\
    {r,  k, -y, -d,  F, -c, -z,  j,  s, -q, -l,  x,  e, -E,  b,  A, -i, -t,  p,  m, -w, -f,  D, -a, -B,  h,  u, -o, -n,  v,  g, -C,},\
    {s,  h, -D,  c,  x, -n, -m,  y,  b, -C,  i,  r, -t, -g,  E, -d, -w,  o,  l, -z, -a,  B, -j, -q,  u,  f, -F,  e,  v, -p, -k,  A,},\
    {t,  e, -D,  j,  o, -y,  0,  y, -o, -j,  D, -e, -t,  t,  e, -D,  j,  o, -y,  0,  y, -o, -j,  D, -e, -t,  t,  e, -D,  j,  o, -y,},\
    {u,  b, -y,  q,  f, -C,  m,  j, -F,  i,  n, -B,  e,  r, -x,  a,  v, -t, -c,  z, -p, -g,  D, -l, -k,  E, -h, -o,  A, -d, -s,  w,},\
    {v, -a, -t,  x, -c, -r,  z, -e, -p,  B, -g, -n,  D, -i, -l,  F, -k, -j,  E, -m, -h,  C, -o, -f,  A, -q, -d,  y, -s, -b,  w, -u,},\
    {w, -d, -o,  E, -l, -g,  z, -t,  a,  r, -B,  i,  j, -C,  q,  b, -u,  y, -f, -m,  F, -n, -e,  x, -v,  c,  p, -D,  k,  h, -A,  s,},\
    {x, -g, -j,  A, -u,  d,  m, -D,  r, -a, -p,  F, -o, -b,  s, -C,  l,  e, -v,  z, -i, -h,  y, -w,  f,  k, -B,  t, -c, -n,  E, -q,},\
    {y, -j, -e,  t, -D,  o,  0, -o,  D, -t,  e,  j, -y,  y, -j, -e,  t, -D,  o,  0, -o,  D, -t,  e,  j, -y,  y, -j, -e,  t, -D,  o,},\
    {z, -m,  0,  m, -z,  z, -m,  0,  m, -z,  z, -m,  0,  m, -z,  z, -m,  0,  m, -z,  z, -m,  0,  m, -z,  z, -m,  0,  m, -z,  z, -m,},\
    {A, -p,  e,  f, -q,  B, -z,  o, -d, -g,  r, -C,  y, -n,  c,  h, -s,  D, -x,  m, -b, -i,  t, -E,  w, -l,  a,  j, -u,  F, -v,  k,},\
    {B, -s,  j, -a, -h,  q, -z,  D, -u,  l, -c, -f,  o, -x,  F, -w,  n, -e, -d,  m, -v,  E, -y,  p, -g, -b,  k, -t,  C, -A,  r, -i,},\
    {C, -v,  o, -h,  a,  f, -m,  t, -A,  E, -x,  q, -j,  c,  d, -k,  r, -y,  F, -z,  s, -l,  e,  b, -i,  p, -w,  D, -B,  u, -n,  g,},\
    {D, -y,  t, -o,  j, -e,  0,  e, -j,  o, -t,  y, -D,  D, -y,  t, -o,  j, -e,  0,  e, -j,  o, -t,  y, -D,  D, -y,  t, -o,  j, -e,},\
    {E, -B,  y, -v,  s, -p,  m, -j,  g, -d,  a,  b, -e,  h, -k,  n, -q,  t, -w,  z, -C,  F, -D,  A, -x,  u, -r,  o, -l,  i, -f,  c,},\
    {F, -E,  D, -C,  B, -A,  z, -y,  x, -w,  v, -u,  t, -s,  r, -q,  p, -o,  n, -m,  l, -k,  j, -i,  h, -g,  f, -e,  d, -c,  b, -a,},\
}

// DST-7
ALIGNED(64) const int16_t uvg_g_dst7_4[4][4] = DEFINE_DST7_P4_MATRIX(29, 55, 74, 84);
ALIGNED(64) const int16_t uvg_g_dst7_8[8][8] = DEFINE_DST7_P8_MATRIX(17, 32, 46, 60, 71, 78, 85, 86);
ALIGNED(64) const int16_t uvg_g_dst7_16[16][16] = DEFINE_DST7_P16_MATRIX(8, 17, 25, 33, 40, 48, 55, 62, 68, 73, 77, 81, 85, 87, 88, 88);
ALIGNED(64) const int16_t uvg_g_dst7_32[32][32] = DEFINE_DST7_P32_MATRIX(4, 9, 13, 17, 21, 26, 30, 34, 38, 42, 46, 50, 53, 56, 60, 63, 66, 68, 72, 74, 77, 78, 80, 82, 84, 85, 86, 87, 88, 89, 90, 90);

ALIGNED(64) const int16_t uvg_g_dst7_4_t[4][4] = DEFINE_DST7_P4_MATRIX_T(29, 55, 74, 84);
ALIGNED(64) const int16_t uvg_g_dst7_8_t[8][8] = DEFINE_DST7_P8_MATRIX_T(17, 32, 46, 60, 71, 78, 85, 86);
ALIGNED(64) const int16_t uvg_g_dst7_16_t[16][16] = DEFINE_DST7_P16_MATRIX_T(8, 17, 25, 33, 40, 48, 55, 62, 68, 73, 77, 81, 85, 87, 88, 88);
ALIGNED(64) const int16_t uvg_g_dst7_32_t[32][32] = DEFINE_DST7_P32_MATRIX_T(4, 9, 13, 17, 21, 26, 30, 34, 38, 42, 46, 50, 53, 56, 60, 63, 66, 68, 72, 74, 77, 78, 80, 82, 84, 85, 86, 87, 88, 89, 90, 90);

// DCT-8
ALIGNED(64) const int16_t uvg_g_dct8_4[4][4] = DEFINE_DCT8_P4_MATRIX(84, 74, 55, 29);
ALIGNED(64) const int16_t uvg_g_dct8_8[8][8] = DEFINE_DCT8_P8_MATRIX(86, 85, 78, 71, 60, 46, 32, 17);
ALIGNED(64) const int16_t uvg_g_dct8_16[16][16] = DEFINE_DCT8_P16_MATRIX(88, 88, 87, 85, 81, 77, 73, 68, 62, 55, 48, 40, 33, 25, 17, 8);
ALIGNED(64) const int16_t uvg_g_dct8_32[32][32] = DEFINE_DCT8_P32_MATRIX(90, 90, 89, 88, 87, 86, 85, 84, 82, 80, 78, 77, 74, 72, 68, 66, 63, 60, 56, 53, 50, 46, 42, 38, 34, 30, 26, 21, 17, 13, 9, 4);

const int16_t* uvg_g_mts_input[2][3][5] = {
  {
    {&uvg_g_dct_4[0][0],  &uvg_g_dct_8[0][0],  &uvg_g_dct_16[0][0],  &uvg_g_dct_32[0][0], NULL},
    {&uvg_g_dct8_4[0][0], &uvg_g_dct8_8[0][0], &uvg_g_dct8_16[0][0], &uvg_g_dct8_32[0][0], NULL},
    {&uvg_g_dst7_4[0][0], &uvg_g_dst7_8[0][0], &uvg_g_dst7_16[0][0], &uvg_g_dst7_32[0][0], NULL}
  },
  {
    {&uvg_g_dct_4_t[0][0],  &uvg_g_dct_8_t[0][0],  &uvg_g_dct_16_t[0][0],  &uvg_g_dct_32_t[0][0], NULL},
    {  &uvg_g_dct8_4[0][0],   &uvg_g_dct8_8[0][0],   &uvg_g_dct8_16[0][0],   &uvg_g_dct8_32[0][0], NULL},
    {&uvg_g_dst7_4_t[0][0], &uvg_g_dst7_8_t[0][0], &uvg_g_dst7_16_t[0][0], &uvg_g_dst7_32_t[0][0], NULL}
  },
};

static void mts_dct_4x4_avx2(const int16_t* input, int16_t* output, tr_type_t type_hor, tr_type_t type_ver, uint8_t bitdepth, uint8_t lfnst_idx)
{
  //const int height = 4;
  const int width = 4;

  const int log2_width_minus2 = uvg_g_convert_to_bit[width];

  const int32_t shift_1st = log2_width_minus2 + bitdepth - 7;
  const int32_t shift_2nd = log2_width_minus2 + 8;

  const int16_t* tdct = uvg_g_mts_input[1][type_hor][0];
  const int16_t* dct = uvg_g_mts_input[0][type_ver][0];

  __m256i tdct_v = _mm256_load_si256((const __m256i*) tdct);
  __m256i  dct_v = _mm256_load_si256((const __m256i*)  dct);
  __m256i   in_v = _mm256_load_si256((const __m256i*)input);

  __m256i tmp = mul_clip_matrix_4x4_avx2(in_v, tdct_v, shift_1st);
  __m256i result = mul_clip_matrix_4x4_avx2(dct_v, tmp, shift_2nd);

  _mm256_store_si256((__m256i*)output, result);
}

static void mts_idct_4x4_avx2(const int16_t* input, int16_t* output, tr_type_t type_hor, tr_type_t type_ver, uint8_t bitdepth, uint8_t lfnst_idx)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);

  const int16_t* tdct = uvg_g_mts_input[1][type_ver][0];
  const int16_t* dct = uvg_g_mts_input[0][type_hor][0];

  __m256i tdct_v = _mm256_load_si256((const __m256i*)tdct);
  __m256i  dct_v = _mm256_load_si256((const __m256i*) dct);
  __m256i   in_v = _mm256_load_si256((const __m256i*)input);

  __m256i tmp = mul_clip_matrix_4x4_avx2(tdct_v, in_v, shift_1st);
  __m256i result = mul_clip_matrix_4x4_avx2(tmp, dct_v, shift_2nd);

  _mm256_store_si256((__m256i*)output, result);
}

static void mts_dct_8x8_avx2(const int16_t* input, int16_t* output, tr_type_t type_hor, tr_type_t type_ver, uint8_t bitdepth, uint8_t lfnst_idx)
{
  int32_t shift_1st = uvg_g_convert_to_bit[8] + 1 + (bitdepth - 8);
  int32_t shift_2nd = uvg_g_convert_to_bit[8] + 8;

  const int16_t* dct1 = uvg_g_mts_input[0][type_hor][1];
  const int16_t* dct2 = uvg_g_mts_input[0][type_ver][1];

  __m256i tmpres[4];

  matmul_8x8_a_bt_t(input, dct1, tmpres, shift_1st);
  matmul_8x8_a_bt(dct2, tmpres, output, shift_2nd);
}

static void mts_idct_8x8_avx2(const int16_t* input, int16_t* output, tr_type_t type_hor, tr_type_t type_ver, uint8_t bitdepth, uint8_t lfnst_idx)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);
  ALIGNED(64) int16_t tmp[8 * 8];

  const int16_t* tdct = uvg_g_mts_input[1][type_ver][1];
  const int16_t* dct = uvg_g_mts_input[0][type_hor][1];

  mul_clip_matrix_8x8_avx2(tdct, input, tmp, shift_1st);
  mul_clip_matrix_8x8_avx2(tmp, dct, output, shift_2nd);
}


static void mts_dct_16x16_avx2(const int16_t* input, int16_t* output, tr_type_t type_hor, tr_type_t type_ver, uint8_t bitdepth, uint8_t lfnst_idx)
{
  int32_t shift_1st = uvg_g_convert_to_bit[16] + 1 + (bitdepth - 8);
  int32_t shift_2nd = uvg_g_convert_to_bit[16] + 8;

  const int16_t* dct1 = uvg_g_mts_input[0][type_hor][2];
  const int16_t* dct2 = uvg_g_mts_input[0][type_ver][2];

  /*
   * Multiply input by the tranpose of DCT matrix into tmpres, and DCT matrix
   * by tmpres - this is then our output matrix
   *
   * It's easier to implement an AVX2 matrix multiplication if you can multiply
   * the left term with the transpose of the right term. Here things are stored
   * row-wise, not column-wise, so we can effectively read DCT_T column-wise
   * into YMM registers by reading DCT row-wise. Also because of this, the
   * first multiplication is hacked to produce the transpose of the result
   * instead, since it will be used in similar fashion as the right operand
   * in the second multiplication.
   */

  const __m256i* d_v = (const __m256i*)dct1;
  const __m256i* d_v2 = (const __m256i*)dct2;
  const __m256i* i_v = (const __m256i*)input;
  __m256i* o_v = (__m256i*)output;
  __m256i tmp[16];

  // Hack! (A * B^T)^T = B * A^T, so we can dispatch the transpose-produciong
  // multiply completely
  matmul_16x16_a_bt(d_v, i_v, tmp, shift_1st);
  matmul_16x16_a_bt(d_v2, tmp, o_v, shift_2nd);

  const int skip_line = lfnst_idx ? 8 : 0;
  const int skip_line2 = lfnst_idx ? 8 : 0;

}

/**********/
//ToDo: This function is not optimised! (DCT8/DST7 inverse 16x16)
/***********/
static void partial_butterfly_inverse_16_mts_avx2(const int16_t* src, int16_t* dst, int32_t shift, tr_type_t type)
{
  int j, k;
  int32_t a[5], b[5], c[5], d[5], t;
  int32_t add = (shift > 0) ? (1 << (shift - 1)) : 0;

  const int16_t* iT = &uvg_g_dst7_16[0][0];

  const int  line = 16;

  if (type == DST7) {
    for (j = 0; j < 16; j++)
    {
      for (k = 0; k < 5; k++)
      {
        a[k] = src[k * line] + src[(10 - k) * line];
        b[k] = src[(11 + k) * line] + src[(10 - k) * line];
        c[k] = src[k * line] - src[(11 + k) * line];
        d[k] = src[k * line] + src[(11 + k) * line] - src[(10 - k) * line];
      }

      t = iT[10] * src[5 * line];

      dst[2] = (short)CLIP(-32768, 32767, (iT[2] * d[0] + iT[8] * d[1] + iT[14] * d[2] + iT[11] * d[3] + iT[5] * d[4] + add) >> shift);
      dst[5] = (short)CLIP(-32768, 32767, (iT[5] * d[0] + iT[14] * d[1] + iT[2] * d[2] - iT[8] * d[3] - iT[11] * d[4] + add) >> shift);
      dst[8] = (short)CLIP(-32768, 32767, (iT[8] * d[0] + iT[5] * d[1] - iT[11] * d[2] - iT[2] * d[3] + iT[14] * d[4] + add) >> shift);
      dst[11] = (short)CLIP(-32768, 32767, (iT[11] * d[0] - iT[2] * d[1] - iT[5] * d[2] + iT[14] * d[3] - iT[8] * d[4] + add) >> shift);
      dst[14] = (short)CLIP(-32768, 32767, (iT[14] * d[0] - iT[11] * d[1] + iT[8] * d[2] - iT[5] * d[3] + iT[2] * d[4] + add) >> shift);

      dst[10] = (short)CLIP(-32768, 32767, (iT[10] * (src[0 * line] - src[2 * line] + src[3 * line] - src[5 * line]
        + src[6 * line] - src[8 * line] + src[9 * line] - src[11 * line]
        + src[12 * line] - src[14 * line] + src[15 * line]) + add) >> shift);

      dst[0] = (short)CLIP(-32768, 32767, (iT[0] * a[0] + iT[9] * b[0] + iT[2] * a[1] + iT[7] * b[1] + iT[4] * a[2] + iT[5] * b[2] + iT[6] * a[3] + iT[3] * b[3] + iT[8] * a[4] + iT[1] * b[4] + t + add) >> shift);
      dst[1] = (short)CLIP(-32768, 32767, (iT[1] * c[0] - iT[8] * b[0] + iT[5] * c[1] - iT[4] * b[1] + iT[9] * c[2] - iT[0] * b[2] + iT[2] * a[3] + iT[7] * c[3] + iT[6] * a[4] + iT[3] * c[4] + t + add) >> shift);
      dst[3] = (short)CLIP(-32768, 32767, (iT[3] * a[0] + iT[6] * b[0] + iT[0] * c[1] + iT[9] * a[1] + iT[1] * a[2] + iT[8] * c[2] + iT[4] * c[3] - iT[5] * b[3] - iT[2] * a[4] - iT[7] * b[4] - t + add) >> shift);
      dst[4] = (short)CLIP(-32768, 32767, (iT[4] * c[0] - iT[5] * b[0] + iT[6] * c[1] + iT[3] * a[1] + iT[7] * a[2] + iT[2] * b[2] - iT[1] * c[3] + iT[8] * b[3] - iT[9] * c[4] - iT[0] * a[4] - t + add) >> shift);
      dst[6] = (short)CLIP(-32768, 32767, (iT[6] * a[0] + iT[3] * b[0] + iT[9] * c[1] + iT[0] * a[1] - iT[1] * a[2] - iT[8] * b[2] - iT[4] * c[3] - iT[5] * a[3] - iT[2] * c[4] + iT[7] * b[4] + t + add) >> shift);
      dst[7] = (short)CLIP(-32768, 32767, (iT[7] * c[0] - iT[2] * b[0] + iT[8] * a[1] + iT[1] * b[1] - iT[6] * c[2] + iT[3] * b[2] - iT[9] * a[3] - iT[0] * b[3] + iT[5] * c[4] - iT[4] * b[4] + t + add) >> shift);
      dst[9] = (short)CLIP(-32768, 32767, (iT[9] * a[0] + iT[0] * b[0] + iT[2] * c[1] - iT[7] * b[1] - iT[5] * c[2] - iT[4] * a[2] + iT[3] * a[3] + iT[6] * b[3] + iT[8] * c[4] - iT[1] * b[4] - t + add) >> shift);
      dst[12] = (short)CLIP(-32768, 32767, (iT[1] * c[0] + iT[8] * a[0] - iT[5] * a[1] - iT[4] * b[1] - iT[0] * c[2] + iT[9] * b[2] + iT[7] * c[3] - iT[2] * b[3] - iT[6] * c[4] - iT[3] * a[4] + t + add) >> shift);
      dst[13] = (short)CLIP(-32768, 32767, (iT[7] * c[0] + iT[2] * a[0] - iT[8] * c[1] + iT[1] * b[1] + iT[3] * c[2] - iT[6] * b[2] + iT[0] * a[3] + iT[9] * b[3] - iT[5] * a[4] - iT[4] * b[4] + t + add) >> shift);
      dst[15] = (short)CLIP(-32768, 32767, (iT[4] * c[0] + iT[5] * a[0] - iT[3] * c[1] - iT[6] * a[1] + iT[2] * c[2] + iT[7] * a[2] - iT[1] * c[3] - iT[8] * a[3] + iT[0] * c[4] + iT[9] * a[4] - t + add) >> shift);
      src++;
      dst += 16;
    }
  }
  else {
    for (j = 0; j < 16; j++)
    {
      for (k = 0; k < 5; k++)
      {
        a[k] = src[(15 - k) * line] + src[(4 - k) * line];
        b[k] = src[(6 + k) * line] + src[(4 - k) * line];
        c[k] = src[(15 - k) * line] - src[(6 + k) * line];
        d[k] = src[(15 - k) * line] + src[(6 + k) * line] - src[(4 - k) * line];
      }

      t = iT[10] * src[5 * line];

      dst[1] = (short)CLIP(-32768, 32767, (-iT[2] * d[0] - iT[5] * d[1] - iT[8] * d[2] - iT[11] * d[3] - iT[14] * d[4] + add) >> shift);
      dst[4] = (short)CLIP(-32768, 32767, (iT[8] * d[0] + iT[14] * d[1] + iT[5] * d[2] - iT[2] * d[3] - iT[11] * d[4] + add) >> shift);
      dst[7] = (short)CLIP(-32768, 32767, (-iT[14] * d[0] - iT[2] * d[1] + iT[11] * d[2] + iT[5] * d[3] - iT[8] * d[4] + add) >> shift);
      dst[10] = (short)CLIP(-32768, 32767, (iT[11] * d[0] - iT[8] * d[1] - iT[2] * d[2] + iT[14] * d[3] - iT[5] * d[4] + add) >> shift);
      dst[13] = (short)CLIP(-32768, 32767, (-iT[5] * d[0] + iT[11] * d[1] - iT[14] * d[2] + iT[8] * d[3] - iT[2] * d[4] + add) >> shift);

      dst[5] = (short)CLIP(-32768, 32767, (-iT[10] * (src[15 * line] + src[14 * line] - src[12 * line] - src[11 * line] + src[9 * line] + src[8 * line] - src[6 * line] - src[5 * line] + src[3 * line] + src[2 * line] - src[0 * line]) + add) >> shift);

      dst[0] = (short)CLIP(-32768, 32767, (iT[0] * a[0] + iT[9] * b[0] + iT[1] * a[1] + iT[8] * b[1] + iT[2] * a[2] + iT[7] * b[2] + iT[3] * a[3] + iT[6] * b[3] + iT[4] * a[4] + iT[5] * b[4] + t + add) >> shift);
      dst[2] = (short)CLIP(-32768, 32767, (iT[4] * c[0] - iT[5] * b[0] + iT[9] * c[1] - iT[0] * b[1] + iT[6] * c[2] + iT[3] * a[2] + iT[1] * c[3] + iT[8] * a[3] + iT[7] * a[4] + iT[2] * b[4] - t + add) >> shift);
      dst[3] = (short)CLIP(-32768, 32767, (-iT[6] * a[0] - iT[3] * b[0] - iT[2] * c[1] - iT[7] * a[1] - iT[9] * c[2] - iT[0] * a[2] - iT[4] * c[3] + iT[5] * b[3] + iT[1] * a[4] + iT[8] * b[4] - t + add) >> shift);
      dst[6] = (short)CLIP(-32768, 32767, (iT[8] * a[0] + iT[1] * c[0] + iT[6] * c[1] - iT[3] * b[1] - iT[5] * a[2] - iT[4] * b[2] - iT[7] * c[3] - iT[2] * a[3] - iT[0] * c[4] + iT[9] * b[4] + t + add) >> shift);
      dst[8] = (short)CLIP(-32768, 32767, (iT[4] * c[0] + iT[5] * a[0] - iT[0] * c[1] + iT[9] * b[1] - iT[3] * c[2] - iT[6] * a[2] + iT[1] * c[3] - iT[8] * b[3] + iT[2] * c[4] + iT[7] * a[4] - t + add) >> shift);
      dst[9] = (short)CLIP(-32768, 32767, (-iT[7] * c[0] - iT[2] * a[0] + iT[4] * a[1] + iT[5] * b[1] + iT[8] * c[2] - iT[1] * b[2] - iT[9] * a[3] - iT[0] * b[3] - iT[3] * c[4] + iT[6] * b[4] - t + add) >> shift);
      dst[11] = (short)CLIP(-32768, 32767, (-iT[9] * a[0] - iT[0] * b[0] + iT[8] * c[1] + iT[1] * a[1] - iT[2] * c[2] + iT[7] * b[2] - iT[6] * a[3] - iT[3] * b[3] + iT[5] * c[4] + iT[4] * a[4] + t + add) >> shift);
      dst[12] = (short)CLIP(-32768, 32767, (iT[7] * c[0] - iT[2] * b[0] - iT[5] * c[1] - iT[4] * a[1] + iT[8] * a[2] + iT[1] * b[2] - iT[0] * a[3] - iT[9] * b[3] - iT[6] * c[4] + iT[3] * b[4] + t + add) >> shift);
      dst[14] = (short)CLIP(-32768, 32767, (iT[3] * a[0] + iT[6] * b[0] - iT[7] * a[1] - iT[2] * b[1] + iT[0] * c[2] + iT[9] * a[2] - iT[4] * c[3] - iT[5] * a[3] + iT[8] * c[4] + iT[1] * a[4] - t + add) >> shift);
      dst[15] = (short)CLIP(-32768, 32767, (-iT[1] * c[0] + iT[8] * b[0] + iT[3] * c[1] - iT[6] * b[1] - iT[5] * c[2] + iT[4] * b[2] + iT[7] * c[3] - iT[2] * b[3] - iT[9] * c[4] + iT[0] * b[4] - t + add) >> shift);
      src++;
      dst += 16;
    }
  }
}

static void mts_idct_16x16_avx2(const int16_t* input, int16_t* output, tr_type_t type_hor, tr_type_t type_ver, uint8_t bitdepth, uint8_t lfnst_idx)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);
  ALIGNED(64) int16_t tmp[16 * 16];

  partial_butterfly_inverse_16_mts_avx2(input, tmp, shift_1st, type_ver);
  partial_butterfly_inverse_16_mts_avx2(tmp, output, shift_2nd, type_hor);
}

// 32x32 matrix multiplication with value clipping.
// Parameters: Two 32x32 matrices containing 16-bit values in consecutive addresses,
//             destination for the result and the shift value for clipping.
static void mul_clip_matrix_32x32_mts_avx2(const int16_t* left,
  const int16_t* right,
  int16_t* dst,
  const int32_t  shift, int skip_line, int skip_line2)
{
  const int32_t add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const int  reduced_line = 32 - skip_line;
  const int  cutoff = 32 - skip_line2;

  const uint32_t* l_32 = (const uint32_t*)left;
  const __m256i* r_v = (const __m256i*)right;
  __m256i* dst_v = (__m256i*)dst;

  __m256i accu[128] = { _mm256_setzero_si256() };
  size_t i, j;

  for (j = 0; j < 64; j += 4) {
    const __m256i r0 = r_v[j + 0];
    const __m256i r1 = r_v[j + 1];
    const __m256i r2 = r_v[j + 2];
    const __m256i r3 = r_v[j + 3];

    __m256i r02l = _mm256_unpacklo_epi16(r0, r2);
    __m256i r02h = _mm256_unpackhi_epi16(r0, r2);
    __m256i r13l = _mm256_unpacklo_epi16(r1, r3);
    __m256i r13h = _mm256_unpackhi_epi16(r1, r3);

    __m256i r02_07 = _mm256_permute2x128_si256(r02l, r02h, 0x20);
    __m256i r02_8f = _mm256_permute2x128_si256(r02l, r02h, 0x31);

    __m256i r13_07 = _mm256_permute2x128_si256(r13l, r13h, 0x20);
    __m256i r13_8f = _mm256_permute2x128_si256(r13l, r13h, 0x31);

    for (i = 0; i < 32; i += 2) {
      size_t acc_base = i << 2;

      uint32_t curr_e = l_32[(i + 0) * (32 / 2) + (j >> 2)];
      uint32_t curr_o = l_32[(i + 1) * (32 / 2) + (j >> 2)];

      __m256i even = _mm256_set1_epi32(curr_e);
      __m256i odd = _mm256_set1_epi32(curr_o);

      __m256i p_e0 = _mm256_madd_epi16(even, r02_07);
      __m256i p_e1 = _mm256_madd_epi16(even, r02_8f);
      __m256i p_e2 = _mm256_madd_epi16(even, r13_07);
      __m256i p_e3 = _mm256_madd_epi16(even, r13_8f);

      __m256i p_o0 = _mm256_madd_epi16(odd, r02_07);
      __m256i p_o1 = _mm256_madd_epi16(odd, r02_8f);
      __m256i p_o2 = _mm256_madd_epi16(odd, r13_07);
      __m256i p_o3 = _mm256_madd_epi16(odd, r13_8f);

      accu[acc_base + 0] = _mm256_add_epi32(p_e0, accu[acc_base + 0]);
      accu[acc_base + 1] = _mm256_add_epi32(p_e1, accu[acc_base + 1]);
      accu[acc_base + 2] = _mm256_add_epi32(p_e2, accu[acc_base + 2]);
      accu[acc_base + 3] = _mm256_add_epi32(p_e3, accu[acc_base + 3]);

      accu[acc_base + 4] = _mm256_add_epi32(p_o0, accu[acc_base + 4]);
      accu[acc_base + 5] = _mm256_add_epi32(p_o1, accu[acc_base + 5]);
      accu[acc_base + 6] = _mm256_add_epi32(p_o2, accu[acc_base + 6]);
      accu[acc_base + 7] = _mm256_add_epi32(p_o3, accu[acc_base + 7]);
    }
  }

  for (i = 0; i < 32; i++) {
    size_t acc_base = i << 2;
    size_t dst_base = i << 1;

    __m256i q0 = truncate_avx2(accu[acc_base + 0], debias, shift);
    __m256i q1 = truncate_avx2(accu[acc_base + 1], debias, shift);
    __m256i q2 = truncate_avx2(accu[acc_base + 2], debias, shift);
    __m256i q3 = truncate_avx2(accu[acc_base + 3], debias, shift);

    __m256i h01 = _mm256_packs_epi32(q0, q1);
    __m256i h23 = _mm256_packs_epi32(q2, q3);

    h01 = _mm256_permute4x64_epi64(h01, _MM_SHUFFLE(3, 1, 2, 0));
    h23 = _mm256_permute4x64_epi64(h23, _MM_SHUFFLE(3, 1, 2, 0));

    _mm256_store_si256(dst_v + dst_base + 0, h01);
    _mm256_store_si256(dst_v + dst_base + 1, h23);
  }
  

}

static void mts_dct_32x32_avx2(const int16_t* input, int16_t* output, tr_type_t type_hor, tr_type_t type_ver, uint8_t bitdepth, uint8_t lfnst_idx)
{
  int32_t shift_1st = uvg_g_convert_to_bit[32] + 1 + (bitdepth - 8); 
  int32_t shift_2nd = uvg_g_convert_to_bit[32] + 8; 
  ALIGNED(64) int16_t tmp[32 * 32];

  const int16_t* tdct = uvg_g_mts_input[1][type_hor][3];
  const int16_t* dct = uvg_g_mts_input[0][type_ver][3];

  int skip_width = (type_hor != DCT2) ? 16 : 0;
  int skip_height = (type_ver != DCT2) ? 16 : 0;
  if(lfnst_idx) {
    skip_width = 24;
    skip_height = 24;
  }

  mul_clip_matrix_32x32_mts_avx2(input, tdct, tmp, shift_1st, skip_width, 0 );
  mul_clip_matrix_32x32_mts_avx2(dct, tmp, output, shift_2nd, skip_width, skip_height);
}


static void mts_idct_32x32_avx2(const int16_t* input, int16_t* output, tr_type_t type_hor, tr_type_t type_ver, uint8_t bitdepth, uint8_t lfnst_idx)
{
  int32_t shift_1st = 7; 
  int32_t shift_2nd = 12 - (bitdepth - 8); 
  ALIGNED(64) int16_t tmp[32 * 32];
  const int16_t* tdct = uvg_g_mts_input[1][type_ver][3];
  const int16_t* dct = uvg_g_mts_input[0][type_hor][3];

  //const int skip_width = (type_hor != DCT2) ? 16 : 0;
  const int skip_height = (type_ver != DCT2) ? 16 : 0;

  mul_clip_matrix_32x32_mts_avx2(tdct, input, tmp, shift_1st, skip_height,0);
  mul_clip_matrix_32x32_mts_avx2(tmp, dct, output, shift_2nd, 0, 0);
}

typedef void tr_func(const int16_t*, int16_t*, tr_type_t , tr_type_t , uint8_t, uint8_t);

// ToDo: Enable MTS 2x2 and 64x64 transforms
static tr_func* dct_table[5] = {
    mts_dct_4x4_avx2, mts_dct_8x8_avx2, mts_dct_16x16_avx2, mts_dct_32x32_avx2, NULL
};

static tr_func* idct_table[5] = {
  mts_idct_4x4_avx2, mts_idct_8x8_avx2, mts_idct_16x16_avx2, mts_idct_32x32_avx2, NULL/*fastInverseDCT2_B64*/
};

typedef void (dct_full_pass)(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver);


// **********************************************
// New tailored functions for each size combination
// **********************************************

static void fast_forward_tr_2xN_avx2_hor(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = (shift > 0) ? (1 << (shift - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);
  __m256i    v_coeff_0 = _mm256_load_si256((__m256i*)coeff);
  __m256i    v_coeff_1 = _mm256_load_si256((__m256i*)(coeff + 16));
  __m256i* v_dst_ptr = dst;

  const int reduced_line = line - skip_line;
  // Handle 8 lines at a time (16 samples, 2 samples per line)
  for (int j = 0; j < reduced_line; j += 8) {
    // src vector: [00 01 02 03 04 05 06 07|08 09 10 11 12 13 14 15]
    __m256i     v_src = _mm256_load_si256((const __m256i*) src);

    // Multiply with a and add together all adjacent elements
    // even vector: [a00+a01 a02+a03 a04+a05 a06+a07|a08+a09 a10+a11 a12+a13 a14+a15]
    __m256i    v_even = _mm256_madd_epi16(v_src, v_coeff_0);
    // odd vector : [a00-a01 a02-a03 a04-a05 a06-a07|a08-a09 a10-a11 a12-a13 a14-a15]
    __m256i     v_odd = _mm256_madd_epi16(v_src, v_coeff_1);

    __m256i v_trunc_0 = truncate_avx2(v_even, debias, shift);
    __m256i v_trunc_1 = truncate_avx2(v_odd, debias, shift);

    v_dst_ptr[0] = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
    
    src += 16;
    v_dst_ptr++;
  }
}

static void fast_forward_tr_2x8_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 2;
  const int height = 8;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_2xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_2x8_coeff_ver;
  if (ver == DST7) {
    ver_coeff = ff_dst7_2x8_coeff_ver;
  }
  // No coeffs for DCT8 and DST7 transforms since they do not exist for this block size

  __m256i v_hor_pass_out;
  fast_forward_tr_2xN_avx2_hor(src, &v_hor_pass_out, hor_coeff, shift_1st, height, 0, 0);

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)ver_coeff;

  // Got data for only 1 vector
  // const __m256i v_shuffle = _mm256_load_si256((const __m256i*)ff_dct2_2x8_shuffle_ver);
  const __m256i v_src_raw = v_hor_pass_out;
  // __m256i           v_src = _mm256_shuffle_epi8(v_src_raw, v_shuffle);
  __m256i           v_src = _mm256_permute4x64_epi64(v_src_raw, _MM_SHUFFLE(3, 1, 2, 0));

  __m256i v_madd[8];
  for (int i = 0; i < 8; ++i) {
    v_madd[i] = _mm256_madd_epi16(v_src, v_coeff[i]);
  }
  __m256i v_hadd_0[4];
  for (int i = 0; i < 4; ++i) {
    const int offset = i * 2;
    v_hadd_0[i] = _mm256_hadd_epi32(v_madd[offset], v_madd[offset + 1]);
  }

  __m256i v_trunc[2];
  for (int i = 0; i < 2; ++i) {
    const int offset = i * 2;
    v_trunc[i] = truncate_avx2(_mm256_hadd_epi32(v_hadd_0[offset], v_hadd_0[offset + 1]), debias, shift_2nd);
  }

  __m256i v_result = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
  const __m256i v_res_shfl = _mm256_load_si256((const __m256i*)ff_dct2_2x8_result_shuffle_ver);
  // Shuffle values to correct order
  v_result = _mm256_permute4x64_epi64(v_result, _MM_SHUFFLE(3, 1, 2, 0));
  v_result = _mm256_shuffle_epi32(v_result, _MM_SHUFFLE(3, 1, 2, 0));
  v_result = _mm256_shuffle_epi8(v_result, v_res_shfl);
  _mm256_store_si256((__m256i*)dst, v_result);
}


static void fast_inverse_tr_2x8_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i v_shuffle = _mm256_load_si256((const __m256i*)fi_tr_8x2_shuffle_hor);

  const __m256i v_src_raw = _mm256_load_si256((const __m256i*)src);

  __m256i v_src = _mm256_shuffle_epi8(v_src_raw, v_shuffle);
  v_src = _mm256_permute4x64_epi64(v_src, _MM_SHUFFLE(3, 1, 2, 0));

  __m256i v_madd_0 = _mm256_madd_epi16(v_src, v_coeff[0]);
  __m256i v_madd_1 = _mm256_madd_epi16(v_src, v_coeff[1]);
  __m256i v_madd_2 = _mm256_madd_epi16(v_src, v_coeff[2]);
  __m256i v_madd_3 = _mm256_madd_epi16(v_src, v_coeff[3]);
  __m256i v_madd_4 = _mm256_madd_epi16(v_src, v_coeff[4]);
  __m256i v_madd_5 = _mm256_madd_epi16(v_src, v_coeff[5]);
  __m256i v_madd_6 = _mm256_madd_epi16(v_src, v_coeff[6]);
  __m256i v_madd_7 = _mm256_madd_epi16(v_src, v_coeff[7]);

  __m256i v_hadd_00 = _mm256_hadd_epi32(v_madd_0, v_madd_1);
  __m256i v_hadd_01 = _mm256_hadd_epi32(v_madd_2, v_madd_3);
  __m256i v_hadd_02 = _mm256_hadd_epi32(v_madd_4, v_madd_5);
  __m256i v_hadd_03 = _mm256_hadd_epi32(v_madd_6, v_madd_7);

  __m256i v_hadd_10 = _mm256_hadd_epi32(v_hadd_00, v_hadd_01);
  __m256i v_hadd_11 = _mm256_hadd_epi32(v_hadd_02, v_hadd_03);

  __m256i v_trunc_0 = truncate_avx2(v_hadd_10, debias, shift);
  __m256i v_trunc_1 = truncate_avx2(v_hadd_11, debias, shift);

  dst[0] = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
}

static void fast_inverse_tr_2x8_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i v_coeff_0 = _mm256_load_si256((const __m256i*) & coeff[0]);
  const __m256i v_coeff_1 = _mm256_load_si256((const __m256i*) & coeff[16]);
  const __m256i v_shuffle = _mm256_load_si256((const __m256i*)fi_tr_8x2_shuffle_ver);
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)fi_tr_8x2_res_shuffle_ver);

  __m256i v_src = _mm256_permute4x64_epi64(src[0], _MM_SHUFFLE(3, 1, 2, 0));
  v_src = _mm256_shuffle_epi8(v_src, v_shuffle);

  __m256i v_madd_0 = _mm256_madd_epi16(v_src, v_coeff_0);
  __m256i v_madd_1 = _mm256_madd_epi16(v_src, v_coeff_1);

  __m256i v_trunc_0 = truncate_avx2(v_madd_0, debias, shift);
  __m256i v_trunc_1 = truncate_avx2(v_madd_1, debias, shift);

  __m256i v_result = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
  v_result = _mm256_shuffle_epi8(v_result, v_res_shuffle);

  _mm256_store_si256((__m256i*)dst, v_result);
}

static void fast_inverse_tr_2x8_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 2;
  const int height = 8;

  int skip_width = 0;
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = fi_dct2_8x2_coeff_hor; // TODO: rename
  if (ver == DST7) {
    ver_coeff = fi_dst7_8x2_coeff_hor;
  }
  const int16_t* hor_coeff = fi_dct2_8x2_coeff_ver; // rename
  // No coeffs for DCT8 and DST7 transforms since they do not exist for this block size

  __m256i v_ver_pass_out;
  fast_inverse_tr_2x8_avx2_ver(src, &v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);
  
  fast_inverse_tr_2x8_avx2_hor(&v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static void fast_forward_tr_2x16_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 2;
  const int height = 16;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_2xN_coeff_hor;
  const int16_t* ver_coeff = &uvg_g_dct_16[0][0];
  if (ver == DST7) {
    ver_coeff = &uvg_g_dst7_16[0][0];
  }
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)ff_dct2_2x16_ver_result_shuffle);
  // No coeffs for DCT8 and DST7 transforms since they do not exist for this block size

  __m256i v_hor_pass_out[2];
  fast_forward_tr_2xN_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, 0);

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  // Permute hor pass output to correct order
  __m256i v_tmp_0 = _mm256_permute4x64_epi64(v_hor_pass_out[0], _MM_SHUFFLE(3, 1, 2, 0));
  __m256i v_tmp_1 = _mm256_permute4x64_epi64(v_hor_pass_out[1], _MM_SHUFFLE(3, 1, 2, 0));
  __m256i v_src_0 = _mm256_permute2x128_si256(v_tmp_0, v_tmp_1, 0x20);
  __m256i v_src_1 = _mm256_permute2x128_si256(v_tmp_0, v_tmp_1, 0x31);

  const __m256i* v_coeff_ptr = (const __m256i*)ver_coeff;

  __m256i v_madd[2][16];
  for (int i = 0; i < 16; ++i) {
    v_madd[0][i] = _mm256_madd_epi16(v_src_0, v_coeff_ptr[i]);
    v_madd[1][i] = _mm256_madd_epi16(v_src_1, v_coeff_ptr[i]);
  }

  __m256i v_hadd_0[2][8];
  for (int dst = 0, src = 0; dst < 8; ++dst, src += 2) {
    v_hadd_0[0][dst] = _mm256_hadd_epi32(v_madd[0][src], v_madd[0][src + 1]);
    v_hadd_0[1][dst] = _mm256_hadd_epi32(v_madd[1][src], v_madd[1][src + 1]);
  }

  __m256i v_hadd_1[2][4];
  for (int dst = 0, src = 0; dst < 4; ++dst, src += 2) {
    v_hadd_1[0][dst] = _mm256_hadd_epi32(v_hadd_0[0][src], v_hadd_0[0][src + 1]);
    v_hadd_1[1][dst] = _mm256_hadd_epi32(v_hadd_0[1][src], v_hadd_0[1][src + 1]);
  }

  __m256i v_tmp_00 = _mm256_permute2x128_si256(v_hadd_1[0][0], v_hadd_1[0][1], 0x20);
  __m256i v_tmp_01 = _mm256_permute2x128_si256(v_hadd_1[0][0], v_hadd_1[0][1], 0x31);
  __m256i v_tmp_02 = _mm256_permute2x128_si256(v_hadd_1[0][2], v_hadd_1[0][3], 0x20);
  __m256i v_tmp_03 = _mm256_permute2x128_si256(v_hadd_1[0][2], v_hadd_1[0][3], 0x31);

  __m256i v_tmp_10 = _mm256_permute2x128_si256(v_hadd_1[1][0], v_hadd_1[1][1], 0x20);
  __m256i v_tmp_11 = _mm256_permute2x128_si256(v_hadd_1[1][0], v_hadd_1[1][1], 0x31);
  __m256i v_tmp_12 = _mm256_permute2x128_si256(v_hadd_1[1][2], v_hadd_1[1][3], 0x20);
  __m256i v_tmp_13 = _mm256_permute2x128_si256(v_hadd_1[1][2], v_hadd_1[1][3], 0x31);

  __m256i v_trunc_00 = truncate_avx2((_mm256_add_epi32(v_tmp_00, v_tmp_01)), debias, shift_2nd);
  __m256i v_trunc_01 = truncate_avx2((_mm256_add_epi32(v_tmp_02, v_tmp_03)), debias, shift_2nd);

  __m256i v_trunc_10 = truncate_avx2((_mm256_add_epi32(v_tmp_10, v_tmp_11)), debias, shift_2nd);
  __m256i v_trunc_11 = truncate_avx2((_mm256_add_epi32(v_tmp_12, v_tmp_13)), debias, shift_2nd);

  __m256i v_result_0 = _mm256_packs_epi32(v_trunc_00, v_trunc_10);
  __m256i v_result_1 = _mm256_packs_epi32(v_trunc_01, v_trunc_11);

  v_result_0 = _mm256_shuffle_epi8(v_result_0, v_res_shuffle);
  v_result_1 = _mm256_shuffle_epi8(v_result_1, v_res_shuffle);

  _mm256_store_si256((__m256i*)&dst[0], v_result_0);
  _mm256_store_si256((__m256i*)&dst[16], v_result_1);
}


static void fast_inverse_tr_2x16_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i* v_src_raw = (const __m256i*)src;
  const __m256i v_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0246);

  __m256i v_src_0 = _mm256_shuffle_epi8(v_src_raw[0], v_shuffle);
  __m256i v_src_1 = _mm256_shuffle_epi8(v_src_raw[1], v_shuffle);

  v_src_0 = _mm256_permute4x64_epi64(v_src_0, _MM_SHUFFLE(3, 1, 2, 0));
  v_src_1 = _mm256_permute4x64_epi64(v_src_1, _MM_SHUFFLE(3, 1, 2, 0));

  __m256i v_madd_0[16];
  __m256i v_madd_1[16];
  for (int c = 0; c < 16; ++c) {
    v_madd_0[c] = _mm256_madd_epi16(v_src_0, v_coeff[0]);
    v_madd_1[c] = _mm256_madd_epi16(v_src_1, v_coeff[1]);
    v_coeff += 2;
  }

  __m256i v_add[16];
  for (int i = 0; i < 16; ++i) {
    v_add[i] = _mm256_add_epi32(v_madd_0[i], v_madd_1[i]);
  }

  __m256i v_hadd_0[8];
  for (int d = 0, s = 0; d < 8; ++d, s += 2) {
    v_hadd_0[d] = _mm256_hadd_epi32(v_add[s + 0], v_add[s + 1]);
  }

  __m256i v_hadd_1[4];
  for (int d = 0, s = 0; d < 4; ++d, s += 2) {
    v_hadd_1[d] = _mm256_hadd_epi32(v_hadd_0[s + 0], v_hadd_0[s + 1]);
  }

  __m256i v_trunc[4];
  for (int i = 0; i < 4; ++i) {
    v_trunc[i] = truncate_avx2(v_hadd_1[i], debias, shift);
  }

  __m256i v_tmp0 = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
  __m256i v_tmp1 = _mm256_packs_epi32(v_trunc[2], v_trunc[3]);

  dst[0] = _mm256_permute2x128_si256(v_tmp0, v_tmp1, 0x20);
  dst[1] = _mm256_permute2x128_si256(v_tmp0, v_tmp1, 0x31);
}

static void fast_inverse_tr_2x16_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0415);

  __m256i v_src_lo = _mm256_unpacklo_epi16(src[0], src[1]);
  __m256i v_src_hi = _mm256_unpackhi_epi16(src[0], src[1]);

  __m256i v_madd_lo_0 = _mm256_madd_epi16(v_src_lo, v_coeff[0]);
  __m256i v_madd_lo_1 = _mm256_madd_epi16(v_src_lo, v_coeff[1]);

  __m256i v_madd_hi_0 = _mm256_madd_epi16(v_src_hi, v_coeff[0]);
  __m256i v_madd_hi_1 = _mm256_madd_epi16(v_src_hi, v_coeff[1]);

  __m256i v_trunc_0 = truncate_avx2(v_madd_lo_0, debias, shift);
  __m256i v_trunc_1 = truncate_avx2(v_madd_lo_1, debias, shift);
  __m256i v_trunc_2 = truncate_avx2(v_madd_hi_0, debias, shift);
  __m256i v_trunc_3 = truncate_avx2(v_madd_hi_1, debias, shift);

  __m256i v_tmp0 = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
  __m256i v_tmp1 = _mm256_packs_epi32(v_trunc_2, v_trunc_3);

  v_tmp0 = _mm256_shuffle_epi8(v_tmp0, v_res_shuffle);
  v_tmp1 = _mm256_shuffle_epi8(v_tmp1, v_res_shuffle);

  __m256i v_result_0 = _mm256_permute2x128_si256(v_tmp0, v_tmp1, 0x20);
  __m256i v_result_1 = _mm256_permute2x128_si256(v_tmp0, v_tmp1, 0x31);

  _mm256_store_si256((__m256i*) & dst[0], v_result_0);
  _mm256_store_si256((__m256i*) & dst[16], v_result_1);
}

static void fast_inverse_tr_2x16_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 2;
  const int height = 16;

  int skip_width = 0;
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = fi_dct2_16x2_coeff_hor; // TODO: rename
  if (ver == DST7) {
    ver_coeff = fi_dst7_16x2_coeff_hor;
  }
  const int16_t* hor_coeff = fi_dct2_16x2_coeff_ver; // rename
  // No coeffs for DCT8 and DST7 transforms since they do not exist for this block size

  __m256i v_ver_pass_out[2];
  fast_inverse_tr_2x16_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);

  fast_inverse_tr_2x16_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static void fast_forward_tr_2x32_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 2;
  const int height = 32;

  int skip_width = 0;
  int skip_height = (ver != DCT2 && height == 32) ? 16 : (height > 32 ? height - 32 : 0);

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_2xN_coeff_hor;
  const int16_t* ver_coeff = &uvg_g_dct_32[0][0];
  // For result shuffling, can use existing shuffle vector
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)ff_dct2_2x16_ver_result_shuffle);
  // No coeffs for DCT8 and DST7 transforms since they do not exist for this block size

  ALIGNED(32) int16_t v_hor_pass_out[2*32];
  fast_forward_tr_2xN_avx2_hor(src, (__m256i *)v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);

  __m256i temp_out[4];
  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  for (int j = 0; j < 2; ++j) {
    __m256i res_0 = _mm256_setzero_si256();
    __m256i res_1 = _mm256_setzero_si256();
    __m256i res_2 = _mm256_setzero_si256();
    __m256i res_3 = _mm256_setzero_si256();
    const int16_t* coeff_start = ff_dct2_32x32_coeff_ver;
    const int32_t* temp_source = (int32_t*)(v_hor_pass_out + j * 4);
    for (int i = 0; i < 16; ++i) {

      __m256i v_src = _mm256_set1_epi32(*temp_source);
      temp_source += i & 1 ? 3 : 1;
      __m256i v_coeff_0 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;
      __m256i v_coeff_1 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;
      __m256i v_coeff_2 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;
      __m256i v_coeff_3 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;

      __m256i madd_0 = _mm256_madd_epi16(v_src, v_coeff_0);
      __m256i madd_1 = _mm256_madd_epi16(v_src, v_coeff_1);
      __m256i madd_2 = _mm256_madd_epi16(v_src, v_coeff_2);
      __m256i madd_3 = _mm256_madd_epi16(v_src, v_coeff_3);

      res_0 = _mm256_add_epi32(res_0, madd_0);
      res_1 = _mm256_add_epi32(res_1, madd_1);
      res_2 = _mm256_add_epi32(res_2, madd_2);
      res_3 = _mm256_add_epi32(res_3, madd_3);
    }
    __m256i v_trunc_0 = truncate_avx2(res_0, debias, shift_2nd);
    __m256i v_trunc_1 = truncate_avx2(res_1, debias, shift_2nd);
    __m256i v_trunc_2 = truncate_avx2(res_2, debias, shift_2nd);
    __m256i v_trunc_3 = truncate_avx2(res_3, debias, shift_2nd);

    v_trunc_0 = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
    v_trunc_1 = _mm256_packs_epi32(v_trunc_2, v_trunc_3);
    v_trunc_0 = _mm256_permute4x64_epi64(v_trunc_0, _MM_SHUFFLE(3, 1, 2, 0));
    v_trunc_1 = _mm256_permute4x64_epi64(v_trunc_1, _MM_SHUFFLE(3, 1, 2, 0));
    _mm256_store_si256(temp_out + 2 * j, v_trunc_0);
    _mm256_store_si256(temp_out + 2 * j + 1, v_trunc_1);
  }
  transpose_avx2(temp_out, (__m256i*) dst, 32, 2);
}


static void fast_inverse_tr_2x32_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const int64_t* c_ptr = (const int64_t*)coeff; // Handle as 64 bit integer to load four coeffs into vector at the same time
  const __m256i* v_src_raw = (const __m256i*)src;
  const __m256i v_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0246);

  __m256i v_src[4];
  for (int i = 0; i < 4; ++i) {
    v_src[i] = _mm256_shuffle_epi8(v_src_raw[i], v_shuffle);
  }
  for (int i = 0; i < 4; ++i) {
    v_src[i] = _mm256_permute4x64_epi64(v_src[i], _MM_SHUFFLE(3, 1, 2, 0));
  }

  __m256i v_add[32];
  for (int c = 0; c < 32; c++) {
    const __m256i v_coeff_0 = _mm256_setr_epi64x(c_ptr[0], c_ptr[1], c_ptr[0], c_ptr[1]);
    const __m256i v_coeff_1 = _mm256_setr_epi64x(c_ptr[2], c_ptr[3], c_ptr[2], c_ptr[3]);
    const __m256i v_coeff_2 = _mm256_setr_epi64x(c_ptr[4], c_ptr[5], c_ptr[4], c_ptr[5]);
    const __m256i v_coeff_3 = _mm256_setr_epi64x(c_ptr[6], c_ptr[7], c_ptr[6], c_ptr[7]);

    __m256i v_madd_0 = _mm256_madd_epi16(v_src[0], v_coeff_0);
    __m256i v_madd_1 = _mm256_madd_epi16(v_src[1], v_coeff_1);
    __m256i v_madd_2 = _mm256_madd_epi16(v_src[2], v_coeff_2);
    __m256i v_madd_3 = _mm256_madd_epi16(v_src[3], v_coeff_3);

    __m256i v_add_00 = _mm256_add_epi32(v_madd_0, v_madd_1);
    __m256i v_add_01 = _mm256_add_epi32(v_madd_2, v_madd_3);

    v_add[c] = _mm256_add_epi32(v_add_00, v_add_01);
    c_ptr += 8;
  }

  __m256i v_hadd_0[16];
  for (int d = 0, s = 0; d < 16; ++d, s += 2) {
    v_hadd_0[d] = _mm256_hadd_epi32(v_add[s + 0], v_add[s + 1]);
  }

  __m256i v_hadd_1[8];
  for (int d = 0, s = 0; d < 8; ++d, s += 2) {
    v_hadd_1[d] = _mm256_hadd_epi32(v_hadd_0[s + 0], v_hadd_0[s + 1]);
  }

  __m256i v_trunc[8];
  for (int i = 0; i < 8; ++i) {
    v_trunc[i] = truncate_avx2(v_hadd_1[i], debias, shift);
  }

  __m256i v_tmp0 = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
  __m256i v_tmp1 = _mm256_packs_epi32(v_trunc[2], v_trunc[3]);
  __m256i v_tmp2 = _mm256_packs_epi32(v_trunc[4], v_trunc[5]);
  __m256i v_tmp3 = _mm256_packs_epi32(v_trunc[6], v_trunc[7]);

  dst[0] = _mm256_permute2x128_si256(v_tmp0, v_tmp1, 0x20);
  dst[1] = _mm256_permute2x128_si256(v_tmp2, v_tmp3, 0x20);
  dst[2] = _mm256_permute2x128_si256(v_tmp0, v_tmp1, 0x31);
  dst[3] = _mm256_permute2x128_si256(v_tmp2, v_tmp3, 0x31);
}

static void fast_inverse_tr_2x32_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i* v_src_raw = src;
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0415);

  const __m256i v_src_lo0 = _mm256_unpacklo_epi16(v_src_raw[0], v_src_raw[2]);
  const __m256i v_src_lo1 = _mm256_unpacklo_epi16(v_src_raw[1], v_src_raw[3]);
  const __m256i v_src_hi0 = _mm256_unpackhi_epi16(v_src_raw[0], v_src_raw[2]);
  const __m256i v_src_hi1 = _mm256_unpackhi_epi16(v_src_raw[1], v_src_raw[3]);

  __m256i v_trunc_lo_00 = truncate_avx2(_mm256_madd_epi16(v_src_lo0, v_coeff[0]), debias, shift);
  __m256i v_trunc_lo_01 = truncate_avx2(_mm256_madd_epi16(v_src_lo0, v_coeff[1]), debias, shift);
  __m256i v_trunc_lo_10 = truncate_avx2(_mm256_madd_epi16(v_src_lo1, v_coeff[0]), debias, shift);
  __m256i v_trunc_lo_11 = truncate_avx2(_mm256_madd_epi16(v_src_lo1, v_coeff[1]), debias, shift);
  __m256i v_trunc_hi_00 = truncate_avx2(_mm256_madd_epi16(v_src_hi0, v_coeff[0]), debias, shift);
  __m256i v_trunc_hi_01 = truncate_avx2(_mm256_madd_epi16(v_src_hi0, v_coeff[1]), debias, shift);
  __m256i v_trunc_hi_10 = truncate_avx2(_mm256_madd_epi16(v_src_hi1, v_coeff[0]), debias, shift);
  __m256i v_trunc_hi_11 = truncate_avx2(_mm256_madd_epi16(v_src_hi1, v_coeff[1]), debias, shift);

  __m256i v_result[4];
  __m256i v_tmp[4];
  v_tmp[0] = _mm256_shuffle_epi8(_mm256_packs_epi32(v_trunc_lo_00, v_trunc_lo_01), v_res_shuffle);
  v_tmp[1] = _mm256_shuffle_epi8(_mm256_packs_epi32(v_trunc_lo_10, v_trunc_lo_11), v_res_shuffle);
  v_tmp[2] = _mm256_shuffle_epi8(_mm256_packs_epi32(v_trunc_hi_00, v_trunc_hi_01), v_res_shuffle);
  v_tmp[3] = _mm256_shuffle_epi8(_mm256_packs_epi32(v_trunc_hi_10, v_trunc_hi_11), v_res_shuffle);

  v_result[0] = _mm256_permute2x128_si256(v_tmp[0], v_tmp[2], 0x20);
  v_result[1] = _mm256_permute2x128_si256(v_tmp[0], v_tmp[2], 0x31);
  v_result[2] = _mm256_permute2x128_si256(v_tmp[1], v_tmp[3], 0x20);
  v_result[3] = _mm256_permute2x128_si256(v_tmp[1], v_tmp[3], 0x31);

  for (int i = 0; i < 4; ++i) {
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }
}

static void fast_inverse_tr_2x32_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 2;
  const int height = 32;

  int skip_width = 0;
  int skip_height = (ver != DCT2 && height == 32) ? 16 : (height > 32 ? height - 32 : 0);

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = &uvg_g_dct_32_t[0][0]; // rename
  const int16_t* hor_coeff = fi_dct2_32x2_coeff_ver; // TODO: rename
  // No coeffs for DCT8 and DST7 transforms since they do not exist for this block size

  __m256i v_ver_pass_out[4];
  fast_inverse_tr_2x32_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);
  
  fast_inverse_tr_2x32_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
  
}


static void fast_forward_tr_4xN_avx2_hor(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = (shift > 0) ? (1 << (shift - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i v_coeff_0 = _mm256_load_si256((const __m256i*) & coeff[0]);
  const __m256i v_coeff_1 = _mm256_load_si256((const __m256i*) & coeff[16]);
  const __m256i v_coeff_2 = _mm256_load_si256((const __m256i*) & coeff[32]);
  const __m256i v_coeff_3 = _mm256_load_si256((const __m256i*) & coeff[48]);

  const __m256i v_permute_0 = _mm256_load_si256((__m256i*)ff_dct2_b4_permute_0);
  const __m256i v_permute_1 = _mm256_load_si256((__m256i*)ff_dct2_b4_permute_1);

  const int reduced_line = line - skip_line;
  // Handle 4 lines at a time (16 samples, 4 samples per line)
  for (int j = 0; j < reduced_line; j += 4) {
    //                  line 0          line 1            line 2          line 3
    // src vector:     [s00 s01 s02 s03 s04 s05 s06 s07 | s08 s09 s10 s11 s12 s13 s14 s15]
    __m256i v_src_raw = _mm256_load_si256((const __m256i*) src);

    // Arrange data for column-wise calculation. Data and coeffs are ordered so no further shuffling
    // or permutes are needed.
    // vec 1 : [s00 s01 s04 s05 s08 s09 s12 s13 | s00 s01 s04 s05 s08 s09 s12 s13]
    // vec 2 : [s02 s03 s06 s07 s10 s11 s14 s15 | s02 s03 s06 s07 s10 s11 s14 s15]
    __m256i v_src_0 = _mm256_permutevar8x32_epi32(v_src_raw, v_permute_0);
    __m256i v_src_1 = _mm256_permutevar8x32_epi32(v_src_raw, v_permute_1);

    __m256i v_madd_0 = _mm256_madd_epi16(v_src_0, v_coeff_0);
    __m256i v_madd_1 = _mm256_madd_epi16(v_src_1, v_coeff_1);
    __m256i v_madd_2 = _mm256_madd_epi16(v_src_0, v_coeff_2);
    __m256i v_madd_3 = _mm256_madd_epi16(v_src_1, v_coeff_3);


    __m256i v_add_0 = _mm256_add_epi32(v_madd_0, v_madd_1);
    __m256i v_add_1 = _mm256_add_epi32(v_madd_2, v_madd_3);

    __m256i v_trunc_0 = truncate_avx2(v_add_0, debias, shift);
    __m256i v_trunc_1 = truncate_avx2(v_add_1, debias, shift);

    dst[0] = _mm256_packs_epi32(v_trunc_0, v_trunc_1);

    src += 16;
    dst += 1;
  }
}

static void fast_forward_tr_4x4_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 4;
  const int height = 4;
  
  int skip_width = 0;
  int skip_height = 0;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  // TODO: coeffs for DST7 and DCT8
  const int16_t* hor_coeff = fast_forward_dct2_b4_coeff;
  const int16_t* ver_coeff = fast_forward_dct2_b4_coeff;
  if (hor == DST7) {
    hor_coeff = fast_forward_dst7_b4_coeff;
  }
  else if (hor == DCT8) {
    hor_coeff = fast_forward_dct8_b4_coeff;
  }
  if (ver == DST7) {
    ver_coeff = fast_forward_dst7_b4_coeff;
  }
  else if (ver == DCT8) {
    ver_coeff = fast_forward_dct8_b4_coeff;
  }

  __m256i v_hor_pass_out;
  fast_forward_tr_4xN_avx2_hor(src, &v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0;
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i v_coeff_0 = _mm256_load_si256((const __m256i*) & ver_coeff[0]);
  const __m256i v_coeff_1 = _mm256_load_si256((const __m256i*) & ver_coeff[16]);
  const __m256i v_coeff_2 = _mm256_load_si256((const __m256i*) & ver_coeff[32]);
  const __m256i v_coeff_3 = _mm256_load_si256((const __m256i*) & ver_coeff[48]);

  const __m256i v_permute_0 = _mm256_load_si256((__m256i*)ff_dct2_b4_permute_0);
  const __m256i v_permute_1 = _mm256_load_si256((__m256i*)ff_dct2_b4_permute_1);

  __m256i v_src_0 = _mm256_permutevar8x32_epi32(v_hor_pass_out, v_permute_0);
  __m256i v_src_1 = _mm256_permutevar8x32_epi32(v_hor_pass_out, v_permute_1);

  __m256i v_madd_0 = _mm256_madd_epi16(v_src_0, v_coeff_0);
  __m256i v_madd_1 = _mm256_madd_epi16(v_src_1, v_coeff_1);
  __m256i v_madd_2 = _mm256_madd_epi16(v_src_0, v_coeff_2);
  __m256i v_madd_3 = _mm256_madd_epi16(v_src_1, v_coeff_3);

  __m256i v_add_0 = _mm256_add_epi32(v_madd_0, v_madd_1);
  __m256i v_add_1 = _mm256_add_epi32(v_madd_2, v_madd_3);

  __m256i v_trunc_0 = truncate_avx2(v_add_0, debias, shift_2nd);
  __m256i v_trunc_1 = truncate_avx2(v_add_1, debias, shift_2nd);

  __m256i v_result = _mm256_packs_epi32(v_trunc_0, v_trunc_1);

  _mm256_store_si256((__m256i*)dst, v_result);
}


static void fast_inverse_tr_4x4_avx2_hor(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i v_shuffle = _mm256_load_si256((const __m256i*)fi_tr_4x4_shuffle_hor);

  const __m256i v_src_raw = _mm256_load_si256((const __m256i*)src);
  __m256i v_src = _mm256_shuffle_epi8(v_src_raw, v_shuffle);
  v_src         = _mm256_permute4x64_epi64(v_src, _MM_SHUFFLE(3, 1, 2, 0));
  v_src         = _mm256_shuffle_epi32(v_src, _MM_SHUFFLE(3, 1, 2, 0));

  __m256i v_madd_0 = _mm256_madd_epi16(v_src, v_coeff[0]);
  __m256i v_madd_1 = _mm256_madd_epi16(v_src, v_coeff[1]);
  __m256i v_madd_2 = _mm256_madd_epi16(v_src, v_coeff[2]);
  __m256i v_madd_3 = _mm256_madd_epi16(v_src, v_coeff[3]);

  __m256i v_trunc_0 = truncate_avx2(_mm256_hadd_epi32(v_madd_0, v_madd_1), debias, shift);
  __m256i v_trunc_1 = truncate_avx2(_mm256_hadd_epi32(v_madd_2, v_madd_3), debias, shift);

  dst[0] = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
}

static void fast_inverse_tr_4x4_avx2_ver(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = (shift > 0) ? (1 << (shift - 1)) : 0;
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)fi_tr_4x4_result_shuffle_ver);

  __m256i v_src = _mm256_permute4x64_epi64(src[0], _MM_SHUFFLE(3, 1, 2, 0));
  v_src         = _mm256_shuffle_epi32(v_src, _MM_SHUFFLE(3, 1, 2, 0));

  __m256i v_madd_0 = _mm256_madd_epi16(v_src, v_coeff[0]);
  __m256i v_madd_1 = _mm256_madd_epi16(v_src, v_coeff[1]);
  __m256i v_madd_2 = _mm256_madd_epi16(v_src, v_coeff[2]);
  __m256i v_madd_3 = _mm256_madd_epi16(v_src, v_coeff[3]);

  __m256i v_trunc_0 = truncate_avx2(_mm256_hadd_epi32(v_madd_0, v_madd_1), debias, shift);
  __m256i v_trunc_1 = truncate_avx2(_mm256_hadd_epi32(v_madd_2, v_madd_3), debias, shift);

  __m256i v_result = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
  v_result         = _mm256_shuffle_epi8(v_result, v_res_shuffle);

  _mm256_store_si256((__m256i*)dst, v_result);
}

static void fast_inverse_tr_4x4_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 4;
  const int height = 4;

  int skip_width = 0;
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* hor_coeff = fi_dct2_4xN_coeff_hor;
  const int16_t* ver_coeff = fi_dct2_4xN_coeff_hor; // Can use same table for both passes
  if (hor == DST7) {
    hor_coeff = fi_dst7_4xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = fi_dct8_4xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = fi_dst7_4xN_coeff_hor;
  } else if (ver == DCT8) {
    ver_coeff = fi_dct8_4xN_coeff_hor;
  }

  __m256i v_hor_pass_out;
  fast_inverse_tr_4x4_avx2_hor(src, &v_hor_pass_out, ver_coeff, shift_1st, height, 0, skip_width);

  fast_inverse_tr_4x4_avx2_ver(&v_hor_pass_out, dst, hor_coeff, shift_2nd, width, skip_width, skip_height);
}


static void fast_forward_tr_4x8_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 4;
  const int height = 8;

  int skip_width = 0;
  int skip_height = 0;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = fast_forward_dct2_b4_coeff;
  const int16_t* ver_coeff = ff_dct2_4x8_coeff_ver;
  if (hor == DST7) {
    hor_coeff = fast_forward_dst7_b4_coeff;
  } else if (hor == DCT8) {
    hor_coeff = fast_forward_dct8_b4_coeff;
  }
  if (ver == DST7) {
    ver_coeff = ff_dst7_4x8_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = ff_dct8_4x8_coeff_ver;
  }

  __m256i v_hor_pass_out[2];
  fast_forward_tr_4xN_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0;
  const __m256i debias = _mm256_set1_epi32(add);
  const __m256i* v_coeff = (const __m256i*)ver_coeff;
  
  __m256i v_madd[2][8];
  for (int i = 0; i < 8; ++i) {
    v_madd[0][i] = _mm256_madd_epi16(v_hor_pass_out[0], v_coeff[0]);
    v_madd[1][i] = _mm256_madd_epi16(v_hor_pass_out[1], v_coeff[1]);
    v_coeff += 2;
  }

  __m256i v_add[8];
  for (int i = 0; i < 8; ++i) {
    v_add[i] = _mm256_add_epi32(v_madd[0][i], v_madd[1][i]);
  }

  __m256i v_trunc[4];
  for (int dst = 0, src = 0; dst < 4; ++dst, src += 2) {
    v_trunc[dst] = truncate_avx2(_mm256_hadd_epi32(v_add[src + 0], v_add[src + 1]), debias, shift_2nd);
  }

  __m256i v_result[2];
  v_result[0] = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
  v_result[1] = _mm256_packs_epi32(v_trunc[2], v_trunc[3]);

  // Order results
  v_result[0] = _mm256_permute4x64_epi64(v_result[0], _MM_SHUFFLE(3, 1, 2, 0));
  v_result[1] = _mm256_permute4x64_epi64(v_result[1], _MM_SHUFFLE(3, 1, 2, 0));

  v_result[0] = _mm256_shuffle_epi32(v_result[0], _MM_SHUFFLE(3, 1, 2, 0));
  v_result[1] = _mm256_shuffle_epi32(v_result[1], _MM_SHUFFLE(3, 1, 2, 0));

  _mm256_store_si256((__m256i*)&dst[0],  v_result[0]);
  _mm256_store_si256((__m256i*)&dst[16], v_result[1]);
}


static void fast_inverse_tr_4x8_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i v_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0415);
  const __m256i v_permute = _mm256_load_si256((const __m256i*)permute_32b_0415);

  const __m256i* v_src_raw = (const __m256i*)src;

  __m256i v_src_0 = _mm256_shuffle_epi8(v_src_raw[0], v_shuffle);
  __m256i v_src_1 = _mm256_shuffle_epi8(v_src_raw[1], v_shuffle);
  v_src_0 = _mm256_permutevar8x32_epi32(v_src_0, v_permute);
  v_src_1 = _mm256_permutevar8x32_epi32(v_src_1, v_permute);

  __m256i v_madd_00 = _mm256_madd_epi16(v_src_0, v_coeff[0]);
  __m256i v_madd_10 = _mm256_madd_epi16(v_src_1, v_coeff[1]);

  __m256i v_madd_01 = _mm256_madd_epi16(v_src_0, v_coeff[2]);
  __m256i v_madd_11 = _mm256_madd_epi16(v_src_1, v_coeff[3]);

  __m256i v_madd_02 = _mm256_madd_epi16(v_src_0, v_coeff[4]);
  __m256i v_madd_12 = _mm256_madd_epi16(v_src_1, v_coeff[5]);

  __m256i v_madd_03 = _mm256_madd_epi16(v_src_0, v_coeff[6]);
  __m256i v_madd_13 = _mm256_madd_epi16(v_src_1, v_coeff[7]);

  __m256i v_madd_04 = _mm256_madd_epi16(v_src_0, v_coeff[8]);
  __m256i v_madd_14 = _mm256_madd_epi16(v_src_1, v_coeff[9]);

  __m256i v_madd_05 = _mm256_madd_epi16(v_src_0, v_coeff[10]);
  __m256i v_madd_15 = _mm256_madd_epi16(v_src_1, v_coeff[11]);

  __m256i v_madd_06 = _mm256_madd_epi16(v_src_0, v_coeff[12]);
  __m256i v_madd_16 = _mm256_madd_epi16(v_src_1, v_coeff[13]);

  __m256i v_madd_07 = _mm256_madd_epi16(v_src_0, v_coeff[14]);
  __m256i v_madd_17 = _mm256_madd_epi16(v_src_1, v_coeff[15]);

  __m256i v_add_0 = _mm256_add_epi32(v_madd_00, v_madd_10);
  __m256i v_add_1 = _mm256_add_epi32(v_madd_01, v_madd_11);
  __m256i v_add_2 = _mm256_add_epi32(v_madd_02, v_madd_12);
  __m256i v_add_3 = _mm256_add_epi32(v_madd_03, v_madd_13);
  __m256i v_add_4 = _mm256_add_epi32(v_madd_04, v_madd_14);
  __m256i v_add_5 = _mm256_add_epi32(v_madd_05, v_madd_15);
  __m256i v_add_6 = _mm256_add_epi32(v_madd_06, v_madd_16);
  __m256i v_add_7 = _mm256_add_epi32(v_madd_07, v_madd_17);

  __m256i v_hadd_0 = _mm256_hadd_epi32(v_add_0, v_add_1);
  __m256i v_hadd_1 = _mm256_hadd_epi32(v_add_2, v_add_3);
  __m256i v_hadd_2 = _mm256_hadd_epi32(v_add_4, v_add_5);
  __m256i v_hadd_3 = _mm256_hadd_epi32(v_add_6, v_add_7);

  __m256i v_trunc_0 = truncate_avx2(v_hadd_0, debias, shift);
  __m256i v_trunc_1 = truncate_avx2(v_hadd_1, debias, shift);
  __m256i v_trunc_2 = truncate_avx2(v_hadd_2, debias, shift);
  __m256i v_trunc_3 = truncate_avx2(v_hadd_3, debias, shift);

  dst[0] = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
  dst[1] = _mm256_packs_epi32(v_trunc_2, v_trunc_3);
}

static void fast_inverse_tr_4x8_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0415);

  __m256i v_src_0 = _mm256_permute2x128_si256(src[0], src[1], 0x20);
  __m256i v_src_1 = _mm256_permute2x128_si256(src[0], src[1], 0x31);

  __m256i v_madd_00 = _mm256_madd_epi16(v_src_0, v_coeff[0]);
  __m256i v_madd_10 = _mm256_madd_epi16(v_src_1, v_coeff[1]);

  __m256i v_madd_01 = _mm256_madd_epi16(v_src_0, v_coeff[2]);
  __m256i v_madd_11 = _mm256_madd_epi16(v_src_1, v_coeff[3]);

  __m256i v_madd_02 = _mm256_madd_epi16(v_src_0, v_coeff[4]);
  __m256i v_madd_12 = _mm256_madd_epi16(v_src_1, v_coeff[5]);

  __m256i v_madd_03 = _mm256_madd_epi16(v_src_0, v_coeff[6]);
  __m256i v_madd_13 = _mm256_madd_epi16(v_src_1, v_coeff[7]);

  __m256i v_trunc_0 = truncate_avx2(_mm256_add_epi32(v_madd_00, v_madd_10), debias, shift);
  __m256i v_trunc_1 = truncate_avx2(_mm256_add_epi32(v_madd_01, v_madd_11), debias, shift);
  __m256i v_trunc_2 = truncate_avx2(_mm256_add_epi32(v_madd_02, v_madd_12), debias, shift);
  __m256i v_trunc_3 = truncate_avx2(_mm256_add_epi32(v_madd_03, v_madd_13), debias, shift);

  __m256i v_tmp0 = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
  __m256i v_tmp1 = _mm256_packs_epi32(v_trunc_2, v_trunc_3);

  __m256i v_result_0 = _mm256_permute2x128_si256(v_tmp0, v_tmp1, 0x20);
  __m256i v_result_1 = _mm256_permute2x128_si256(v_tmp0, v_tmp1, 0x31);

  v_result_0 = _mm256_shuffle_epi8(v_result_0, v_res_shuffle);
  v_result_1 = _mm256_shuffle_epi8(v_result_1, v_res_shuffle);

  v_result_0 = _mm256_permute4x64_epi64(v_result_0, _MM_SHUFFLE(3, 1, 2, 0));
  v_result_1 = _mm256_permute4x64_epi64(v_result_1, _MM_SHUFFLE(3, 1, 2, 0));

  v_result_0 = _mm256_shuffle_epi32(v_result_0, _MM_SHUFFLE(3, 1, 2, 0));
  v_result_1 = _mm256_shuffle_epi32(v_result_1, _MM_SHUFFLE(3, 1, 2, 0));

  _mm256_store_si256((__m256i*) & dst[0], v_result_0);
  _mm256_store_si256((__m256i*) & dst[16], v_result_1);
}

static void fast_inverse_tr_4x8_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 4;
  const int height = 8;

  int skip_width = 0;
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = fi_dct2_8x4_coeff_hor; // TODO: rename coeff tables
  const int16_t* hor_coeff = fi_dct2_8x4_coeff_ver;
  if (hor == DST7) {
    hor_coeff = fi_dst7_8x4_coeff_ver;
  } else if (hor == DCT8) {
    hor_coeff = fi_dct8_8x4_coeff_ver;
  }
  if (ver == DST7) {
    ver_coeff = fi_dst7_8x4_coeff_hor;
  } else if (ver == DCT8) {
    ver_coeff = fi_dct8_8x4_coeff_hor;
  }

  __m256i v_ver_pass_out[2];
  fast_inverse_tr_4x8_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);
  
  fast_inverse_tr_4x8_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static void fast_forward_tr_4x16_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 4;
  const int height = 16;

  int skip_width = 0;
  int skip_height = 0;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = fast_forward_dct2_b4_coeff;
  const int16_t* ver_coeff = &uvg_g_dct_16[0][0];
  if (hor == DST7) {
    hor_coeff = fast_forward_dst7_b4_coeff;
  } else if (hor == DCT8) {
    hor_coeff = fast_forward_dct8_b4_coeff;
  }
  if (ver == DST7) {
    ver_coeff = &uvg_g_dst7_16[0][0];
  } else if (ver == DCT8) {
    ver_coeff = &uvg_g_dct8_16[0][0];
  }

  __m256i v_hor_pass_out[4];
  fast_forward_tr_4xN_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0;
  const __m256i debias = _mm256_set1_epi32(add);
  const int64_t* coeff_ptr = (const int64_t*)ver_coeff; // Read four coeffs at once by casting into 64 bit integer

  __m256i v_madd[4][16];
  for (int i = 0; i < 16; ++i) {
    const __m256i v_coeff_0 = _mm256_set1_epi64x(coeff_ptr[0]);
    const __m256i v_coeff_1 = _mm256_set1_epi64x(coeff_ptr[1]);
    const __m256i v_coeff_2 = _mm256_set1_epi64x(coeff_ptr[2]);
    const __m256i v_coeff_3 = _mm256_set1_epi64x(coeff_ptr[3]);
    v_madd[0][i] = _mm256_madd_epi16(v_hor_pass_out[0], v_coeff_0);
    v_madd[1][i] = _mm256_madd_epi16(v_hor_pass_out[1], v_coeff_1);
    v_madd[2][i] = _mm256_madd_epi16(v_hor_pass_out[2], v_coeff_2);
    v_madd[3][i] = _mm256_madd_epi16(v_hor_pass_out[3], v_coeff_3);
    coeff_ptr += 4;
  }

  __m256i v_add[16];
  for (int i = 0; i < 16; ++i) {
    __m256i v_tmp0 = _mm256_add_epi32(v_madd[0][i], v_madd[1][i]);
    __m256i v_tmp1 = _mm256_add_epi32(v_madd[2][i], v_madd[3][i]);

    v_add[i] = _mm256_add_epi32(v_tmp0, v_tmp1);
  }

  __m256i v_trunc[8];
  for (int dst = 0, src = 0; dst < 8; ++dst, src += 2) {
    v_trunc[dst] = truncate_avx2(_mm256_hadd_epi32(v_add[src + 0], v_add[src + 1]), debias, shift_2nd);
  }

  __m256i v_result[4];
  v_result[0] = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
  v_result[1] = _mm256_packs_epi32(v_trunc[2], v_trunc[3]);
  v_result[2] = _mm256_packs_epi32(v_trunc[4], v_trunc[5]);
  v_result[3] = _mm256_packs_epi32(v_trunc[6], v_trunc[7]);

  for (int i = 0; i < 4; ++i) {
    v_result[i] = _mm256_permute4x64_epi64(v_result[i], _MM_SHUFFLE(3, 1, 2, 0));
  }

  for (int i = 0; i < 4; ++i) {
    v_result[i] = _mm256_shuffle_epi32(v_result[i], _MM_SHUFFLE(3, 1, 2, 0));
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }
}


static void fast_inverse_tr_4x16_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i* v_src_raw = (const __m256i*)src;
  const __m256i v_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0415);

  __m256i v_src_0 = _mm256_shuffle_epi8(v_src_raw[0], v_shuffle);
  __m256i v_src_1 = _mm256_shuffle_epi8(v_src_raw[1], v_shuffle);
  __m256i v_src_2 = _mm256_shuffle_epi8(v_src_raw[2], v_shuffle);
  __m256i v_src_3 = _mm256_shuffle_epi8(v_src_raw[3], v_shuffle);

  v_src_0 = _mm256_permute4x64_epi64(v_src_0, _MM_SHUFFLE(3, 1, 2, 0));
  v_src_1 = _mm256_permute4x64_epi64(v_src_1, _MM_SHUFFLE(3, 1, 2, 0));
  v_src_2 = _mm256_permute4x64_epi64(v_src_2, _MM_SHUFFLE(3, 1, 2, 0));
  v_src_3 = _mm256_permute4x64_epi64(v_src_3, _MM_SHUFFLE(3, 1, 2, 0));

  v_src_0 = _mm256_shuffle_epi32(v_src_0, _MM_SHUFFLE(3, 1, 2, 0));
  v_src_1 = _mm256_shuffle_epi32(v_src_1, _MM_SHUFFLE(3, 1, 2, 0));
  v_src_2 = _mm256_shuffle_epi32(v_src_2, _MM_SHUFFLE(3, 1, 2, 0));
  v_src_3 = _mm256_shuffle_epi32(v_src_3, _MM_SHUFFLE(3, 1, 2, 0));

  __m256i v_madd_0[16];
  __m256i v_madd_1[16];
  __m256i v_madd_2[16];
  __m256i v_madd_3[16];
  for (int c = 0; c < 16; c++) {
    v_madd_0[c] = _mm256_madd_epi16(v_src_0, v_coeff[0]);
    v_madd_1[c] = _mm256_madd_epi16(v_src_1, v_coeff[1]);
    v_madd_2[c] = _mm256_madd_epi16(v_src_2, v_coeff[2]);
    v_madd_3[c] = _mm256_madd_epi16(v_src_3, v_coeff[3]);
    v_coeff += 4;
  }

  __m256i v_add[16];
  for (int i = 0; i < 16; ++i) {
    __m256i v_add_0 = _mm256_add_epi32(v_madd_0[i], v_madd_1[i]);
    __m256i v_add_1 = _mm256_add_epi32(v_madd_2[i], v_madd_3[i]);

    v_add[i] = _mm256_add_epi32(v_add_0, v_add_1);
  }

  __m256i v_hadd[8];
  for (int d = 0, s = 0; d < 8; ++d, s += 2) {
    v_hadd[d] = _mm256_hadd_epi32(v_add[s + 0], v_add[s + 1]);
  }

  __m256i v_trunc[8];
  for (int i = 0; i < 8; ++i) {
    v_trunc[i] = truncate_avx2(v_hadd[i], debias, shift);
  }

  dst[0] = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
  dst[1] = _mm256_packs_epi32(v_trunc[2], v_trunc[3]);
  dst[2] = _mm256_packs_epi32(v_trunc[4], v_trunc[5]);
  dst[3] = _mm256_packs_epi32(v_trunc[6], v_trunc[7]);
}

static void fast_inverse_tr_4x16_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0415);

  __m256i v_src_0 = _mm256_permute2x128_si256(src[0], src[1], 0x20);
  __m256i v_src_1 = _mm256_permute2x128_si256(src[0], src[1], 0x31);
  __m256i v_src_2 = _mm256_permute2x128_si256(src[2], src[3], 0x20);
  __m256i v_src_3 = _mm256_permute2x128_si256(src[2], src[3], 0x31);

  __m256i v_madd_0[4];
  __m256i v_madd_1[4];
  __m256i v_madd_2[4];
  __m256i v_madd_3[4];
  for (int c = 0; c < 4; ++c) {
    v_madd_0[c] = _mm256_madd_epi16(v_src_0, v_coeff[0]);
    v_madd_1[c] = _mm256_madd_epi16(v_src_1, v_coeff[1]);
    v_madd_2[c] = _mm256_madd_epi16(v_src_2, v_coeff[0]);
    v_madd_3[c] = _mm256_madd_epi16(v_src_3, v_coeff[1]);
    v_coeff += 2;
  }

  __m256i v_trunc_0[4];
  __m256i v_trunc_1[4];
  for (int i = 0; i < 4; ++i) {
    v_trunc_0[i] = truncate_avx2(_mm256_add_epi32(v_madd_0[i], v_madd_1[i]), debias, shift);
    v_trunc_1[i] = truncate_avx2(_mm256_add_epi32(v_madd_2[i], v_madd_3[i]), debias, shift);
  }

  __m256i v_result[4];
  __m256i v_tmp0 = _mm256_packs_epi32(v_trunc_0[0], v_trunc_0[1]);
  __m256i v_tmp1 = _mm256_packs_epi32(v_trunc_0[2], v_trunc_0[3]);
  __m256i v_tmp2 = _mm256_packs_epi32(v_trunc_1[0], v_trunc_1[1]);
  __m256i v_tmp3 = _mm256_packs_epi32(v_trunc_1[2], v_trunc_1[3]);

  v_tmp0 = _mm256_shuffle_epi8(v_tmp0, v_res_shuffle);
  v_tmp1 = _mm256_shuffle_epi8(v_tmp1, v_res_shuffle);
  v_tmp2 = _mm256_shuffle_epi8(v_tmp2, v_res_shuffle);
  v_tmp3 = _mm256_shuffle_epi8(v_tmp3, v_res_shuffle);

  __m256i v_tmp32_0 = _mm256_unpacklo_epi32(v_tmp0, v_tmp1);
  __m256i v_tmp32_1 = _mm256_unpackhi_epi32(v_tmp0, v_tmp1);
  __m256i v_tmp32_2 = _mm256_unpacklo_epi32(v_tmp2, v_tmp3);
  __m256i v_tmp32_3 = _mm256_unpackhi_epi32(v_tmp2, v_tmp3);

  v_result[0] = _mm256_permute2x128_si256(v_tmp32_0, v_tmp32_1, 0x20);
  v_result[1] = _mm256_permute2x128_si256(v_tmp32_0, v_tmp32_1, 0x31);
  v_result[2] = _mm256_permute2x128_si256(v_tmp32_2, v_tmp32_3, 0x20);
  v_result[3] = _mm256_permute2x128_si256(v_tmp32_2, v_tmp32_3, 0x31);

  for (int i = 0; i < 4; ++i) {
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }
}

static void fast_inverse_tr_4x16_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 4;
  const int height = 16;

  int skip_width = 0;
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = fi_dct2_16x4_coeff_hor; // TODO: rename coeff tables
  const int16_t* hor_coeff = fi_dct2_16x4_coeff_ver;
  if (hor == DST7) {
    hor_coeff = fi_dst7_16x4_coeff_ver;
  } else if (hor == DCT8) {
    hor_coeff = fi_dct8_16x4_coeff_ver;
  }
  if (ver == DST7) {
    ver_coeff = fi_dst7_16x4_coeff_hor;
  } else if (ver == DCT8) {
    ver_coeff = fi_dct8_16x4_coeff_hor;
  }

  __m256i v_ver_pass_out[4];
  fast_inverse_tr_4x16_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);

  fast_inverse_tr_4x16_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static void fast_forward_tr_4x32_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 4;
  const int height = 32;

  int skip_width = 0;
  int skip_height = (ver != DCT2 && height == 32) ? 16 : (height > 32 ? height - 32 : 0);

  const int  reduced_line = width - skip_width;
  const int  cutoff = height - skip_height;
  int16_t* p_dst = dst;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = fast_forward_dct2_b4_coeff;
  const int16_t* ver_coeff = ff_dct2_32xN_coeff_hor;
  if (hor == DST7) {
    hor_coeff = fast_forward_dst7_b4_coeff;
  } else if (hor == DCT8) {
    hor_coeff = fast_forward_dct8_b4_coeff;
  }
  if (ver == DST7) {
    ver_coeff = ff_dst7_4x32_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = ff_dct8_4x32_coeff_ver;
  }

  ALIGNED(32) int16_t v_hor_pass_out[4*32];
  fast_forward_tr_4xN_avx2_hor(src, (__m256i*)v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);


  __m256i temp_out[8];
  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  for (int j = 0; j < 4; ++j) {
    __m256i res_0 = _mm256_setzero_si256();
    __m256i res_1 = _mm256_setzero_si256();
    __m256i res_2 = _mm256_setzero_si256();
    __m256i res_3 = _mm256_setzero_si256();
    const int16_t* coeff_start = ver_coeff;
    const int32_t* temp_source = (int32_t*)(v_hor_pass_out + j * 4);
    for (int i = 0; i < 16; ++i) {

      __m256i v_src = _mm256_set1_epi32(*temp_source);
      temp_source += i & 1 ? 7 : 1;
      __m256i v_coeff_0 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;
      __m256i v_coeff_1 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;
      __m256i v_coeff_2 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;
      __m256i v_coeff_3 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;

      __m256i madd_0 = _mm256_madd_epi16(v_src, v_coeff_0);
      __m256i madd_1 = _mm256_madd_epi16(v_src, v_coeff_1);
      __m256i madd_2 = _mm256_madd_epi16(v_src, v_coeff_2);
      __m256i madd_3 = _mm256_madd_epi16(v_src, v_coeff_3);

      res_0 = _mm256_add_epi32(res_0, madd_0);
      res_1 = _mm256_add_epi32(res_1, madd_1);
      res_2 = _mm256_add_epi32(res_2, madd_2);
      res_3 = _mm256_add_epi32(res_3, madd_3);
    }
    __m256i v_trunc_0 = truncate_avx2(res_0, debias, shift_2nd);
    __m256i v_trunc_1 = truncate_avx2(res_1, debias, shift_2nd);
    __m256i v_trunc_2 = truncate_avx2(res_2, debias, shift_2nd);
    __m256i v_trunc_3 = truncate_avx2(res_3, debias, shift_2nd);

    v_trunc_0 = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
    v_trunc_1 = _mm256_packs_epi32(v_trunc_2, v_trunc_3);
    v_trunc_0 = _mm256_permute4x64_epi64(v_trunc_0, _MM_SHUFFLE(3, 1, 2, 0));
    v_trunc_1 = _mm256_permute4x64_epi64(v_trunc_1, _MM_SHUFFLE(3, 1, 2, 0));
    _mm256_store_si256(temp_out + 2 * j, v_trunc_0);
    _mm256_store_si256(temp_out + 2 * j + 1, v_trunc_1);
  }
  transpose_avx2(temp_out, (__m256i*) dst, 32, 4);


}


static void fast_inverse_tr_4x32_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const int64_t* c_ptr = (const int64_t*)coeff; // Handle as 64 bit integer to load four coeffs into vector at the same time
  const __m256i* v_src_raw = (const __m256i*)src;
  const __m256i v_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0415);

  __m256i v_src[8];
  for (int i = 0; i < 8; ++i) {
    v_src[i] = _mm256_shuffle_epi8(v_src_raw[i], v_shuffle);
  }
  for (int i = 0; i < 8; ++i) {
    v_src[i] = _mm256_permute4x64_epi64(v_src[i], _MM_SHUFFLE(3, 1, 2, 0));
    v_src[i] = _mm256_shuffle_epi32(v_src[i], _MM_SHUFFLE(3, 1, 2, 0));
  }

  __m256i v_add[32];
  for (int c = 0; c < 32; c++) {
    __m256i v_madd[8];
    for (int i = 0; i < 8; ++i) {
      const __m256i v_coeff = _mm256_set1_epi64x(*c_ptr);
      v_madd[i] = _mm256_madd_epi16(v_src[i], v_coeff);
      c_ptr++;
    }

    __m256i v_add_0[4];
    for (int d = 0, s = 0; d < 4; ++d, s += 2) {
      v_add_0[d] = _mm256_add_epi32(v_madd[s + 0], v_madd[s + 1]);
    }

    __m256i v_add_10 = _mm256_add_epi32(v_add_0[0], v_add_0[1]);
    __m256i v_add_11 = _mm256_add_epi32(v_add_0[2], v_add_0[3]);

    v_add[c] = _mm256_add_epi32(v_add_10, v_add_11);
  }

  __m256i v_hadd[16];
  for (int d = 0, s = 0; d < 16; ++d, s += 2) {
    v_hadd[d] = _mm256_hadd_epi32(v_add[s + 0], v_add[s + 1]);
  }

  __m256i v_trunc[16];
  for (int i = 0; i < 16; ++i) {
    v_trunc[i] = truncate_avx2(v_hadd[i], debias, shift);
  }

  for (int d = 0, s = 0; d < 8; ++d, s += 2) {
    dst[d] = _mm256_packs_epi32(v_trunc[s + 0], v_trunc[s + 1]);
  }
  // TODO: cutoff for dct8 and dst7
}

static void fast_inverse_tr_4x32_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i* v_src_raw = src;

  __m256i v_src[8];
  __m256i v_tmp[8];
  v_src[0] = _mm256_permute2x128_si256(v_src_raw[0], v_src_raw[1], 0x20);
  v_src[1] = _mm256_permute2x128_si256(v_src_raw[0], v_src_raw[1], 0x31);
  v_src[2] = _mm256_permute2x128_si256(v_src_raw[2], v_src_raw[3], 0x20);
  v_src[3] = _mm256_permute2x128_si256(v_src_raw[2], v_src_raw[3], 0x31);
  v_src[4] = _mm256_permute2x128_si256(v_src_raw[4], v_src_raw[5], 0x20);
  v_src[5] = _mm256_permute2x128_si256(v_src_raw[4], v_src_raw[5], 0x31);
  v_src[6] = _mm256_permute2x128_si256(v_src_raw[6], v_src_raw[7], 0x20);
  v_src[7] = _mm256_permute2x128_si256(v_src_raw[6], v_src_raw[7], 0x31);

  for (int d = 0, c = 0; c < 4; ++c, d += 2) {
    __m256i v_madd_00 = _mm256_madd_epi16(v_src[0], v_coeff[0]);
    __m256i v_madd_01 = _mm256_madd_epi16(v_src[1], v_coeff[1]);
    __m256i v_madd_10 = _mm256_madd_epi16(v_src[2], v_coeff[0]);
    __m256i v_madd_11 = _mm256_madd_epi16(v_src[3], v_coeff[1]);
    __m256i v_madd_20 = _mm256_madd_epi16(v_src[4], v_coeff[0]);
    __m256i v_madd_21 = _mm256_madd_epi16(v_src[5], v_coeff[1]);
    __m256i v_madd_30 = _mm256_madd_epi16(v_src[6], v_coeff[0]);
    __m256i v_madd_31 = _mm256_madd_epi16(v_src[7], v_coeff[1]);
    v_coeff += 2;

    __m256i v_trunc_0 = truncate_avx2(_mm256_add_epi32(v_madd_00, v_madd_01), debias, shift);
    __m256i v_trunc_1 = truncate_avx2(_mm256_add_epi32(v_madd_10, v_madd_11), debias, shift);
    __m256i v_trunc_2 = truncate_avx2(_mm256_add_epi32(v_madd_20, v_madd_21), debias, shift);
    __m256i v_trunc_3 = truncate_avx2(_mm256_add_epi32(v_madd_30, v_madd_31), debias, shift);

    v_tmp[d + 0] = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
    v_tmp[d + 1] = _mm256_packs_epi32(v_trunc_2, v_trunc_3);

    v_tmp[d + 0] = _mm256_permute4x64_epi64(v_tmp[d + 0], _MM_SHUFFLE(3, 1, 2, 0));
    v_tmp[d + 1] = _mm256_permute4x64_epi64(v_tmp[d + 1], _MM_SHUFFLE(3, 1, 2, 0));
  }

  __m256i v_result[8];
  transpose_avx2(v_tmp, v_result, 32, 4);

  for (int i = 0; i < 8; ++i) {
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }
}

static void fast_inverse_tr_4x32_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 4;
  const int height = 32;

  int skip_width = 0;
  int skip_height = (ver != DCT2 && height == 32) ? 16 : (height > 32 ? height - 32 : 0);

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = &uvg_g_dct_32_t[0][0];
  const int16_t* hor_coeff = fi_dct2_32x4_coeff_ver; // TODO: rename
  if (hor == DST7) {
    hor_coeff = fi_dst7_32x4_coeff_ver; // TODO: rename
  } else if (hor == DCT8) {
    hor_coeff = fi_dct8_32x4_coeff_ver; // TODO: rename
  }
  if (ver == DST7) {
    ver_coeff = &uvg_g_dst7_32_t[0][0];
  } else if (ver == DCT8) {
    ver_coeff = &uvg_g_dct8_32[0][0];
  }

  __m256i v_ver_pass_out[8];
  fast_inverse_tr_4x32_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);

  fast_inverse_tr_4x32_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static void fast_forward_tr_8xN_avx2_hor(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = (shift > 0) ? (1 << (shift - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;

  const int reduced_line = line - skip_line;
  // Handle 2 lines at a time (16 samples, 8 samples per line)
  for (int j = 0; j < reduced_line; j += 2) {
    //                    line 1                    line 2
    // src vector:       [s0 s1 s2 s3 s4 s5 s6 s7 | s0 s1 s2 s3 s4 s5 s6 s7]
    __m256i    v_src = _mm256_load_si256((const __m256i*)src);

    // Rearrange source in a way samples can be added together column-wise using add
    // after first round of madd operations.
    // Need 4 source vectors arranged as follows. High 128 lanes are the same as low:
    // vec_01 = [s0 s1 s0 s1 s0 s1 s0 s1 |...]
    // vec_02 = [s2 s3 s2 s3 s2 s3 s2 s3 |...]
    // vec_03 = [s4 s5 s4 s5 s4 s5 s4 s5 |...]
    // vec_04 = [s6 s7 s6 s7 s6 s7 s6 s7 |...]

    __m256i  v_src_0 = _mm256_shuffle_epi32(v_src, _MM_SHUFFLE(0, 0, 0, 0));
    __m256i  v_src_1 = _mm256_shuffle_epi32(v_src, _MM_SHUFFLE(1, 1, 1, 1));
    __m256i  v_src_2 = _mm256_shuffle_epi32(v_src, _MM_SHUFFLE(2, 2, 2, 2));
    __m256i  v_src_3 = _mm256_shuffle_epi32(v_src, _MM_SHUFFLE(3, 3, 3, 3));

    // Lane 1
    __m256i v_madd_0 = _mm256_madd_epi16(v_src_0, v_coeff[0]);
    __m256i v_madd_1 = _mm256_madd_epi16(v_src_1, v_coeff[1]);
    __m256i v_madd_2 = _mm256_madd_epi16(v_src_2, v_coeff[2]);
    __m256i v_madd_3 = _mm256_madd_epi16(v_src_3, v_coeff[3]);

    __m256i v_add_00 = _mm256_add_epi32(v_madd_0, v_madd_1);
    __m256i v_add_01 = _mm256_add_epi32(v_madd_2, v_madd_3);

    __m256i v_add_10 = _mm256_add_epi32(v_add_00, v_add_01);

    // Lane 2
    __m256i v_madd_4 = _mm256_madd_epi16(v_src_0, v_coeff[4]);
    __m256i v_madd_5 = _mm256_madd_epi16(v_src_1, v_coeff[5]);
    __m256i v_madd_6 = _mm256_madd_epi16(v_src_2, v_coeff[6]);
    __m256i v_madd_7 = _mm256_madd_epi16(v_src_3, v_coeff[7]);

    __m256i v_add_02 = _mm256_add_epi32(v_madd_4, v_madd_5);
    __m256i v_add_03 = _mm256_add_epi32(v_madd_6, v_madd_7);

    __m256i v_add_11 = _mm256_add_epi32(v_add_02, v_add_03);

    // Trunc results from both lanes
    __m256i v_trunc_0 = truncate_avx2(v_add_10, debias, shift);
    __m256i v_trunc_1 = truncate_avx2(v_add_11, debias, shift);

    dst[0] = _mm256_packs_epi32(v_trunc_0, v_trunc_1);

    src += 16;
    dst += 1;
  }
}

static void fast_forward_tr_8x2_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 8;
  const int height = 2;
  
  int skip_width  = 0;
  int skip_height = 0;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_8xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_8x2_coeff_ver;
  // Only DCT2 is defined for 8x2 block
  if (hor == DST7) {
    hor_coeff = ff_dst7_8xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = ff_dct8_8xN_coeff_hor;
  }

  __m256i v_hor_pass_out;
  fast_forward_tr_8xN_avx2_hor(src, &v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  // TODO: coeffs for DST7 and DCT8 transforms
  const __m256i* v_coeff = (const __m256i*)ver_coeff;
  const __m256i v_shuffle = _mm256_load_si256((const __m256i*)ff_dct2_8x2_ver_pass_shuffle);

  // 8x2, only 16 samples, handle all at once
  __m256i v_src_per = _mm256_permute4x64_epi64(v_hor_pass_out, _MM_SHUFFLE(3, 1, 2, 0));
  // Weave lo and hi halfs of each 128 bit lane
  __m256i     v_src = _mm256_shuffle_epi8(v_src_per, v_shuffle);
  //            v_src = _mm256_unpackhi_epi16(v_src_raw, v_src_swp);

  __m256i v_madd_0 = _mm256_madd_epi16(v_src, v_coeff[0]);
  __m256i v_madd_1 = _mm256_madd_epi16(v_src, v_coeff[1]);

  __m256i v_trunc_0 = truncate_avx2(v_madd_0, debias, shift_2nd);
  __m256i v_trunc_1 = truncate_avx2(v_madd_1, debias, shift_2nd);

  __m256i v_result = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
          v_result = _mm256_permute4x64_epi64(v_result, _MM_SHUFFLE(3, 1, 2, 0)); // TODO: this permute can probably be optimized away

  _mm256_store_si256((__m256i*)dst, v_result);
}


static void fast_inverse_tr_8x2_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i v_shuffle = _mm256_load_si256((const __m256i*)fi_tr_2x8_shuffle_hor);
  const __m256i v_coeff_0 = _mm256_load_si256((const __m256i*) & coeff[0]);
  const __m256i v_coeff_1 = _mm256_load_si256((const __m256i*) & coeff[16]);

  // Got data for one vector
  const __m256i v_src_raw = _mm256_load_si256((const __m256i*)src);

  __m256i v_src = _mm256_permute4x64_epi64(v_src_raw, _MM_SHUFFLE(3, 1, 2, 0));
  v_src = _mm256_shuffle_epi8(v_src, v_shuffle);

  __m256i    v_even = _mm256_madd_epi16(v_src, v_coeff_0);
  // odd vector : [a00-a01 a02-a03 a04-a05 a06-a07|a08-a09 a10-a11 a12-a13 a14-a15]
  __m256i     v_odd = _mm256_madd_epi16(v_src, v_coeff_1);

  __m256i v_trunc_0 = truncate_avx2(v_even, debias, shift);
  __m256i v_trunc_1 = truncate_avx2(v_odd, debias, shift);

  dst[0] = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
}

static void fast_inverse_tr_8x2_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = (shift > 0) ? (1 << (shift - 1)) : 0;
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i v_shuffle1 = _mm256_load_si256((const __m256i*)fi_tr_2x8_result_shuffle1_ver);
  const __m256i v_shuffle2 = _mm256_load_si256((const __m256i*)fi_tr_2x8_result_shuffle2_ver);

  // Duplicate sources to enable vertical addition
  __m256i v_src_0 = _mm256_permute4x64_epi64(*src, _MM_SHUFFLE(1, 1, 0, 0));
  __m256i v_src_1 = _mm256_permute4x64_epi64(*src, _MM_SHUFFLE(3, 3, 2, 2));

  __m256i v_madd_00 = _mm256_madd_epi16(v_src_0, v_coeff[0]);
  __m256i v_madd_01 = _mm256_madd_epi16(v_src_1, v_coeff[1]);
  
  __m256i v_madd_10 = _mm256_madd_epi16(v_src_0, v_coeff[2]);
  __m256i v_madd_11 = _mm256_madd_epi16(v_src_1, v_coeff[3]);

  __m256i v_madd_20 = _mm256_madd_epi16(v_src_0, v_coeff[4]);
  __m256i v_madd_21 = _mm256_madd_epi16(v_src_1, v_coeff[5]);
  
  __m256i v_madd_30 = _mm256_madd_epi16(v_src_0, v_coeff[6]);
  __m256i v_madd_31 = _mm256_madd_epi16(v_src_1, v_coeff[7]);

  __m256i v_add_0 = _mm256_add_epi32(v_madd_00, v_madd_01);
  __m256i v_add_1 = _mm256_add_epi32(v_madd_10, v_madd_11);
  __m256i v_add_2 = _mm256_add_epi32(v_madd_20, v_madd_21);
  __m256i v_add_3 = _mm256_add_epi32(v_madd_30, v_madd_31);

  __m256i v_trunc_0 = truncate_avx2(_mm256_hadd_epi32(v_add_0, v_add_1), debias, shift);
  __m256i v_trunc_1 = truncate_avx2(_mm256_hadd_epi32(v_add_2, v_add_3), debias, shift);

  __m256i v_result = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
  //v_result = _mm256_shuffle_epi8(v_result, v_shuffle1);
  //v_result = _mm256_permute4x64_epi64(v_result, _MM_SHUFFLE(3, 1, 2, 0));
  //v_result = _mm256_shuffle_epi8(v_result, v_shuffle2);

  _mm256_store_si256((__m256i*)dst, v_result);
}

static void fast_inverse_tr_8x2_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 8;
  const int height = 2;

  int skip_width = 0;
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = ff_dct2_2xN_coeff_hor; // TODO: rename
  const int16_t* hor_coeff = fi_dct2_2x8_coeff_ver; // rename
  if (hor == DST7) {
    hor_coeff = fi_dst7_2x8_coeff_ver;
  }
  // Only dct2 transform is defined for this block size

  __m256i v_ver_pass_out;
  fast_inverse_tr_8x2_avx2_ver(src, &v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);

  fast_inverse_tr_8x2_avx2_hor(&v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}

static void fast_forward_tr_8x4_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 8;
  const int height = 4;
  
  int skip_width = 0;
  int skip_height = 0;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_8xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_8x4_coeff_ver;
  if (hor == DST7) {
    hor_coeff = ff_dst7_8xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = ff_dct8_8xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = ff_dst7_8x4_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = ff_dct8_8x4_coeff_ver;
  }
  
  __m256i v_hor_pass_out[2];
  fast_forward_tr_8xN_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i     v_shuffle = _mm256_load_si256((const __m256i*)ff_dct2_8x4_ver_pass_shuffle);
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)ff_dct2_8x4_ver_pass_result_shuffle);
  const __m256i*      v_coeff = (const __m256i*)ver_coeff;

  // 32 samples, process in two steps
  __m256i v_src_per_0 = _mm256_permute4x64_epi64(v_hor_pass_out[0], _MM_SHUFFLE(3, 1, 2, 0));
  __m256i v_src_per_1 = _mm256_permute4x64_epi64(v_hor_pass_out[1], _MM_SHUFFLE(3, 1, 2, 0));
  // Weave lo and hi halfs of each 128 bit lane
  __m256i     v_src_0 = _mm256_shuffle_epi8(v_src_per_0, v_shuffle);
  __m256i     v_src_1 = _mm256_shuffle_epi8(v_src_per_1, v_shuffle);

  __m256i   v_madd_00 = _mm256_madd_epi16(v_src_0, v_coeff[0]);
  __m256i   v_madd_01 = _mm256_madd_epi16(v_src_0, v_coeff[2]);
  __m256i   v_madd_02 = _mm256_madd_epi16(v_src_0, v_coeff[4]);
  __m256i   v_madd_03 = _mm256_madd_epi16(v_src_0, v_coeff[6]);

  __m256i   v_madd_10 = _mm256_madd_epi16(v_src_1, v_coeff[1]);
  __m256i   v_madd_11 = _mm256_madd_epi16(v_src_1, v_coeff[3]);
  __m256i   v_madd_12 = _mm256_madd_epi16(v_src_1, v_coeff[5]);
  __m256i   v_madd_13 = _mm256_madd_epi16(v_src_1, v_coeff[7]);

  __m256i     v_add_0 = _mm256_add_epi32(v_madd_00, v_madd_10);
  __m256i     v_add_1 = _mm256_add_epi32(v_madd_01, v_madd_11);
  __m256i     v_add_2 = _mm256_add_epi32(v_madd_02, v_madd_12);
  __m256i     v_add_3 = _mm256_add_epi32(v_madd_03, v_madd_13);

  __m256i   v_trunc_0 = truncate_avx2(v_add_0, debias, shift_2nd);
  __m256i   v_trunc_1 = truncate_avx2(v_add_1, debias, shift_2nd);
  __m256i   v_trunc_2 = truncate_avx2(v_add_2, debias, shift_2nd);
  __m256i   v_trunc_3 = truncate_avx2(v_add_3, debias, shift_2nd);
           
  __m256i  v_result_0 = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
  __m256i  v_result_1 = _mm256_packs_epi32(v_trunc_2, v_trunc_3);

  // Swap each middle 64 bit chunk in both 128 bit lanes
  v_result_0 = _mm256_permute4x64_epi64(v_result_0, _MM_SHUFFLE(3, 1, 2, 0));
  v_result_1 = _mm256_permute4x64_epi64(v_result_1, _MM_SHUFFLE(3, 1, 2, 0));

  // Swap each middle 16 bit value in each 64 bit chunk
  v_result_0 = _mm256_shuffle_epi8(v_result_0, v_res_shuffle);
  v_result_1 = _mm256_shuffle_epi8(v_result_1, v_res_shuffle);

  _mm256_store_si256((__m256i*)dst, v_result_0);
  _mm256_store_si256((__m256i*)(dst + 16), v_result_1);
}


static void fast_inverse_tr_8x4_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;

  const __m256i v_src_raw_0 = _mm256_load_si256((const __m256i*) & src[0]);
  const __m256i v_src_raw_1 = _mm256_load_si256((const __m256i*) & src[16]);

  __m256i v_src_lo = _mm256_unpacklo_epi16(v_src_raw_0, v_src_raw_1);
  __m256i v_src_hi = _mm256_unpackhi_epi16(v_src_raw_0, v_src_raw_1);

  __m256i v_src_0 = _mm256_permute2x128_si256(v_src_lo, v_src_hi, 0x20);
  __m256i v_src_1 = _mm256_permute2x128_si256(v_src_lo, v_src_hi, 0x31);

  __m256i v_madd_00 = _mm256_madd_epi16(v_src_0, v_coeff[0]);
  __m256i v_madd_01 = _mm256_madd_epi16(v_src_0, v_coeff[1]);
  __m256i v_madd_02 = _mm256_madd_epi16(v_src_0, v_coeff[2]);
  __m256i v_madd_03 = _mm256_madd_epi16(v_src_0, v_coeff[3]);

  __m256i v_madd_10 = _mm256_madd_epi16(v_src_1, v_coeff[4]);
  __m256i v_madd_11 = _mm256_madd_epi16(v_src_1, v_coeff[5]);
  __m256i v_madd_12 = _mm256_madd_epi16(v_src_1, v_coeff[6]);
  __m256i v_madd_13 = _mm256_madd_epi16(v_src_1, v_coeff[7]);

  __m256i v_trunc_0 = truncate_avx2(_mm256_add_epi32(v_madd_00, v_madd_10), debias, shift);
  __m256i v_trunc_1 = truncate_avx2(_mm256_add_epi32(v_madd_01, v_madd_11), debias, shift);
  __m256i v_trunc_2 = truncate_avx2(_mm256_add_epi32(v_madd_02, v_madd_12), debias, shift);
  __m256i v_trunc_3 = truncate_avx2(_mm256_add_epi32(v_madd_03, v_madd_13), debias, shift);

  dst[0] = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
  dst[1] = _mm256_packs_epi32(v_trunc_2, v_trunc_3);
}

static void fast_inverse_tr_8x4_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)fi_tr_4x8_result_shuffle_ver);

  __m256i v_src_0 = _mm256_permute2x128_si256(src[0], src[1], 0x20);
  __m256i v_src_1 = _mm256_permute2x128_si256(src[0], src[1], 0x31);

  __m256i v_madd_0[8];
  __m256i v_madd_1[8];
  for (int i = 0; i < 8; ++i) {
    v_madd_0[i] = _mm256_madd_epi16(v_src_0, v_coeff[0]);
    v_madd_1[i] = _mm256_madd_epi16(v_src_1, v_coeff[1]);

    v_coeff += 2;
  }

  __m256i v_add[8];
  for (int i = 0; i < 8; ++i) {
    v_add[i] = _mm256_add_epi32(v_madd_0[i], v_madd_1[i]);
  }

  __m256i v_hadd[4];
  for (int dst = 0, src = 0; dst < 4; ++dst, src += 2) {
    v_hadd[dst] = _mm256_hadd_epi32(v_add[src + 0], v_add[src + 1]);
  }

  __m256i v_trunc[4];
  for (int i = 0; i < 4; ++i) {
    v_trunc[i] = truncate_avx2(v_hadd[i], debias, shift);
  }

  __m256i v_result[2];
  v_result[0] = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
  v_result[1] = _mm256_packs_epi32(v_trunc[2], v_trunc[3]);

  v_result[0] = _mm256_shuffle_epi8(v_result[0], v_res_shuffle);
  v_result[1] = _mm256_shuffle_epi8(v_result[1], v_res_shuffle);

  __m256i v_tmp0 = _mm256_permute2x128_si256(v_result[0], v_result[1], 0x20);
  __m256i v_tmp1 = _mm256_permute2x128_si256(v_result[0], v_result[1], 0x31);

  v_result[0] = _mm256_permute4x64_epi64(v_tmp0, _MM_SHUFFLE(3, 1, 2, 0));
  v_result[1] = _mm256_permute4x64_epi64(v_tmp1, _MM_SHUFFLE(3, 1, 2, 0));

  _mm256_store_si256((__m256i*) & dst[0], v_result[0]);
  _mm256_store_si256((__m256i*) & dst[16], v_result[1]);
}

static void fast_inverse_tr_8x4_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 8;
  const int height = 4;

  int skip_width = 0;
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = fi_dct2_4x8_coeff_hor; // TODO: rename coeff tables
  const int16_t* hor_coeff = fi_dct2_4x8_coeff_ver;
  if (hor == DST7) {
    hor_coeff = fi_dst7_4x8_coeff_ver;
  } else if (hor == DCT8) {
    hor_coeff = fi_dct8_4x8_coeff_ver;
  }
  if (ver == DST7) {
    ver_coeff = fi_dst7_4x8_coeff_hor;
  } else if (ver == DCT8) {
    ver_coeff = fi_dct8_4x8_coeff_hor;
  }

  __m256i v_ver_pass_out[2];
  fast_inverse_tr_8x4_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);
  
  fast_inverse_tr_8x4_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static void fast_forward_tr_8x8_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 8;
  const int height = 8;
  
  int skip_width = 0;
  int skip_height = 0;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_8xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_8x8_coeff_ver;
  if (hor == DST7) {
    hor_coeff = ff_dst7_8xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = ff_dct8_8xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = ff_dst7_8x8_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = ff_dct8_8x8_coeff_ver;
  }

  __m256i v_hor_pass_out[4];
  fast_forward_tr_8xN_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);
  const int32_t* coeff_ptr = (const int32_t*)ver_coeff; // Cast into 32 bit integer to read two coeffs at a time

  __m256i v_src_lo_0 = _mm256_unpacklo_epi16(v_hor_pass_out[0], v_hor_pass_out[1]);
  __m256i v_src_lo_1 = _mm256_unpacklo_epi16(v_hor_pass_out[2], v_hor_pass_out[3]);
  __m256i v_src_hi_0 = _mm256_unpackhi_epi16(v_hor_pass_out[0], v_hor_pass_out[1]);
  __m256i v_src_hi_1 = _mm256_unpackhi_epi16(v_hor_pass_out[2], v_hor_pass_out[3]);

  __m256i v_trunc[8];

  __m256i v_src_0 = _mm256_permute2x128_si256(v_src_lo_0, v_src_hi_0, 0x20);
  __m256i v_src_1 = _mm256_permute2x128_si256(v_src_lo_0, v_src_hi_0, 0x31);
  __m256i v_src_2 = _mm256_permute2x128_si256(v_src_lo_1, v_src_hi_1, 0x20);
  __m256i v_src_3 = _mm256_permute2x128_si256(v_src_lo_1, v_src_hi_1, 0x31);

  for (int i = 0; i < 8; ++i) {
    __m256i v_coeff_0 = _mm256_set1_epi32(coeff_ptr[0]);
    __m256i v_coeff_1 = _mm256_set1_epi32(coeff_ptr[1]);
    __m256i v_coeff_2 = _mm256_set1_epi32(coeff_ptr[2]);
    __m256i v_coeff_3 = _mm256_set1_epi32(coeff_ptr[3]);

    __m256i v_madd_0 = _mm256_madd_epi16(v_src_0, v_coeff_0);
    __m256i v_madd_1 = _mm256_madd_epi16(v_src_1, v_coeff_1);
    __m256i v_madd_2 = _mm256_madd_epi16(v_src_2, v_coeff_2);
    __m256i v_madd_3 = _mm256_madd_epi16(v_src_3, v_coeff_3);

    __m256i v_add_0 = _mm256_add_epi32(v_madd_0, v_madd_1);
    __m256i v_add_1 = _mm256_add_epi32(v_madd_2, v_madd_3);

    v_trunc[i] = truncate_avx2(_mm256_add_epi32(v_add_0, v_add_1), debias, shift_2nd);
    coeff_ptr += 4;
  }

  __m256i v_result[4];
  v_result[0] = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
  v_result[1] = _mm256_packs_epi32(v_trunc[2], v_trunc[3]);
  v_result[2] = _mm256_packs_epi32(v_trunc[4], v_trunc[5]);
  v_result[3] = _mm256_packs_epi32(v_trunc[6], v_trunc[7]);

  for (int i = 0; i < 4; ++i) {
    v_result[i] = _mm256_permute4x64_epi64(v_result[i], _MM_SHUFFLE(3, 1, 2, 0));
  }
  
  for (int i = 0; i < 4; ++i) {
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }
}


static void fast_inverse_tr_8x8_avx2_hor(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i v_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0415);

  const __m256i* v_src_raw = (const __m256i*)src;

  __m256i v_src[4];
  v_src[0] = _mm256_permute4x64_epi64(v_src_raw[0], _MM_SHUFFLE(3, 1, 2, 0));
  v_src[1] = _mm256_permute4x64_epi64(v_src_raw[1], _MM_SHUFFLE(3, 1, 2, 0));
  v_src[2] = _mm256_permute4x64_epi64(v_src_raw[2], _MM_SHUFFLE(3, 1, 2, 0));
  v_src[3] = _mm256_permute4x64_epi64(v_src_raw[3], _MM_SHUFFLE(3, 1, 2, 0));
  
  v_src[0] = _mm256_shuffle_epi8(v_src[0], v_shuffle);
  v_src[1] = _mm256_shuffle_epi8(v_src[1], v_shuffle);
  v_src[2] = _mm256_shuffle_epi8(v_src[2], v_shuffle);
  v_src[3] = _mm256_shuffle_epi8(v_src[3], v_shuffle);

  const __m256i* v_c_ptr = v_coeff;
  __m256i v_madd_0[8];
  __m256i v_madd_1[8];
  __m256i v_madd_2[8];
  __m256i v_madd_3[8];
  for (int i = 0; i < 8; ++i) {
    v_madd_0[i] = _mm256_madd_epi16(v_src[0], v_c_ptr[0]);
    v_madd_1[i] = _mm256_madd_epi16(v_src[1], v_c_ptr[1]);
    v_madd_2[i] = _mm256_madd_epi16(v_src[2], v_c_ptr[2]);
    v_madd_3[i] = _mm256_madd_epi16(v_src[3], v_c_ptr[3]);
    v_c_ptr += 4;
  }

  __m256i v_add[8];
  for (int i = 0; i < 8; ++i) {
    __m256i v_add_0 = _mm256_add_epi32(v_madd_0[i], v_madd_1[i]);
    __m256i v_add_1 = _mm256_add_epi32(v_madd_2[i], v_madd_3[i]);

    v_add[i] = _mm256_add_epi32(v_add_0, v_add_1);
  }
  
  __m256i v_trunc[8];
  for (int i = 0; i < 8; ++i) {
    v_trunc[i] = truncate_avx2(v_add[i], debias, shift);
  }

  dst[0] = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
  dst[1] = _mm256_packs_epi32(v_trunc[2], v_trunc[3]);
  dst[2] = _mm256_packs_epi32(v_trunc[4], v_trunc[5]);
  dst[3] = _mm256_packs_epi32(v_trunc[6], v_trunc[7]);
}

static void fast_inverse_tr_8x8_avx2_ver(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0415);

  __m256i v_src[4];
  v_src[0] = _mm256_shuffle_epi32(src[0], _MM_SHUFFLE(3, 1, 2, 0));
  v_src[1] = _mm256_shuffle_epi32(src[1], _MM_SHUFFLE(3, 1, 2, 0));
  v_src[2] = _mm256_shuffle_epi32(src[2], _MM_SHUFFLE(3, 1, 2, 0));
  v_src[3] = _mm256_shuffle_epi32(src[3], _MM_SHUFFLE(3, 1, 2, 0));

  __m256i v_tmp0 = _mm256_permute2x128_si256(v_src[0], v_src[1], 0x20);
  __m256i v_tmp1 = _mm256_permute2x128_si256(v_src[0], v_src[1], 0x31);
  __m256i v_tmp2 = _mm256_permute2x128_si256(v_src[2], v_src[3], 0x20);
  __m256i v_tmp3 = _mm256_permute2x128_si256(v_src[2], v_src[3], 0x31);

  v_src[0] = _mm256_unpacklo_epi64(v_tmp0, v_tmp2);
  v_src[1] = _mm256_unpackhi_epi64(v_tmp0, v_tmp2);
  v_src[2] = _mm256_unpacklo_epi64(v_tmp1, v_tmp3);
  v_src[3] = _mm256_unpackhi_epi64(v_tmp1, v_tmp3);
  

  const __m256i* v_c_ptr = v_coeff;
  __m256i v_madd_0[8];
  __m256i v_madd_1[8];
  __m256i v_madd_2[8];
  __m256i v_madd_3[8];
  for (int i = 0; i < 8; ++i) {
    v_madd_0[i] = _mm256_madd_epi16(v_src[0], v_c_ptr[0]);
    v_madd_1[i] = _mm256_madd_epi16(v_src[1], v_c_ptr[1]);
    v_madd_2[i] = _mm256_madd_epi16(v_src[2], v_c_ptr[2]);
    v_madd_3[i] = _mm256_madd_epi16(v_src[3], v_c_ptr[3]);
    v_c_ptr += 4;
  }

  __m256i v_add[8];
  for (int i = 0; i < 8; ++i) {
    __m256i v_add_0 = _mm256_add_epi32(v_madd_0[i], v_madd_1[i]);
    __m256i v_add_1 = _mm256_add_epi32(v_madd_2[i], v_madd_3[i]);

    v_add[i] = _mm256_add_epi32(v_add_0, v_add_1);
  }

  __m256i v_trunc[8];
  for (int i = 0; i < 8; ++i) {
    v_trunc[i] = truncate_avx2(v_add[i], debias, shift);
  }

  __m256i v_result[4];
  v_result[0] = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
  v_result[1] = _mm256_packs_epi32(v_trunc[2], v_trunc[3]);
  v_result[2] = _mm256_packs_epi32(v_trunc[4], v_trunc[5]);
  v_result[3] = _mm256_packs_epi32(v_trunc[6], v_trunc[7]);

  v_result[0] = _mm256_shuffle_epi8(v_result[0], v_res_shuffle);
  v_result[1] = _mm256_shuffle_epi8(v_result[1], v_res_shuffle);
  v_result[2] = _mm256_shuffle_epi8(v_result[2], v_res_shuffle);
  v_result[3] = _mm256_shuffle_epi8(v_result[3], v_res_shuffle);

  __m256i v_rtmp0 = _mm256_unpacklo_epi32(v_result[0], v_result[1]);
  __m256i v_rtmp1 = _mm256_unpackhi_epi32(v_result[0], v_result[1]);
  __m256i v_rtmp2 = _mm256_unpacklo_epi32(v_result[2], v_result[3]);
  __m256i v_rtmp3 = _mm256_unpackhi_epi32(v_result[2], v_result[3]);

  __m256i v_tmp20 = _mm256_unpacklo_epi64(v_rtmp0, v_rtmp2);
  __m256i v_tmp21 = _mm256_unpackhi_epi64(v_rtmp0, v_rtmp2);
  __m256i v_tmp22 = _mm256_unpacklo_epi64(v_rtmp1, v_rtmp3);
  __m256i v_tmp23 = _mm256_unpackhi_epi64(v_rtmp1, v_rtmp3);

  v_result[0] = _mm256_permute2x128_si256(v_tmp20, v_tmp21, 0x20);
  v_result[1] = _mm256_permute2x128_si256(v_tmp20, v_tmp21, 0x31);
  v_result[2] = _mm256_permute2x128_si256(v_tmp22, v_tmp23, 0x20);
  v_result[3] = _mm256_permute2x128_si256(v_tmp22, v_tmp23, 0x31);

  for (int i = 0; i < 4; ++i) {
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }
}

static void fast_inverse_tr_8x8_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 8;
  const int height = 8;

  int skip_width = 0;
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* hor_coeff = fi_dct2_8x8_coeff_hor;
  const int16_t* ver_coeff = fi_dct2_8x8_coeff_hor;
  if (hor == DST7) {
    hor_coeff = fi_dst7_8x8_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = fi_dct8_8x8_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = fi_dst7_8x8_coeff_hor;
  } else if (ver == DCT8) {
    ver_coeff = fi_dct8_8x8_coeff_hor;
  }

  __m256i v_hor_pass_out[4];
  fast_inverse_tr_8x8_avx2_hor(src, v_hor_pass_out, ver_coeff, shift_1st, height, 0, skip_width);

  fast_inverse_tr_8x8_avx2_ver(v_hor_pass_out, dst, hor_coeff, shift_2nd, width, skip_width, skip_height);
}


static void fast_forward_tr_8x16_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 8;
  const int height = 16;
  // TODO: might be able to get rid of skips in these tailored solutions
  int skip_width = 0;
  int skip_height = 0; // This is not used anywhere

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_8xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_8x16_coeff_ver;
  if (hor == DST7) {
    hor_coeff = ff_dst7_8xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = ff_dct8_8xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = ff_dst7_8x16_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = ff_dct8_8x16_coeff_ver;
  }

  __m256i v_hor_pass_out[8];
  fast_forward_tr_8xN_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  // Can use same shuffles as 8x4
  const __m256i     v_shuffle = _mm256_load_si256((const __m256i*)ff_dct2_8x4_ver_pass_shuffle);
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)ff_dct2_8x4_ver_pass_result_shuffle);
  //const __m256i* v_coeff = (const __m256i*)ver_coeff;
  const int32_t *line_coeff = (const int32_t*)ver_coeff;

  // Multiply+add all source vectors with coeff vectors
  __m256i v_madd[8][16];
  __m256i* v_src_ptr = v_hor_pass_out;
  for (int i = 0; i < 8; ++i) {
    __m256i v_src_per = _mm256_permute4x64_epi64(v_src_ptr[0], _MM_SHUFFLE(3, 1, 2, 0));
    // Weave lo and hi halfs of each 128 bit lane
    __m256i     v_src = _mm256_shuffle_epi8(v_src_per, v_shuffle);
    
    for (int ii = 0; ii < 16; ++ii) {
      //int coeff_row = ii * 8 + i;
      const int32_t coeff = line_coeff[ii];
      const __m256i v_coeff = _mm256_set1_epi32(coeff);
      v_madd[i][ii] = _mm256_madd_epi16(v_src, v_coeff);
    }
    line_coeff += 16;
    v_src_ptr += 1;
  }

  // Add vectors
  __m256i v_add_0[4][16];
  for (int i = 0; i < 4; ++i) {
    for (int ii = 0; ii < 16; ++ii) {
      int offset = i * 2;
      v_add_0[i][ii] = _mm256_add_epi32(v_madd[offset][ii], v_madd[offset + 1][ii]);
    }
  }
  // Second round of additions
  __m256i v_add_1[2][16];
  for (int i = 0; i < 2; ++i) {
    for (int ii = 0; ii < 16; ++ii) {
      int offset = i * 2;
      v_add_1[i][ii] = _mm256_add_epi32(v_add_0[offset][ii], v_add_0[offset + 1][ii]);
    }
  }
  // Third round of additions
  __m256i v_trunc[16];
  for (int ii = 0; ii < 16; ++ii) {
    v_trunc[ii] = _mm256_add_epi32(v_add_1[0][ii], v_add_1[1][ii]);
    v_trunc[ii] = truncate_avx2(v_trunc[ii], debias, shift_2nd);
  }


  for (int i = 0; i < 16; i += 2) {
    __m256i v_result = _mm256_packs_epi32(v_trunc[i], v_trunc[i + 1]);

    // Swap each middle 64 bit chunk in both 128 bit lanes
    v_result = _mm256_permute4x64_epi64(v_result, _MM_SHUFFLE(3, 1, 2, 0));
    // Swap each middle 16 bit value in each 64 bit chunk
    v_result = _mm256_shuffle_epi8(v_result, v_res_shuffle);
    
    _mm256_store_si256((__m256i*)dst, v_result);
    dst += 16;
  }
}


static void fast_inverse_tr_8x16_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i* v_src_raw = (const __m256i*)src;
  const __m256i v_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0415);

  __m256i v_tmp[8];
  for (int i = 0; i < 8; ++i) {
    v_tmp[i] = _mm256_permute4x64_epi64(v_src_raw[i], _MM_SHUFFLE(3, 1, 2, 0));
  }

  __m256i v_src[8];
  for (int i = 0; i < 8; ++i) {
    v_src[i] = _mm256_shuffle_epi8(v_tmp[i], v_shuffle);
  }

  __m256i v_trunc[16];
  for (int c = 0; c < 16; c++) {
    __m256i v_madd_0 = _mm256_madd_epi16(v_src[0], v_coeff[0]);
    __m256i v_madd_1 = _mm256_madd_epi16(v_src[1], v_coeff[1]);
    __m256i v_madd_2 = _mm256_madd_epi16(v_src[2], v_coeff[2]);
    __m256i v_madd_3 = _mm256_madd_epi16(v_src[3], v_coeff[3]);
    __m256i v_madd_4 = _mm256_madd_epi16(v_src[4], v_coeff[4]);
    __m256i v_madd_5 = _mm256_madd_epi16(v_src[5], v_coeff[5]);
    __m256i v_madd_6 = _mm256_madd_epi16(v_src[6], v_coeff[6]);
    __m256i v_madd_7 = _mm256_madd_epi16(v_src[7], v_coeff[7]);

    v_coeff += 8;

    __m256i v_add_00 = _mm256_add_epi32(v_madd_0, v_madd_1);
    __m256i v_add_01 = _mm256_add_epi32(v_madd_2, v_madd_3);
    __m256i v_add_02 = _mm256_add_epi32(v_madd_4, v_madd_5);
    __m256i v_add_03 = _mm256_add_epi32(v_madd_6, v_madd_7);

    __m256i v_add_10 = _mm256_add_epi32(v_add_00, v_add_01);
    __m256i v_add_11 = _mm256_add_epi32(v_add_02, v_add_03);

    v_trunc[c] = truncate_avx2(_mm256_add_epi32(v_add_10, v_add_11), debias, shift);
  }

  for (int d = 0, s = 0; d < 8; ++d, s += 2) {
    dst[d] = _mm256_packs_epi32(v_trunc[s + 0], v_trunc[s + 1]);
  }
}

static void fast_inverse_tr_8x16_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0415);

  __m256i v_src[8];
  for (int i = 0; i < 8; ++i) {
    v_src[i] = _mm256_shuffle_epi32(src[i], _MM_SHUFFLE(3, 1, 2, 0));
  }

  __m256i v_tmp[8];
  v_tmp[0] = _mm256_permute2x128_si256(v_src[0], v_src[1], 0x20);
  v_tmp[1] = _mm256_permute2x128_si256(v_src[2], v_src[3], 0x20);
  v_tmp[2] = _mm256_permute2x128_si256(v_src[4], v_src[5], 0x20);
  v_tmp[3] = _mm256_permute2x128_si256(v_src[6], v_src[7], 0x20);
  v_tmp[4] = _mm256_permute2x128_si256(v_src[0], v_src[1], 0x31);
  v_tmp[5] = _mm256_permute2x128_si256(v_src[2], v_src[3], 0x31);
  v_tmp[6] = _mm256_permute2x128_si256(v_src[4], v_src[5], 0x31);
  v_tmp[7] = _mm256_permute2x128_si256(v_src[6], v_src[7], 0x31);

  v_src[0] = _mm256_unpacklo_epi32(v_tmp[0], v_tmp[1]);
  v_src[1] = _mm256_unpackhi_epi32(v_tmp[0], v_tmp[1]);
  v_src[2] = _mm256_unpacklo_epi32(v_tmp[4], v_tmp[5]);
  v_src[3] = _mm256_unpackhi_epi32(v_tmp[4], v_tmp[5]);
  v_src[4] = _mm256_unpacklo_epi32(v_tmp[2], v_tmp[3]);
  v_src[5] = _mm256_unpackhi_epi32(v_tmp[2], v_tmp[3]);
  v_src[6] = _mm256_unpacklo_epi32(v_tmp[6], v_tmp[7]);
  v_src[7] = _mm256_unpackhi_epi32(v_tmp[6], v_tmp[7]);

  __m256i v_trunc[2][8];
  for (int d = 0, s = 0; d < 2; ++d, s += 4) {
    const __m256i* v_c_ptr = v_coeff;
    __m256i v_madd_0[8];
    __m256i v_madd_1[8];
    __m256i v_madd_2[8];
    __m256i v_madd_3[8];
    for (int c = 0; c < 8; ++c) {
      v_madd_0[c] = _mm256_madd_epi16(v_src[s + 0], v_c_ptr[0]);
      v_madd_1[c] = _mm256_madd_epi16(v_src[s + 1], v_c_ptr[1]);
      v_madd_2[c] = _mm256_madd_epi16(v_src[s + 2], v_c_ptr[2]);
      v_madd_3[c] = _mm256_madd_epi16(v_src[s + 3], v_c_ptr[3]);
      v_c_ptr += 4;
    }

    for (int i = 0; i < 8; ++i) {
      __m256i v_add_0 = _mm256_add_epi32(v_madd_0[i], v_madd_1[i]);
      __m256i v_add_1 = _mm256_add_epi32(v_madd_2[i], v_madd_3[i]);

      v_trunc[d][i] = truncate_avx2(_mm256_add_epi32(v_add_0, v_add_1), debias, shift);
    }
  }

  __m256i v_rtmp[8];
  v_rtmp[0] = _mm256_packs_epi32(v_trunc[0][0], v_trunc[0][1]);
  v_rtmp[1] = _mm256_packs_epi32(v_trunc[0][2], v_trunc[0][3]);
  v_rtmp[2] = _mm256_packs_epi32(v_trunc[0][4], v_trunc[0][5]);
  v_rtmp[3] = _mm256_packs_epi32(v_trunc[0][6], v_trunc[0][7]);
  v_rtmp[4] = _mm256_packs_epi32(v_trunc[1][0], v_trunc[1][1]);
  v_rtmp[5] = _mm256_packs_epi32(v_trunc[1][2], v_trunc[1][3]);
  v_rtmp[6] = _mm256_packs_epi32(v_trunc[1][4], v_trunc[1][5]);
  v_rtmp[7] = _mm256_packs_epi32(v_trunc[1][6], v_trunc[1][7]);

  for (int i = 0; i < 8; ++i) {
    v_rtmp[i] = _mm256_shuffle_epi8(v_rtmp[i], v_res_shuffle);
  }

  __m256i v_tmp32_lo0 = _mm256_unpacklo_epi32(v_rtmp[0], v_rtmp[1]);
  __m256i v_tmp32_lo1 = _mm256_unpacklo_epi32(v_rtmp[2], v_rtmp[3]);
  __m256i v_tmp32_lo2 = _mm256_unpacklo_epi32(v_rtmp[4], v_rtmp[5]);
  __m256i v_tmp32_lo3 = _mm256_unpacklo_epi32(v_rtmp[6], v_rtmp[7]);

  __m256i v_tmp32_hi0 = _mm256_unpackhi_epi32(v_rtmp[0], v_rtmp[1]);
  __m256i v_tmp32_hi1 = _mm256_unpackhi_epi32(v_rtmp[2], v_rtmp[3]);
  __m256i v_tmp32_hi2 = _mm256_unpackhi_epi32(v_rtmp[4], v_rtmp[5]);
  __m256i v_tmp32_hi3 = _mm256_unpackhi_epi32(v_rtmp[6], v_rtmp[7]);

  __m256i v_tmp64_lo0 = _mm256_unpacklo_epi64(v_tmp32_lo0, v_tmp32_lo1);
  __m256i v_tmp64_lo1 = _mm256_unpacklo_epi64(v_tmp32_hi0, v_tmp32_hi1);
  __m256i v_tmp64_lo2 = _mm256_unpacklo_epi64(v_tmp32_lo2, v_tmp32_lo3);
  __m256i v_tmp64_lo3 = _mm256_unpacklo_epi64(v_tmp32_hi2, v_tmp32_hi3);

  __m256i v_tmp64_hi0 = _mm256_unpackhi_epi64(v_tmp32_lo0, v_tmp32_lo1);
  __m256i v_tmp64_hi1 = _mm256_unpackhi_epi64(v_tmp32_hi0, v_tmp32_hi1);
  __m256i v_tmp64_hi2 = _mm256_unpackhi_epi64(v_tmp32_lo2, v_tmp32_lo3);
  __m256i v_tmp64_hi3 = _mm256_unpackhi_epi64(v_tmp32_hi2, v_tmp32_hi3);

  __m256i v_result[8];
  v_result[0] = _mm256_permute2x128_si256(v_tmp64_lo0, v_tmp64_lo1, 0x20);
  v_result[1] = _mm256_permute2x128_si256(v_tmp64_lo0, v_tmp64_lo1, 0x31);
  v_result[2] = _mm256_permute2x128_si256(v_tmp64_hi0, v_tmp64_hi1, 0x20);
  v_result[3] = _mm256_permute2x128_si256(v_tmp64_hi0, v_tmp64_hi1, 0x31);
  v_result[4] = _mm256_permute2x128_si256(v_tmp64_lo2, v_tmp64_lo3, 0x20);
  v_result[5] = _mm256_permute2x128_si256(v_tmp64_lo2, v_tmp64_lo3, 0x31);
  v_result[6] = _mm256_permute2x128_si256(v_tmp64_hi2, v_tmp64_hi3, 0x20);
  v_result[7] = _mm256_permute2x128_si256(v_tmp64_hi2, v_tmp64_hi3, 0x31);

  for (int i = 0; i < 8; ++i) {
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }
}

static void fast_inverse_tr_8x16_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 8;
  const int height = 16;

  int skip_width = 0;
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = fi_dct2_16x8_coeff_hor; // TODO: rename coeff tables
  const int16_t* hor_coeff = fi_dct2_16x8_coeff_ver;
  if (hor == DST7) {
    hor_coeff = fi_dst7_16x8_coeff_ver;
  } else if (hor == DCT8) {
    hor_coeff = fi_dct8_16x8_coeff_ver;
  }
  if (ver == DST7) {
    ver_coeff = fi_dst7_16x8_coeff_hor;
  } else if (ver == DCT8) {
    ver_coeff = fi_dct8_16x8_coeff_hor;
  }

  __m256i v_ver_pass_out[8];
  fast_inverse_tr_8x16_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);
  
  fast_inverse_tr_8x16_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static void fast_forward_tr_8x32_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 8;
  const int height = 32;

  int skip_width = 0;
  int skip_height = (ver != DCT2 && height == 32) ? 16 : (height > 32 ? height - 32 : 0);

  const int  reduced_line = width - skip_width;
  const int  cutoff = height - skip_height;
  int16_t* p_dst = dst;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_8xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_8x32_coeff_ver;
  if (hor == DST7) {
    hor_coeff = ff_dst7_8xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = ff_dct8_8xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = ff_dst7_8x32_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = ff_dct8_8x32_coeff_ver;
  }

  ALIGNED(32) int16_t v_hor_pass_out[8 * 32];
  fast_forward_tr_8xN_avx2_hor(src, (__m256i *)v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);

  __m256i temp_out[16];
  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  for (int j = 0; j < 8; ++j) {
    __m256i res_0 = _mm256_setzero_si256();
    __m256i res_1 = _mm256_setzero_si256();
    __m256i res_2 = _mm256_setzero_si256();
    __m256i res_3 = _mm256_setzero_si256();
    const int16_t* coeff_start = ver_coeff;
    for (int i = 0; i < 16; ++i) {
      int16_t source[2];
      source[0] = v_hor_pass_out[j + i * 16];
      source[1] = v_hor_pass_out[j + i * 16 + 8];
      int32_t paired_source;
      memcpy(&paired_source, source, sizeof(int32_t));

      __m256i v_src = _mm256_set1_epi32(paired_source);
      __m256i v_coeff_0 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;
      __m256i v_coeff_1 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;
      __m256i v_coeff_2 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;
      __m256i v_coeff_3 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;

      __m256i madd_0 = _mm256_madd_epi16(v_src, v_coeff_0);
      __m256i madd_1 = _mm256_madd_epi16(v_src, v_coeff_1);
      __m256i madd_2 = _mm256_madd_epi16(v_src, v_coeff_2);
      __m256i madd_3 = _mm256_madd_epi16(v_src, v_coeff_3);

      res_0 = _mm256_add_epi32(res_0, madd_0);
      res_1 = _mm256_add_epi32(res_1, madd_1);
      res_2 = _mm256_add_epi32(res_2, madd_2);
      res_3 = _mm256_add_epi32(res_3, madd_3);
    }
    __m256i v_trunc_0 = truncate_avx2(res_0, debias, shift_2nd);
    __m256i v_trunc_1 = truncate_avx2(res_1, debias, shift_2nd);
    __m256i v_trunc_2 = truncate_avx2(res_2, debias, shift_2nd);
    __m256i v_trunc_3 = truncate_avx2(res_3, debias, shift_2nd);

    v_trunc_0 = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
    v_trunc_1 = _mm256_packs_epi32(v_trunc_2, v_trunc_3);
    v_trunc_0 = _mm256_permute4x64_epi64(v_trunc_0, _MM_SHUFFLE(3, 1, 2, 0));
    v_trunc_1 = _mm256_permute4x64_epi64(v_trunc_1, _MM_SHUFFLE(3, 1, 2, 0));
    _mm256_store_si256(temp_out + 2 * j, v_trunc_0);
    _mm256_store_si256(temp_out + 2 * j + 1, v_trunc_1);
  }
  transpose_avx2(temp_out, (__m256i*) dst, 32, 8);
#undef NUM_PARTS
#undef PART_DIMENSION

}


static void fast_inverse_tr_8x32_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const int32_t* c_ptr = (const int32_t*)coeff; // Handle as 32 bit integer to load two coeffs into vector at the same time
  const __m256i* v_src_raw = (const __m256i*)src;

  __m256i v_tmp[16];
  for (int i = 0; i < 16; i += 2) {
    v_tmp[i + 0] = _mm256_permute2x128_si256(v_src_raw[i + 0], v_src_raw[i + 1], 0x20);
    v_tmp[i + 1] = _mm256_permute2x128_si256(v_src_raw[i + 0], v_src_raw[i + 1], 0x31);
  }

  __m256i v_tmp16_lo[8];
  __m256i v_tmp16_hi[8];
  for (int d = 0, s = 0; d < 8; ++d, s += 2) {
    v_tmp16_lo[d] = _mm256_unpacklo_epi16(v_tmp[s + 0], v_tmp[s + 1]);
    v_tmp16_hi[d] = _mm256_unpackhi_epi16(v_tmp[s + 0], v_tmp[s + 1]);
  }

  __m256i v_src[16];
  for (int d = 0, s = 0; d < 16; d += 2, ++s) {
    v_src[d + 0] = _mm256_permute2x128_si256(v_tmp16_lo[s], v_tmp16_hi[s], 0x20);
    v_src[d + 1] = _mm256_permute2x128_si256(v_tmp16_lo[s], v_tmp16_hi[s], 0x31);
  }

  __m256i v_trunc[32];

  for (int row = 0; row < 32; ++row) {
    __m256i v_res = _mm256_setzero_si256();
    for (int i = 0; i < 16; ++i) {
      __m256i v_coeff = _mm256_set1_epi32(*c_ptr);
      __m256i v_madd = _mm256_madd_epi16(v_src[i], v_coeff);
      v_res = _mm256_add_epi32(v_res, v_madd);
      c_ptr++;
    }

    v_trunc[row] = truncate_avx2(v_res, debias, shift);
  }

  for (int d = 0, s = 0; d < 16; ++d, s += 2) {
    dst[d] = _mm256_packs_epi32(v_trunc[s + 0], v_trunc[s + 1]);
  }
}

static void fast_inverse_tr_8x32_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff; // Handle as 32 bit integer to load two coeffs into vector at the same time
  const __m256i* v_src_raw = src;
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0246);

  __m256i v_src[16];
  for (int i = 0; i < 16; i += 2) {
    v_src[i + 0] = _mm256_permute2x128_si256(v_src_raw[i + 0], v_src_raw[i + 1], 0x20);
    v_src[i + 1] = _mm256_permute2x128_si256(v_src_raw[i + 0], v_src_raw[i + 1], 0x31);
  }

  __m256i v_tmp[16];
  for (int s = 0; s < 16; s += 2) {
    __m256i v_add[8];
    for (int d = 0, c = 0; d < 8; ++d, c += 2) {
      __m256i v_madd_0 = _mm256_madd_epi16(v_src[s + 0], v_coeff[c + 0]);
      __m256i v_madd_1 = _mm256_madd_epi16(v_src[s + 1], v_coeff[c + 1]);

      v_add[d] = _mm256_add_epi32(v_madd_0, v_madd_1);
    }

    __m256i v_hadd[4];
    v_hadd[0] = _mm256_hadd_epi32(v_add[0], v_add[1]);
    v_hadd[1] = _mm256_hadd_epi32(v_add[2], v_add[3]);
    v_hadd[2] = _mm256_hadd_epi32(v_add[4], v_add[5]);
    v_hadd[3] = _mm256_hadd_epi32(v_add[6], v_add[7]);

    __m256i v_trunc[4];
    v_trunc[0] = truncate_avx2(v_hadd[0], debias, shift);
    v_trunc[1] = truncate_avx2(v_hadd[1], debias, shift);
    v_trunc[2] = truncate_avx2(v_hadd[2], debias, shift);
    v_trunc[3] = truncate_avx2(v_hadd[3], debias, shift);

    v_tmp[s + 0] = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
    v_tmp[s + 1] = _mm256_packs_epi32(v_trunc[2], v_trunc[3]);
  }

  for (int i = 0; i < 16; ++i) {
    v_tmp[i] = _mm256_shuffle_epi8(v_tmp[i], v_res_shuffle);
  }

  __m256i v_tmp64_lo[8];
  __m256i v_tmp64_hi[8];
  for (int d = 0, s = 0; d < 8; ++d, s += 2) {
    v_tmp64_lo[d] = _mm256_unpacklo_epi64(v_tmp[s + 0], v_tmp[s + 1]);
    v_tmp64_hi[d] = _mm256_unpackhi_epi64(v_tmp[s + 0], v_tmp[s + 1]);
  }

  __m256i v_result[16];
  for (int d = 0, s = 0; d < 16; d += 2, ++s) {
    v_result[d + 0] = _mm256_permute2x128_si256(v_tmp64_lo[s], v_tmp64_hi[s], 0x20);
    v_result[d + 1] = _mm256_permute2x128_si256(v_tmp64_lo[s], v_tmp64_hi[s], 0x31);
  }

  for (int i = 0; i < 16; ++i) {
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }

  // TODO: mts cutoff
}

static void fast_inverse_tr_8x32_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 8;
  const int height = 32;

  int skip_width = 0;
  int skip_height = (ver != DCT2 && height == 32) ? 16 : (height > 32 ? height - 32 : 0);

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = &uvg_g_dct_32_t[0][0];
  const int16_t* hor_coeff = fi_dct2_32x8_coeff_ver; // TODO: rename table
  if (hor == DST7) {
    hor_coeff = fi_dst7_32x8_coeff_ver; // TODO: rename
  } else if (hor == DCT8) {
    hor_coeff = fi_dct8_32x8_coeff_ver; // TODO: rename
  }
  if (ver == DST7) {
    ver_coeff = &uvg_g_dst7_32_t[0][0];
  } else if (ver == DCT8) {
    ver_coeff = &uvg_g_dct8_32[0][0];
  }

  __m256i v_ver_pass_out[16];
  fast_inverse_tr_8x32_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);

  fast_inverse_tr_8x32_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static void fast_forward_DCT2_B16_avx2_hor(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = (shift > 0) ? (1 << (shift - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  // ISP_TODO: might be faster to load these from arrays
  const __m256i v_permute_0 = _mm256_set1_epi32(0);
  const __m256i v_permute_1 = _mm256_set1_epi32(1);
  const __m256i v_permute_2 = _mm256_set1_epi32(2);
  const __m256i v_permute_3 = _mm256_set1_epi32(3);
  const __m256i v_permute_4 = _mm256_set1_epi32(4);
  const __m256i v_permute_5 = _mm256_set1_epi32(5);
  const __m256i v_permute_6 = _mm256_set1_epi32(6);
  const __m256i v_permute_7 = _mm256_set1_epi32(7);

  const __m256i* v_coeff = (const __m256i*)coeff;

  const int reduced_line = line - skip_line;
  // Handle 1 line at a time, 16 samples per line
  for (int j = 0; j < reduced_line; ++j) {
    //                    line 1
    // src vector:       [s00 s01 s02 s03 s04 s05 s06 s07 | s08 s09 s10 s11 s12 s13 s14 s15]
    __m256i v_src_raw = _mm256_load_si256((const __m256i*)src);

    // Arrange data so calculations can be done column-wise (to avoid using hadds).
    // Need 8 source vectors. First will be filled with s00 and s01 pairs. Second with s02 and s03 pairs and so on
    __m256i v_src_0 = _mm256_permutevar8x32_epi32(v_src_raw, v_permute_0);
    __m256i v_src_1 = _mm256_permutevar8x32_epi32(v_src_raw, v_permute_1);
    __m256i v_src_2 = _mm256_permutevar8x32_epi32(v_src_raw, v_permute_2);
    __m256i v_src_3 = _mm256_permutevar8x32_epi32(v_src_raw, v_permute_3);
    __m256i v_src_4 = _mm256_permutevar8x32_epi32(v_src_raw, v_permute_4);
    __m256i v_src_5 = _mm256_permutevar8x32_epi32(v_src_raw, v_permute_5);
    __m256i v_src_6 = _mm256_permutevar8x32_epi32(v_src_raw, v_permute_6);
    __m256i v_src_7 = _mm256_permutevar8x32_epi32(v_src_raw, v_permute_7);

    __m256i v_madd_0_00 = _mm256_madd_epi16(v_src_0, v_coeff[0]);
    __m256i v_madd_0_01 = _mm256_madd_epi16(v_src_1, v_coeff[1]);
    __m256i v_madd_0_02 = _mm256_madd_epi16(v_src_2, v_coeff[2]);
    __m256i v_madd_0_03 = _mm256_madd_epi16(v_src_3, v_coeff[3]);
    __m256i v_madd_0_04 = _mm256_madd_epi16(v_src_4, v_coeff[4]);
    __m256i v_madd_0_05 = _mm256_madd_epi16(v_src_5, v_coeff[5]);
    __m256i v_madd_0_06 = _mm256_madd_epi16(v_src_6, v_coeff[6]);
    __m256i v_madd_0_07 = _mm256_madd_epi16(v_src_7, v_coeff[7]);

    __m256i v_madd_0_08 = _mm256_madd_epi16(v_src_0, v_coeff[8]);
    __m256i v_madd_0_09 = _mm256_madd_epi16(v_src_1, v_coeff[9]);
    __m256i v_madd_0_10 = _mm256_madd_epi16(v_src_2, v_coeff[10]);
    __m256i v_madd_0_11 = _mm256_madd_epi16(v_src_3, v_coeff[11]);
    __m256i v_madd_0_12 = _mm256_madd_epi16(v_src_4, v_coeff[12]);
    __m256i v_madd_0_13 = _mm256_madd_epi16(v_src_5, v_coeff[13]);
    __m256i v_madd_0_14 = _mm256_madd_epi16(v_src_6, v_coeff[14]);
    __m256i v_madd_0_15 = _mm256_madd_epi16(v_src_7, v_coeff[15]);

    __m256i  v_madd_1_0 = _mm256_add_epi32(v_madd_0_00, v_madd_0_01);
    __m256i  v_madd_1_1 = _mm256_add_epi32(v_madd_0_02, v_madd_0_03);
    __m256i  v_madd_1_2 = _mm256_add_epi32(v_madd_0_04, v_madd_0_05);
    __m256i  v_madd_1_3 = _mm256_add_epi32(v_madd_0_06, v_madd_0_07);
    __m256i  v_madd_1_4 = _mm256_add_epi32(v_madd_0_08, v_madd_0_09);
    __m256i  v_madd_1_5 = _mm256_add_epi32(v_madd_0_10, v_madd_0_11);
    __m256i  v_madd_1_6 = _mm256_add_epi32(v_madd_0_12, v_madd_0_13);
    __m256i  v_madd_1_7 = _mm256_add_epi32(v_madd_0_14, v_madd_0_15);

    __m256i  v_madd_2_0 = _mm256_add_epi32(v_madd_1_0, v_madd_1_1);
    __m256i  v_madd_2_1 = _mm256_add_epi32(v_madd_1_2, v_madd_1_3);
    __m256i  v_madd_2_2 = _mm256_add_epi32(v_madd_1_4, v_madd_1_5);
    __m256i  v_madd_2_3 = _mm256_add_epi32(v_madd_1_6, v_madd_1_7);

    __m256i  v_madd_3_0 = _mm256_add_epi32(v_madd_2_0, v_madd_2_1);
    __m256i  v_madd_3_1 = _mm256_add_epi32(v_madd_2_2, v_madd_2_3);

    __m256i v_trunc_0 = truncate_avx2(v_madd_3_0, debias, shift);
    __m256i v_trunc_1 = truncate_avx2(v_madd_3_1, debias, shift);

    __m256i v_result = _mm256_packs_epi32(v_trunc_0, v_trunc_1);

    dst[0] = v_result;
    
    src += 16;
    dst++;
  }
}

static void fast_forward_tr_16x2_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 16;
  const int height = 2;
  // TODO: might be able to get rid of skips in these tailored solutions
  int skip_width = 0;
  int skip_height = 0; // This is not used anywhere

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_16xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_16x2_coeff_ver;
  if (hor == DST7) {
    hor_coeff = ff_dst7_16xN_coeff_hor;    
  }

  __m256i v_hor_pass_out[2];
  fast_forward_DCT2_B16_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)ver_coeff;

  // Got samples for 2 source vectors
  // Unpack -> samples to be added are adjacent
  __m256i v_src_hi = _mm256_unpackhi_epi16(v_hor_pass_out[0], v_hor_pass_out[1]);
  __m256i v_src_lo = _mm256_unpacklo_epi16(v_hor_pass_out[0], v_hor_pass_out[1]);

  __m256i v_madd_hi_0 = _mm256_madd_epi16(v_src_hi, v_coeff[0]);
  __m256i v_madd_hi_1 = _mm256_madd_epi16(v_src_hi, v_coeff[1]);
  __m256i v_madd_lo_0 = _mm256_madd_epi16(v_src_lo, v_coeff[0]);
  __m256i v_madd_lo_1 = _mm256_madd_epi16(v_src_lo, v_coeff[1]);

  __m256i v_trunc_hi_0 = truncate_avx2(v_madd_hi_0, debias, shift_2nd);
  __m256i v_trunc_hi_1 = truncate_avx2(v_madd_hi_1, debias, shift_2nd);
  __m256i v_trunc_lo_0 = truncate_avx2(v_madd_lo_0, debias, shift_2nd);
  __m256i v_trunc_lo_1 = truncate_avx2(v_madd_lo_1, debias, shift_2nd);

  __m256i v_result_0 = _mm256_packs_epi32(v_trunc_lo_0, v_trunc_hi_0);
  __m256i v_result_1 = _mm256_packs_epi32(v_trunc_lo_1, v_trunc_hi_1);

  _mm256_store_si256((__m256i*)dst, v_result_0);
  _mm256_store_si256((__m256i*)(dst + 16), v_result_1);
}


static void fast_inverse_tr_16x2_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i v_coeff_0 = _mm256_load_si256((const __m256i*) & coeff[0]);
  const __m256i v_coeff_1 = _mm256_load_si256((const __m256i*) & coeff[16]);

  const __m256i v_src_0 = _mm256_load_si256((const __m256i*) & src[0]);
  const __m256i v_src_1 = _mm256_load_si256((const __m256i*) & src[16]);

  const __m256i v_src_lo = _mm256_unpacklo_epi16(v_src_0, v_src_1);
  const __m256i v_src_hi = _mm256_unpackhi_epi16(v_src_0, v_src_1);

  __m256i v_trunc_0 = truncate_avx2(_mm256_madd_epi16(v_src_lo, v_coeff_0), debias, shift);
  __m256i v_trunc_1 = truncate_avx2(_mm256_madd_epi16(v_src_lo, v_coeff_1), debias, shift);
  __m256i v_trunc_2 = truncate_avx2(_mm256_madd_epi16(v_src_hi, v_coeff_0), debias, shift);
  __m256i v_trunc_3 = truncate_avx2(_mm256_madd_epi16(v_src_hi, v_coeff_1), debias, shift);

  dst[0] = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
  dst[1] = _mm256_packs_epi32(v_trunc_2, v_trunc_3);
}

static void fast_inverse_tr_16x2_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = (shift > 0) ? (1 << (shift - 1)) : 0;
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;

  __m256i v_madd_e[16];
  __m256i v_madd_o[16];
  for (int i = 0, c = 0; i < 16; ++i, c += 2) {
    v_madd_e[i] = _mm256_madd_epi16(src[0], v_coeff[c + 0]);
    v_madd_o[i] = _mm256_madd_epi16(src[1], v_coeff[c + 1]);
  }

  __m256i v_add[16];
  for (int i = 0; i < 16; ++i) {
    v_add[i] = _mm256_add_epi32(v_madd_e[i], v_madd_o[i]);
  }

  for (int i = 0; i < 16; ++i) {
    v_add[i] = _mm256_permute4x64_epi64(v_add[i], _MM_SHUFFLE(3, 1, 2, 0));
  }

  __m256i v_hadd_0[8];
  for (int src = 0, dst = 0; dst < 8; ++dst, src += 2) {
    v_hadd_0[dst] = _mm256_hadd_epi32(v_add[src + 0], v_add[src + 1]);
  }

  __m256i v_trunc[4];
  for (int src = 0, dst = 0; dst < 4; ++dst, src += 2) {
    v_trunc[dst] = truncate_avx2(_mm256_hadd_epi32(v_hadd_0[src + 0], v_hadd_0[src + 1]), debias, shift);
  }

  __m256i v_tmp0 = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
  __m256i v_tmp1 = _mm256_packs_epi32(v_trunc[2], v_trunc[3]);

  __m256i v_result_0 = _mm256_permute2x128_si256(v_tmp0, v_tmp1, 0x20);
  __m256i v_result_1 = _mm256_permute2x128_si256(v_tmp0, v_tmp1, 0x31);

  _mm256_store_si256((__m256i*) & dst[0], v_result_0);
  _mm256_store_si256((__m256i*) & dst[16], v_result_1);
}

static void fast_inverse_tr_16x2_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 16;
  const int height = 2;

  int skip_width = 0;
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = ff_dct2_2xN_coeff_hor; // TODO: rename
  const int16_t* hor_coeff = fi_dct2_2x16_coeff_ver; // rename
  if (hor == DST7) {
    hor_coeff = fi_dst7_2x16_coeff_ver;
  }
  // DST7 and DCT8 are not defined for this block size

  __m256i v_ver_pass_out[2];
  fast_inverse_tr_16x2_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);

  fast_inverse_tr_16x2_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static void fast_forward_tr_16x4_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 16;
  const int height = 4;
  // TODO: might be able to get rid of skips in these tailored solutions
  int skip_width = 0;
  int skip_height = 0; // This is not used anywhere

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_16xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_16x4_coeff_ver;
  if (hor == DST7) {
    hor_coeff = ff_dst7_16xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = ff_dct8_16xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = ff_dst7_16x4_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = ff_dct8_16x4_coeff_ver;
  }

  __m256i v_hor_pass_out[4];
  fast_forward_DCT2_B16_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)ver_coeff;

  // Got samples for 4 vectors
  __m256i v_src_hi_0 = _mm256_unpackhi_epi16(v_hor_pass_out[0], v_hor_pass_out[1]);
  __m256i v_src_hi_1 = _mm256_unpackhi_epi16(v_hor_pass_out[2], v_hor_pass_out[3]);
  __m256i v_src_lo_0 = _mm256_unpacklo_epi16(v_hor_pass_out[0], v_hor_pass_out[1]);
  __m256i v_src_lo_1 = _mm256_unpacklo_epi16(v_hor_pass_out[2], v_hor_pass_out[3]);

  __m256i v_madd_hi_00 = _mm256_madd_epi16(v_src_hi_0, v_coeff[0]);
  __m256i v_madd_hi_01 = _mm256_madd_epi16(v_src_hi_0, v_coeff[2]);
  __m256i v_madd_hi_02 = _mm256_madd_epi16(v_src_hi_0, v_coeff[4]);
  __m256i v_madd_hi_03 = _mm256_madd_epi16(v_src_hi_0, v_coeff[6]);
  __m256i v_madd_hi_10 = _mm256_madd_epi16(v_src_hi_1, v_coeff[1]);
  __m256i v_madd_hi_11 = _mm256_madd_epi16(v_src_hi_1, v_coeff[3]);
  __m256i v_madd_hi_12 = _mm256_madd_epi16(v_src_hi_1, v_coeff[5]);
  __m256i v_madd_hi_13 = _mm256_madd_epi16(v_src_hi_1, v_coeff[7]);

  __m256i v_madd_lo_00 = _mm256_madd_epi16(v_src_lo_0, v_coeff[0]);
  __m256i v_madd_lo_01 = _mm256_madd_epi16(v_src_lo_0, v_coeff[2]);
  __m256i v_madd_lo_02 = _mm256_madd_epi16(v_src_lo_0, v_coeff[4]);
  __m256i v_madd_lo_03 = _mm256_madd_epi16(v_src_lo_0, v_coeff[6]);
  __m256i v_madd_lo_10 = _mm256_madd_epi16(v_src_lo_1, v_coeff[1]);
  __m256i v_madd_lo_11 = _mm256_madd_epi16(v_src_lo_1, v_coeff[3]);
  __m256i v_madd_lo_12 = _mm256_madd_epi16(v_src_lo_1, v_coeff[5]);
  __m256i v_madd_lo_13 = _mm256_madd_epi16(v_src_lo_1, v_coeff[7]);

  __m256i v_add_hi_0 = _mm256_add_epi32(v_madd_hi_00, v_madd_hi_10);
  __m256i v_add_hi_1 = _mm256_add_epi32(v_madd_hi_01, v_madd_hi_11);
  __m256i v_add_hi_2 = _mm256_add_epi32(v_madd_hi_02, v_madd_hi_12);
  __m256i v_add_hi_3 = _mm256_add_epi32(v_madd_hi_03, v_madd_hi_13);

  __m256i v_add_lo_0 = _mm256_add_epi32(v_madd_lo_00, v_madd_lo_10);
  __m256i v_add_lo_1 = _mm256_add_epi32(v_madd_lo_01, v_madd_lo_11);
  __m256i v_add_lo_2 = _mm256_add_epi32(v_madd_lo_02, v_madd_lo_12);
  __m256i v_add_lo_3 = _mm256_add_epi32(v_madd_lo_03, v_madd_lo_13);

  __m256i v_trunc_hi_0 = truncate_avx2(v_add_hi_0, debias, shift_2nd);
  __m256i v_trunc_hi_1 = truncate_avx2(v_add_hi_1, debias, shift_2nd);
  __m256i v_trunc_hi_2 = truncate_avx2(v_add_hi_2, debias, shift_2nd);
  __m256i v_trunc_hi_3 = truncate_avx2(v_add_hi_3, debias, shift_2nd);

  __m256i v_trunc_lo_0 = truncate_avx2(v_add_lo_0, debias, shift_2nd);
  __m256i v_trunc_lo_1 = truncate_avx2(v_add_lo_1, debias, shift_2nd);
  __m256i v_trunc_lo_2 = truncate_avx2(v_add_lo_2, debias, shift_2nd);
  __m256i v_trunc_lo_3 = truncate_avx2(v_add_lo_3, debias, shift_2nd);
  
  __m256i v_result_0 = _mm256_packs_epi32(v_trunc_lo_0, v_trunc_hi_0);
  __m256i v_result_1 = _mm256_packs_epi32(v_trunc_lo_1, v_trunc_hi_1);
  __m256i v_result_2 = _mm256_packs_epi32(v_trunc_lo_2, v_trunc_hi_2);
  __m256i v_result_3 = _mm256_packs_epi32(v_trunc_lo_3, v_trunc_hi_3);

  _mm256_store_si256((__m256i*)dst, v_result_0);
  _mm256_store_si256((__m256i*)(dst + 16), v_result_1);
  _mm256_store_si256((__m256i*)(dst + 32), v_result_2);
  _mm256_store_si256((__m256i*)(dst + 48), v_result_3);
}


static void fast_inverse_tr_16x4_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;

  const __m256i* v_src_raw = (const __m256i*)src;

  __m256i v_src_lo_0 = _mm256_unpacklo_epi16(v_src_raw[0], v_src_raw[1]);
  __m256i v_src_lo_1 = _mm256_unpacklo_epi16(v_src_raw[2], v_src_raw[3]);
  __m256i v_src_hi_0 = _mm256_unpackhi_epi16(v_src_raw[0], v_src_raw[1]);
  __m256i v_src_hi_1 = _mm256_unpackhi_epi16(v_src_raw[2], v_src_raw[3]);

  __m256i v_madd_lo_0[4];
  __m256i v_madd_lo_1[4];
  __m256i v_madd_hi_0[4];
  __m256i v_madd_hi_1[4];
  for (int i = 0; i < 4; i++) {
    v_madd_lo_0[i] = _mm256_madd_epi16(v_src_lo_0, v_coeff[0]);
    v_madd_lo_1[i] = _mm256_madd_epi16(v_src_lo_1, v_coeff[1]);

    v_madd_hi_0[i] = _mm256_madd_epi16(v_src_hi_0, v_coeff[0]);
    v_madd_hi_1[i] = _mm256_madd_epi16(v_src_hi_1, v_coeff[1]);

    v_coeff += 2;
  }

  __m256i v_trunc_lo[4];
  __m256i v_trunc_hi[4];
  for (int i = 0; i < 4; ++i) {
    v_trunc_lo[i] = truncate_avx2(_mm256_add_epi32(v_madd_lo_0[i], v_madd_lo_1[i]), debias, shift);
    v_trunc_hi[i] = truncate_avx2(_mm256_add_epi32(v_madd_hi_0[i], v_madd_hi_1[i]), debias, shift);
  }

  dst[0] = _mm256_packs_epi32(v_trunc_lo[0], v_trunc_hi[0]);
  dst[1] = _mm256_packs_epi32(v_trunc_lo[1], v_trunc_hi[1]);
  dst[2] = _mm256_packs_epi32(v_trunc_lo[2], v_trunc_hi[2]);
  dst[3] = _mm256_packs_epi32(v_trunc_lo[3], v_trunc_hi[3]);
}

static void fast_inverse_tr_16x4_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)fi_tr_4x8_result_shuffle_ver); // Can use existing shuffle vector

  __m256i v_src_0 = _mm256_permute2x128_si256(src[0], src[2], 0x20);
  __m256i v_src_1 = _mm256_permute2x128_si256(src[0], src[2], 0x31);
  __m256i v_src_2 = _mm256_permute2x128_si256(src[1], src[3], 0x20);
  __m256i v_src_3 = _mm256_permute2x128_si256(src[1], src[3], 0x31);

  __m256i v_madd_0[16];
  __m256i v_madd_1[16];
  __m256i v_madd_2[16];
  __m256i v_madd_3[16];
  for (int i = 0; i < 16; ++i) {
    v_madd_0[i] = _mm256_madd_epi16(v_src_0, v_coeff[0]);
    v_madd_1[i] = _mm256_madd_epi16(v_src_1, v_coeff[1]);
    v_madd_2[i] = _mm256_madd_epi16(v_src_2, v_coeff[0]);
    v_madd_3[i] = _mm256_madd_epi16(v_src_3, v_coeff[1]);

    v_coeff += 2;
  }

  __m256i v_add_0[16];
  __m256i v_add_1[16];
  for (int i = 0; i < 16; ++i) {
    v_add_0[i] = _mm256_add_epi32(v_madd_0[i], v_madd_1[i]);
    v_add_1[i] = _mm256_add_epi32(v_madd_2[i], v_madd_3[i]);

  }

  __m256i v_hadd_0[16];
  for (int i = 0; i < 16; ++i) {
    v_hadd_0[i] = _mm256_hadd_epi32(v_add_0[i], v_add_1[i]);
  }

  __m256i v_hadd_1[8];
  for (int dst = 0, src = 0; dst < 8; ++dst, src += 2) {
    v_hadd_1[dst] = _mm256_hadd_epi32(v_hadd_0[src + 0], v_hadd_0[src + 1]);
  }

  __m256i v_trunc[8];
  for (int i = 0; i < 8; ++i) {
    v_trunc[i] = truncate_avx2(v_hadd_1[i], debias, shift);
  }

  __m256i v_result[4];
  __m256i v_tmp0 = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
  __m256i v_tmp1 = _mm256_packs_epi32(v_trunc[2], v_trunc[3]);
  __m256i v_tmp2 = _mm256_packs_epi32(v_trunc[4], v_trunc[5]);
  __m256i v_tmp3 = _mm256_packs_epi32(v_trunc[6], v_trunc[7]);

  v_tmp0 = _mm256_shuffle_epi8(v_tmp0, v_res_shuffle);
  v_tmp1 = _mm256_shuffle_epi8(v_tmp1, v_res_shuffle);
  v_tmp2 = _mm256_shuffle_epi8(v_tmp2, v_res_shuffle);
  v_tmp3 = _mm256_shuffle_epi8(v_tmp3, v_res_shuffle);

  __m256i v_tmp_lo_0 = _mm256_unpacklo_epi64(v_tmp0, v_tmp1);
  __m256i v_tmp_lo_1 = _mm256_unpacklo_epi64(v_tmp2, v_tmp3);
  __m256i v_tmp_hi_0 = _mm256_unpackhi_epi64(v_tmp0, v_tmp1);
  __m256i v_tmp_hi_1 = _mm256_unpackhi_epi64(v_tmp2, v_tmp3);

  v_result[0] = _mm256_permute2x128_si256(v_tmp_lo_0, v_tmp_lo_1, 0x20);
  v_result[1] = _mm256_permute2x128_si256(v_tmp_hi_0, v_tmp_hi_1, 0x20);
  v_result[2] = _mm256_permute2x128_si256(v_tmp_lo_0, v_tmp_lo_1, 0x31);
  v_result[3] = _mm256_permute2x128_si256(v_tmp_hi_0, v_tmp_hi_1, 0x31);

  for (int i = 0; i < 4; ++i) {
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }
}

static void fast_inverse_tr_16x4_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 16;
  const int height = 4;

  int skip_width = 0;
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = fi_dct2_4x16_coeff_hor; // TODO: rename coeff tables
  const int16_t* hor_coeff = fi_dct2_4x16_coeff_ver;
  if (hor == DST7) {
    hor_coeff = fi_dst7_4x16_coeff_ver;
  } else if (hor == DCT8) {
    hor_coeff = fi_dct8_4x16_coeff_ver;
  }
  if (ver == DST7) {
    ver_coeff = fi_dst7_4x16_coeff_hor;
  } else if (ver == DCT8) {
    ver_coeff = fi_dct8_4x16_coeff_hor;
  }

  __m256i v_ver_pass_out[4];
  fast_inverse_tr_16x4_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);
  
  fast_inverse_tr_16x4_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static void fast_forward_tr_16x8_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 16;
  const int height = 8;
  
  int skip_width = 0;
  int skip_height = 0;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_16xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_16x8_coeff_ver;
  if (hor == DST7) {
    hor_coeff = ff_dst7_16xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = ff_dct8_16xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = ff_dst7_16x8_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = ff_dct8_16x8_coeff_ver;
  }

  __m256i v_hor_pass_out[8];
  fast_forward_DCT2_B16_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  const int32_t* line_coeff = (const int32_t*)ver_coeff;

  // Got 8 lines of samples. Handle two lines at a time (beacuse of unpack)
  __m256i v_madd_hi[4][8];
  __m256i v_madd_lo[4][8];
  __m256i* v_src_ptr = v_hor_pass_out;
  for (int i = 0; i < 4; ++i) {
    __m256i v_src_hi = _mm256_unpackhi_epi16(v_src_ptr[0], v_src_ptr[1]);
    __m256i v_src_lo = _mm256_unpacklo_epi16(v_src_ptr[0], v_src_ptr[1]);

    // Apply coefficients
    for (int ii = 0; ii < 8; ++ii) {
      const int32_t coeff = line_coeff[ii];
      const __m256i v_coeff = _mm256_set1_epi32(coeff);
      v_madd_hi[i][ii] = _mm256_madd_epi16(v_src_hi, v_coeff);
      v_madd_lo[i][ii] = _mm256_madd_epi16(v_src_lo, v_coeff);
    }

    line_coeff += 8;
    v_src_ptr += 2;
  }

  // First round of additions
  __m256i v_add_hi[2][8];
  __m256i v_add_lo[2][8];
  for (int i = 0; i < 2; ++i) {
    for (int ii = 0; ii < 8; ++ii) {
      const int offset = i * 2;
      v_add_hi[i][ii] = _mm256_add_epi32(v_madd_hi[offset][ii], v_madd_hi[offset + 1][ii]);
      v_add_lo[i][ii] = _mm256_add_epi32(v_madd_lo[offset][ii], v_madd_lo[offset + 1][ii]);
    }
  }

  // Final round of additions, truncation and store
  for (int ii = 0; ii < 8; ++ii) {
    __m256i v_trunc_hi = truncate_avx2(_mm256_add_epi32(v_add_hi[0][ii], v_add_hi[1][ii]), debias, shift_2nd);
    __m256i v_trunc_lo = truncate_avx2(_mm256_add_epi32(v_add_lo[0][ii], v_add_lo[1][ii]), debias, shift_2nd);
    __m256i v_result = _mm256_packs_epi32(v_trunc_lo, v_trunc_hi);

    _mm256_store_si256((__m256i*)dst, v_result);
    dst += 16;
  }
}


static void fast_inverse_tr_16x8_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i* v_src_raw = (const __m256i*)src;

  __m256i v_src_lo[4];
  __m256i v_src_hi[4];
  for (int dst = 0, src = 0; dst < 4; ++dst, src += 2) {
    v_src_lo[dst] = _mm256_unpacklo_epi16(v_src_raw[src + 0], v_src_raw[src + 1]);
    v_src_hi[dst] = _mm256_unpackhi_epi16(v_src_raw[src + 0], v_src_raw[src + 1]);
  }

  __m256i v_trunc_lo[8];
  __m256i v_trunc_hi[8];

  for (int c = 0; c < 8; c++) {
    __m256i v_madd_lo[4];
    __m256i v_madd_hi[4];
    for (int i = 0; i < 4; ++i) {
      v_madd_lo[i] = _mm256_madd_epi16(v_src_lo[i], v_coeff[i]);
      v_madd_hi[i] = _mm256_madd_epi16(v_src_hi[i], v_coeff[i]);
    }
    v_coeff += 4;

    __m256i v_add_lo_0 = _mm256_add_epi32(v_madd_lo[0], v_madd_lo[1]);
    __m256i v_add_lo_1 = _mm256_add_epi32(v_madd_lo[2], v_madd_lo[3]);

    __m256i v_add_hi_0 = _mm256_add_epi32(v_madd_hi[0], v_madd_hi[1]);
    __m256i v_add_hi_1 = _mm256_add_epi32(v_madd_hi[2], v_madd_hi[3]);

    v_trunc_lo[c] = truncate_avx2(_mm256_add_epi32(v_add_lo_0, v_add_lo_1), debias, shift);
    v_trunc_hi[c] = truncate_avx2(_mm256_add_epi32(v_add_hi_0, v_add_hi_1), debias, shift);
  }

  for (int i = 0; i < 8; ++i) {
    dst[i] = _mm256_packs_epi32(v_trunc_lo[i], v_trunc_hi[i]);
  }
}

static void fast_inverse_tr_16x8_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0415);

  __m256i v_tmp32_lo_0 = _mm256_unpacklo_epi32(src[0], src[1]);
  __m256i v_tmp32_lo_1 = _mm256_unpacklo_epi32(src[2], src[3]);
  __m256i v_tmp32_lo_2 = _mm256_unpacklo_epi32(src[4], src[5]);
  __m256i v_tmp32_lo_3 = _mm256_unpacklo_epi32(src[6], src[7]);

  __m256i v_tmp32_hi_0 = _mm256_unpackhi_epi32(src[0], src[1]);
  __m256i v_tmp32_hi_1 = _mm256_unpackhi_epi32(src[2], src[3]);
  __m256i v_tmp32_hi_2 = _mm256_unpackhi_epi32(src[4], src[5]);
  __m256i v_tmp32_hi_3 = _mm256_unpackhi_epi32(src[6], src[7]);

  __m256i v_tmp64_lo_0 = _mm256_unpacklo_epi64(v_tmp32_lo_0, v_tmp32_lo_1);
  __m256i v_tmp64_lo_1 = _mm256_unpacklo_epi64(v_tmp32_lo_2, v_tmp32_lo_3);
  __m256i v_tmp64_lo_2 = _mm256_unpacklo_epi64(v_tmp32_hi_0, v_tmp32_hi_1);
  __m256i v_tmp64_lo_3 = _mm256_unpacklo_epi64(v_tmp32_hi_2, v_tmp32_hi_3);

  __m256i v_tmp64_hi_0 = _mm256_unpackhi_epi64(v_tmp32_lo_0, v_tmp32_lo_1);
  __m256i v_tmp64_hi_1 = _mm256_unpackhi_epi64(v_tmp32_lo_2, v_tmp32_lo_3);
  __m256i v_tmp64_hi_2 = _mm256_unpackhi_epi64(v_tmp32_hi_0, v_tmp32_hi_1);
  __m256i v_tmp64_hi_3 = _mm256_unpackhi_epi64(v_tmp32_hi_2, v_tmp32_hi_3);

  __m256i v_src[8];
  v_src[0] = _mm256_permute2x128_si256(v_tmp64_lo_0, v_tmp64_lo_1, 0x20);
  v_src[1] = _mm256_permute2x128_si256(v_tmp64_hi_0, v_tmp64_hi_1, 0x20);
  v_src[2] = _mm256_permute2x128_si256(v_tmp64_lo_2, v_tmp64_lo_3, 0x20);
  v_src[3] = _mm256_permute2x128_si256(v_tmp64_hi_2, v_tmp64_hi_3, 0x20);
  v_src[4] = _mm256_permute2x128_si256(v_tmp64_lo_0, v_tmp64_lo_1, 0x31);
  v_src[5] = _mm256_permute2x128_si256(v_tmp64_hi_0, v_tmp64_hi_1, 0x31);
  v_src[6] = _mm256_permute2x128_si256(v_tmp64_lo_2, v_tmp64_lo_3, 0x31);
  v_src[7] = _mm256_permute2x128_si256(v_tmp64_hi_2, v_tmp64_hi_3, 0x31);


  __m256i v_trunc[16];
  for (int c = 0; c < 16; ++c) {
    __m256i v_madd[8];
    for (int i = 0; i < 8; ++i) {
      v_madd[i] = _mm256_madd_epi16(v_src[i], v_coeff[i]);
    }
    v_coeff += 8;

    __m256i v_add_0[4];
    for (int dst = 0, src = 0; dst < 4; ++dst, src += 2) {
      v_add_0[dst] = _mm256_add_epi32(v_madd[src + 0], v_madd[src + 1]);
    }

    __m256i v_add_10 = _mm256_add_epi32(v_add_0[0], v_add_0[1]);
    __m256i v_add_11 = _mm256_add_epi32(v_add_0[2], v_add_0[3]);

    v_trunc[c] = truncate_avx2(_mm256_add_epi32(v_add_10, v_add_11), debias, shift);
  }

  __m256i v_result[8];
  for (int dst = 0, src = 0; dst < 8; ++dst, src += 2) {
    v_result[dst] = _mm256_packs_epi32(v_trunc[src + 0], v_trunc[src + 1]);
  }

  for (int i = 0; i < 8; ++i) {
    v_result[i] = _mm256_shuffle_epi8(v_result[i], v_res_shuffle);
  }

  __m256i v_rtmp32_lo_0 = _mm256_unpacklo_epi32(v_result[0], v_result[1]);
  __m256i v_rtmp32_lo_1 = _mm256_unpacklo_epi32(v_result[2], v_result[3]);
  __m256i v_rtmp32_lo_2 = _mm256_unpacklo_epi32(v_result[4], v_result[5]);
  __m256i v_rtmp32_lo_3 = _mm256_unpacklo_epi32(v_result[6], v_result[7]);

  __m256i v_rtmp32_hi_0 = _mm256_unpackhi_epi32(v_result[0], v_result[1]);
  __m256i v_rtmp32_hi_1 = _mm256_unpackhi_epi32(v_result[2], v_result[3]);
  __m256i v_rtmp32_hi_2 = _mm256_unpackhi_epi32(v_result[4], v_result[5]);
  __m256i v_rtmp32_hi_3 = _mm256_unpackhi_epi32(v_result[6], v_result[7]);

  __m256i v_rtmp64_lo_0 = _mm256_unpacklo_epi64(v_rtmp32_lo_0, v_rtmp32_lo_1);
  __m256i v_rtmp64_lo_1 = _mm256_unpacklo_epi64(v_rtmp32_lo_2, v_rtmp32_lo_3);
  __m256i v_rtmp64_lo_2 = _mm256_unpacklo_epi64(v_rtmp32_hi_0, v_rtmp32_hi_1);
  __m256i v_rtmp64_lo_3 = _mm256_unpacklo_epi64(v_rtmp32_hi_2, v_rtmp32_hi_3);

  __m256i v_rtmp64_hi_0 = _mm256_unpackhi_epi64(v_rtmp32_lo_0, v_rtmp32_lo_1);
  __m256i v_rtmp64_hi_1 = _mm256_unpackhi_epi64(v_rtmp32_lo_2, v_rtmp32_lo_3);
  __m256i v_rtmp64_hi_2 = _mm256_unpackhi_epi64(v_rtmp32_hi_0, v_rtmp32_hi_1);
  __m256i v_rtmp64_hi_3 = _mm256_unpackhi_epi64(v_rtmp32_hi_2, v_rtmp32_hi_3);

  v_result[0] = _mm256_permute2x128_si256(v_rtmp64_lo_0, v_rtmp64_lo_1, 0x20);
  v_result[1] = _mm256_permute2x128_si256(v_rtmp64_hi_0, v_rtmp64_hi_1, 0x20);
  v_result[2] = _mm256_permute2x128_si256(v_rtmp64_lo_2, v_rtmp64_lo_3, 0x20);
  v_result[3] = _mm256_permute2x128_si256(v_rtmp64_hi_2, v_rtmp64_hi_3, 0x20);

  v_result[4] = _mm256_permute2x128_si256(v_rtmp64_lo_0, v_rtmp64_lo_1, 0x31);
  v_result[5] = _mm256_permute2x128_si256(v_rtmp64_hi_0, v_rtmp64_hi_1, 0x31);
  v_result[6] = _mm256_permute2x128_si256(v_rtmp64_lo_2, v_rtmp64_lo_3, 0x31);
  v_result[7] = _mm256_permute2x128_si256(v_rtmp64_hi_2, v_rtmp64_hi_3, 0x31);

  for (int i = 0; i < 8; ++i) {
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }
}

static void fast_inverse_tr_16x8_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 16;
  const int height = 8;

  int skip_width = 0;
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = fi_dct2_8x16_coeff_hor;
  const int16_t* hor_coeff = fi_dct2_8x16_coeff_ver;
  if (hor == DST7) {
    hor_coeff = fi_dst7_8x16_coeff_ver;
  } else if (hor == DCT8) {
    hor_coeff = fi_dct8_8x16_coeff_ver;
  }
  if (ver == DST7) {
    ver_coeff = fi_dst7_8x16_coeff_hor;
  } else if (ver == DCT8) {
    ver_coeff = fi_dct8_8x16_coeff_hor;
  }

  __m256i v_ver_pass_out[8];
  fast_inverse_tr_16x8_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);
  
  fast_inverse_tr_16x8_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static void fast_forward_tr_16x16_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 16;
  const int height = 16;
  
  int skip_width = 0;
  int skip_height = 0; 

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_16xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_16x16_coeff_ver;
  if (hor == DST7) {
    hor_coeff = ff_dst7_16xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = ff_dct8_16xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = ff_dst7_16x16_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = ff_dct8_16x16_coeff_ver;
  }

  __m256i v_hor_pass_out[16];
  fast_forward_DCT2_B16_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

#define NUM_PARTS 4
#define PART_DIMENSION (16 / NUM_PARTS)
  for (int part = 0; part < NUM_PARTS; ++part) {
    const int32_t* coeff_ptr = (const int32_t*)ver_coeff + part * PART_DIMENSION; // Cast into 32 bit integer to read two coeffs at a time
    const __m256i* v_src_ptr = v_hor_pass_out;

    __m256i v_madd_lo[8][PART_DIMENSION];
    __m256i v_madd_hi[8][PART_DIMENSION];
    for (int i = 0; i < 8; ++i) {
      __m256i v_src_lo = _mm256_unpacklo_epi16(v_src_ptr[0], v_src_ptr[1]);
      __m256i v_src_hi = _mm256_unpackhi_epi16(v_src_ptr[0], v_src_ptr[1]);

      for (int c = 0; c < PART_DIMENSION; ++c) {
        const __m256i v_coeff = _mm256_set1_epi32(coeff_ptr[c]);
        v_madd_lo[i][c] = _mm256_madd_epi16(v_src_lo, v_coeff);
        v_madd_hi[i][c] = _mm256_madd_epi16(v_src_hi, v_coeff);
      }
      v_src_ptr += 2;
      coeff_ptr += 16;
    }

    __m256i v_trunc_lo[PART_DIMENSION];
    __m256i v_trunc_hi[PART_DIMENSION];
    for (int i = 0; i < PART_DIMENSION; ++i) {
      __m256i v_add_lo_0[4];
      __m256i v_add_hi_0[4];
      for (int dst = 0, src = 0; dst < 4; ++dst, src += 2) {
        v_add_lo_0[dst] = _mm256_add_epi32(v_madd_lo[src + 0][i], v_madd_lo[src + 1][i]);
        v_add_hi_0[dst] = _mm256_add_epi32(v_madd_hi[src + 0][i], v_madd_hi[src + 1][i]);
      }

      __m256i v_add_lo_1[2];
      __m256i v_add_hi_1[2];
      for (int dst = 0, src = 0; dst < 2; ++dst, src += 2) {
        v_add_lo_1[dst] = _mm256_add_epi32(v_add_lo_0[src + 0], v_add_lo_0[src + 1]);
        v_add_hi_1[dst] = _mm256_add_epi32(v_add_hi_0[src + 0], v_add_hi_0[src + 1]);
      }

      v_trunc_lo[i] = truncate_avx2(_mm256_add_epi32(v_add_lo_1[0], v_add_lo_1[1]), debias, shift_2nd);
      v_trunc_hi[i] = truncate_avx2(_mm256_add_epi32(v_add_hi_1[0], v_add_hi_1[1]), debias, shift_2nd);
    }
    __m256i v_result[PART_DIMENSION];
    for (int i = 0; i < PART_DIMENSION; ++i) {
      v_result[i] = _mm256_packs_epi32(v_trunc_lo[i], v_trunc_hi[i]);
    }

    for (int i = 0; i < PART_DIMENSION; ++i) {
      _mm256_store_si256((__m256i*)dst, v_result[i]);
      dst += 16;
    }
  }
  
#undef NUM_PARTS
#undef PART_DIMENSION

}


static void fast_inverse_tr_16x16_avx2_hor(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  //const int32_t* c_ptr = (const int32_t*)coeff; // Handle as 32 bit integer to load two coeffs into vector at the same time
  //const __m256i* v_src_raw = (const __m256i*)src;

  //__m256i v_madd_lo[8][16];
  //__m256i v_madd_hi[8][16];
  //for (int s = 0; s < 8; ++s) {
  //  __m256i v_src_lo = _mm256_unpacklo_epi16(v_src_raw[0], v_src_raw[1]);
  //  __m256i v_src_hi = _mm256_unpackhi_epi16(v_src_raw[0], v_src_raw[1]);
  //  v_src_raw += 2;

  //  for (int c = 0; c < 16; ++c) {
  //    const __m256i v_coeff = _mm256_set1_epi32(*c_ptr);
  //    v_madd_lo[s][c] = _mm256_madd_epi16(v_src_lo, v_coeff);
  //    v_madd_hi[s][c] = _mm256_madd_epi16(v_src_hi, v_coeff);
  //    c_ptr++;
  //  }
  //}

  //__m256i v_add_lo_0[4][16];
  //__m256i v_add_hi_0[4][16];
  //for (int s = 0, d = 0; d < 4; ++d, s += 2) {
  //  for (int c = 0; c < 16; ++c) {
  //    v_add_lo_0[d][c] = _mm256_add_epi32(v_madd_lo[s + 0][c], v_madd_lo[s + 1][c]);
  //    v_add_hi_0[d][c] = _mm256_add_epi32(v_madd_hi[s + 0][c], v_madd_hi[s + 1][c]);
  //  }
  //}

  //__m256i v_add_lo_1[2][16];
  //__m256i v_add_hi_1[2][16];
  //for (int s = 0, d = 0; d < 2; ++d, s += 2) {
  //  for (int c = 0; c < 16; ++c) {
  //    v_add_lo_1[d][c] = _mm256_add_epi32(v_add_lo_0[s + 0][c], v_add_lo_0[s + 1][c]);
  //    v_add_hi_1[d][c] = _mm256_add_epi32(v_add_hi_0[s + 0][c], v_add_hi_0[s + 1][c]);
  //  }
  //}

  //__m256i v_trunc_lo[16];
  //__m256i v_trunc_hi[16];
  //for (int c = 0; c < 16; ++c) {
  //  v_trunc_lo[c] = truncate_avx2(_mm256_add_epi32(v_add_lo_1[0][c], v_add_lo_1[1][c]), debias, shift);
  //  v_trunc_hi[c] = truncate_avx2(_mm256_add_epi32(v_add_hi_1[0][c], v_add_hi_1[1][c]), debias, shift);
  //}

  //for (int i = 0; i < 16; ++i) {
  //  dst[i] = _mm256_packs_epi32(v_trunc_lo[i], v_trunc_hi[i]);
  //}

  for (int j = 0; j < line; ++j) {
    __m256i res_0 = _mm256_setzero_si256();
    __m256i res_1 = _mm256_setzero_si256();

    __m256i *coeff_start = (__m256i*)coeff;
    for (int i = 0; i < 8; ++i) {
      int16_t source[2];
      source[0] = src[j + i * 32];
      source[1] = src[j + i * 32 + 16];
      int32_t paired_source;
      memcpy(&paired_source, source, sizeof(int32_t));

      __m256i v_src = _mm256_set1_epi32(paired_source);

      __m256i v_coeff0 = _mm256_load_si256(coeff_start);
      coeff_start++;
      __m256i v_coeff1 = _mm256_load_si256(coeff_start);
      coeff_start++;

      __m256i v_madd0 = _mm256_madd_epi16(v_src, v_coeff0);
      __m256i v_madd1 = _mm256_madd_epi16(v_src, v_coeff1);

      res_0 = _mm256_add_epi32(res_0, v_madd0);
      res_1 = _mm256_add_epi32(res_1, v_madd1);
    }

    __m256i v_trunc0 = truncate_avx2(res_0, debias, shift);
    __m256i v_trunc1 = truncate_avx2(res_1, debias, shift);

    __m256i packed = _mm256_packs_epi32(v_trunc0, v_trunc1);
    packed = _mm256_permute4x64_epi64(packed, _MM_SHUFFLE(3, 1, 2, 0));
    dst[j] = packed;
  }
}

static void fast_inverse_tr_16x16_avx2_ver(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  __m256i v_result[16];
  int16_t *src_p = (int16_t*)src;
  for (int j = 0; j < 16; ++j) {
    __m256i res_0 = _mm256_setzero_si256();
    __m256i res_1 = _mm256_setzero_si256();
    __m256i* coeff_start = (__m256i*)coeff;
    for (int i = 0; i < 8; ++i) {
      int16_t source[2];
      source[0] = src_p[j + i * 32];
      source[1] = src_p[j + i * 32 + 16];
      int32_t paired_source;
      memcpy(&paired_source, source, sizeof(int32_t));

      __m256i v_src = _mm256_set1_epi32(paired_source);

      __m256i coeff_0 = _mm256_load_si256(coeff_start);
      coeff_start++;
      __m256i coeff_1 = _mm256_load_si256(coeff_start);
      coeff_start++;

      __m256i madd0 = _mm256_madd_epi16(v_src, coeff_0);
      __m256i madd1 = _mm256_madd_epi16(v_src, coeff_1);

      res_0 = _mm256_add_epi32(res_0, madd0);
      res_1 = _mm256_add_epi32(res_1, madd1);
    }

    __m256i v_trunc0 = truncate_avx2(res_0, debias, shift);
    __m256i v_trunc1 = truncate_avx2(res_1, debias, shift);

    __m256i packed = _mm256_packs_epi32(v_trunc0, v_trunc1);
    packed = _mm256_permute4x64_epi64(packed, _MM_SHUFFLE(3, 1, 2, 0));
    _mm256_store_si256((__m256i *)dst, packed);
    dst += 16;
  }
  //const int32_t* c_ptr = (const int32_t*)coeff; // Handle as 32 bit integer to load two coeffs into vector at the same time
  //const __m256i* v_src_raw = src;

  //// Do a 32-bit transpose to arrange result from previous pass
  //__m256i v_tmp32_lo[8];
  //__m256i v_tmp32_hi[8];
  //for (int d = 0, s = 0; d < 8; ++d, s += 2) {
  //  v_tmp32_lo[d] = _mm256_unpacklo_epi32(v_src_raw[s + 0], v_src_raw[s + 1]);
  //  v_tmp32_hi[d] = _mm256_unpackhi_epi32(v_src_raw[s + 0], v_src_raw[s + 1]);
  //}

  //__m256i v_tmp64_lo[8];
  //__m256i v_tmp64_hi[8];
  //for (int d = 0, s = 0; d < 4; ++d, s += 2) {
  //  v_tmp64_lo[0 + d] = _mm256_unpacklo_epi64(v_tmp32_lo[s + 0], v_tmp32_lo[s + 1]);
  //  v_tmp64_lo[4 + d] = _mm256_unpacklo_epi64(v_tmp32_hi[s + 0], v_tmp32_hi[s + 1]);

  //  v_tmp64_hi[0 + d] = _mm256_unpackhi_epi64(v_tmp32_lo[s + 0], v_tmp32_lo[s + 1]);
  //  v_tmp64_hi[4 + d] = _mm256_unpackhi_epi64(v_tmp32_hi[s + 0], v_tmp32_hi[s + 1]);
  //}
  //
  //__m256i v_src[16];
  //v_src[ 0] = _mm256_permute2x128_si256(v_tmp64_lo[0], v_tmp64_lo[1], 0x20);
  //v_src[ 1] = _mm256_permute2x128_si256(v_tmp64_hi[0], v_tmp64_hi[1], 0x20);
  //v_src[ 2] = _mm256_permute2x128_si256(v_tmp64_lo[4], v_tmp64_lo[5], 0x20);
  //v_src[ 3] = _mm256_permute2x128_si256(v_tmp64_hi[4], v_tmp64_hi[5], 0x20);
  //v_src[ 4] = _mm256_permute2x128_si256(v_tmp64_lo[0], v_tmp64_lo[1], 0x31);
  //v_src[ 5] = _mm256_permute2x128_si256(v_tmp64_hi[0], v_tmp64_hi[1], 0x31);
  //v_src[ 6] = _mm256_permute2x128_si256(v_tmp64_lo[4], v_tmp64_lo[5], 0x31);
  //v_src[ 7] = _mm256_permute2x128_si256(v_tmp64_hi[4], v_tmp64_hi[5], 0x31);

  //v_src[ 8] = _mm256_permute2x128_si256(v_tmp64_lo[2], v_tmp64_lo[3], 0x20);
  //v_src[ 9] = _mm256_permute2x128_si256(v_tmp64_hi[2], v_tmp64_hi[3], 0x20);
  //v_src[10] = _mm256_permute2x128_si256(v_tmp64_lo[6], v_tmp64_lo[7], 0x20);
  //v_src[11] = _mm256_permute2x128_si256(v_tmp64_hi[6], v_tmp64_hi[7], 0x20);
  //v_src[12] = _mm256_permute2x128_si256(v_tmp64_lo[2], v_tmp64_lo[3], 0x31);
  //v_src[13] = _mm256_permute2x128_si256(v_tmp64_hi[2], v_tmp64_hi[3], 0x31);
  //v_src[14] = _mm256_permute2x128_si256(v_tmp64_lo[6], v_tmp64_lo[7], 0x31);
  //v_src[15] = _mm256_permute2x128_si256(v_tmp64_hi[6], v_tmp64_hi[7], 0x31);

  //__m256i v_madd_0[8][16];
  //__m256i v_madd_1[8][16];
  //for (int s = 0; s < 8; ++s) {
  //  for (int c = 0; c < 16; ++c) {
  //    const __m256i v_coeff = _mm256_set1_epi32(*c_ptr);
  //    v_madd_0[s][c] = _mm256_madd_epi16(v_src[0 + s], v_coeff);
  //    v_madd_1[s][c] = _mm256_madd_epi16(v_src[8 + s], v_coeff);
  //    c_ptr++;
  //  }
  //}

  //__m256i v_add_00[4][16];
  //__m256i v_add_01[4][16];
  //for (int s = 0, d = 0; d < 4; ++d, s += 2) {
  //  for (int c = 0; c < 16; ++c) {
  //    v_add_00[d][c] = _mm256_add_epi32(v_madd_0[s + 0][c], v_madd_0[s + 1][c]);
  //    v_add_01[d][c] = _mm256_add_epi32(v_madd_1[s + 0][c], v_madd_1[s + 1][c]);
  //  }
  //}

  //__m256i v_add_10[2][16];
  //__m256i v_add_11[2][16];
  //for (int s = 0, d = 0; d < 2; ++d, s += 2) {
  //  for (int c = 0; c < 16; ++c) {
  //    v_add_10[d][c] = _mm256_add_epi32(v_add_00[s + 0][c], v_add_00[s + 1][c]);
  //    v_add_11[d][c] = _mm256_add_epi32(v_add_01[s + 0][c], v_add_01[s + 1][c]);
  //  }
  //}

  //__m256i v_trunc_0[16];
  //__m256i v_trunc_1[16];
  //for (int c = 0; c < 16; ++c) {
  //  v_trunc_0[c] = truncate_avx2(_mm256_add_epi32(v_add_10[0][c], v_add_10[1][c]), debias, shift);
  //  v_trunc_1[c] = truncate_avx2(_mm256_add_epi32(v_add_11[0][c], v_add_11[1][c]), debias, shift);
  //}

  //__m256i v_result[16];
  //for (int d = 0; d < 16; ++d) {
  //  v_result[d] = _mm256_packs_epi32(v_trunc_0[d], v_trunc_1[d]);
  //}
  //for (int d = 0; d < 16; ++d) {
  //  v_result[d] = _mm256_permute4x64_epi64(v_result[d], _MM_SHUFFLE(3, 1, 2, 0));
  //}
  
  //transpose_avx2(v_result, (__m256i*)dst, 16, 16);
}

static void fast_inverse_tr_16x16_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 16;
  const int height = 16;

  int skip_width = 0;
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* hor_coeff = fi_dct2_16x16_coeff_hor;
  const int16_t* ver_coeff = fi_dct2_16x16_coeff_ver;
  if (hor == DST7) {
    hor_coeff = fi_dst7_16x16_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = fi_dct8_16x16_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = fi_dst7_16x16_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = fi_dct8_16x16_coeff_ver;
  }

  __m256i v_hor_pass_out[16];
  fast_inverse_tr_16x16_avx2_hor(src, v_hor_pass_out, ver_coeff, shift_1st, height, 0, skip_width);

  fast_inverse_tr_16x16_avx2_ver(v_hor_pass_out, dst, hor_coeff, shift_2nd, width, skip_width, skip_height);
}


static void fast_forward_tr_16x32_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 16;
  const int height = 32;
  
  int skip_width = 0;
  int skip_height = (ver != DCT2 && height == 32) ? 16 : (height > 32 ? height - 32 : 0);

  const int  reduced_line = width - skip_width;
  const int  cutoff = height - skip_height;
  int16_t* p_dst = dst;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_16xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_16x32_coeff_ver;
  if (hor == DST7) {
    hor_coeff = ff_dst7_16xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = ff_dct8_16xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = ff_dst7_16x32_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = ff_dct8_16x32_coeff_ver;
  }

  ALIGNED(32) int16_t v_hor_pass_out[32*16];
  fast_forward_DCT2_B16_avx2_hor(src, (__m256i*)v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);


  __m256i temp_out[32];
  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);
  if(ver == DCT2) {
    for (int j = 0; j < 16; ++j) {
      __m256i res_0 = _mm256_setzero_si256();
      __m256i res_1 = _mm256_setzero_si256();
      __m256i res_2 = _mm256_setzero_si256();
      __m256i res_3 = _mm256_setzero_si256();
      const int16_t* coeff_start = ver_coeff;
      for (int i = 0; i < 16; ++i) {
        int16_t source[2];
        source[0] = v_hor_pass_out[j + i * 32];
        source[1] = v_hor_pass_out[j + i * 32 + 16];
        int32_t paired_source;
        memcpy(&paired_source, source, sizeof(int32_t));

        __m256i v_src = _mm256_set1_epi32(paired_source);
        __m256i v_coeff_0 = _mm256_load_si256((__m256i*) coeff_start);
        coeff_start += 16;
        __m256i v_coeff_1 = _mm256_load_si256((__m256i*) coeff_start);
        coeff_start += 16;
        __m256i v_coeff_2 = _mm256_load_si256((__m256i*) coeff_start);
        coeff_start += 16;
        __m256i v_coeff_3 = _mm256_load_si256((__m256i*) coeff_start);
        coeff_start += 16;

        __m256i madd_0 = _mm256_madd_epi16(v_src, v_coeff_0);
        __m256i madd_1 = _mm256_madd_epi16(v_src, v_coeff_1);
        __m256i madd_2 = _mm256_madd_epi16(v_src, v_coeff_2);
        __m256i madd_3 = _mm256_madd_epi16(v_src, v_coeff_3);

        res_0 = _mm256_add_epi32(res_0, madd_0);
        res_1 = _mm256_add_epi32(res_1, madd_1);
        res_2 = _mm256_add_epi32(res_2, madd_2);
        res_3 = _mm256_add_epi32(res_3, madd_3);
      }
      __m256i v_trunc_0 = truncate_avx2(res_0, debias, shift_2nd);
      __m256i v_trunc_1 = truncate_avx2(res_1, debias, shift_2nd);
      __m256i v_trunc_2 = truncate_avx2(res_2, debias, shift_2nd);
      __m256i v_trunc_3 = truncate_avx2(res_3, debias, shift_2nd);

      v_trunc_0 = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
      v_trunc_1 = _mm256_packs_epi32(v_trunc_2, v_trunc_3);
      v_trunc_0 = _mm256_permute4x64_epi64(v_trunc_0, _MM_SHUFFLE(3, 1, 2, 0));
      v_trunc_1 = _mm256_permute4x64_epi64(v_trunc_1, _MM_SHUFFLE(3, 1, 2, 0));
      _mm256_store_si256(temp_out + 2 * j, v_trunc_0);
      _mm256_store_si256(temp_out + 2 * j + 1, v_trunc_1);
    }
    transpose_avx2(temp_out, (__m256i*) dst, 32, 16);
  }
  else {
    for (int j = 0; j < 16; ++j) {
    __m256i res_0 = _mm256_setzero_si256();
    __m256i res_1 = _mm256_setzero_si256();
    const int16_t* coeff_start = ver_coeff;
    for (int i = 0; i < 16; ++i) {
      int16_t source[2];
      source[0] = v_hor_pass_out[j + i * 32];
      source[1] = v_hor_pass_out[j + i * 32 + 16];
      int32_t paired_source;
      memcpy(&paired_source, source, sizeof(int32_t));

      __m256i v_src = _mm256_set1_epi32(paired_source);
      __m256i v_coeff_0 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;
      __m256i v_coeff_1 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 48;

      __m256i madd_0 = _mm256_madd_epi16(v_src, v_coeff_0);
      __m256i madd_1 = _mm256_madd_epi16(v_src, v_coeff_1);

      res_0 = _mm256_add_epi32(res_0, madd_0);
      res_1 = _mm256_add_epi32(res_1, madd_1);
    }
    __m256i v_trunc_0 = truncate_avx2(res_0, debias, shift_2nd);
    __m256i v_trunc_1 = truncate_avx2(res_1, debias, shift_2nd);

    v_trunc_0 = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
    v_trunc_0 = _mm256_permute4x64_epi64(v_trunc_0, _MM_SHUFFLE(3, 1, 2, 0));
    _mm256_store_si256(temp_out + 2 * j, v_trunc_0);
  }
  transpose_avx2(temp_out, (__m256i*) dst, 32, 16);
  }
#if 0
  // To how many parts the vertical pass should be split.
  // At least on my testing it seems that there is no further gain by splitting to more than 4 parts.
#define NUM_PARTS 4
#define PART_DIMENSION (32/NUM_PARTS)
  for (int part = 0; part < NUM_PARTS; ++part) {
    // Got 32 / NUM_PARTS lines of samples. Handle two lines at a time (beacuse of unpack)
    __m256i v_madd_hi[16][PART_DIMENSION];
    __m256i v_madd_lo[16][PART_DIMENSION];
    // Samples are the same between the parts
    __m256i* v_src_ptr = v_hor_pass_out;
    // However for coefficients, the starting point needs to be adjusted
    const int32_t* line_coeff = (const int32_t*)ver_coeff + PART_DIMENSION * part;
    for (int i = 0; i < 16; ++i) {
      __m256i v_src_hi = _mm256_unpackhi_epi16(v_src_ptr[0], v_src_ptr[1]);
      __m256i v_src_lo = _mm256_unpacklo_epi16(v_src_ptr[0], v_src_ptr[1]);

      // Apply coefficients
      // TODO: Here try loading the coefficient directly instead of set1
      for (int ii = 0; ii < PART_DIMENSION; ++ii) {
        const int32_t coeff = line_coeff[ii];
        const __m256i v_coeff = _mm256_set1_epi32(coeff);
        v_madd_hi[i][ii] = _mm256_madd_epi16(v_src_hi, v_coeff);
        v_madd_lo[i][ii] = _mm256_madd_epi16(v_src_lo, v_coeff);
      }

      line_coeff += 32;
      v_src_ptr += 2;
    }

    for (int ii = 0; ii < PART_DIMENSION; ++ii) {
      // First round of additions
      __m256i v_add_hi_0[8];
      __m256i v_add_lo_0[8];
      for (int i = 0; i < 8; ++i) {
        const int offset = i * 2;
        v_add_hi_0[i] = _mm256_add_epi32(v_madd_hi[offset][ii], v_madd_hi[offset + 1][ii]);
        v_add_lo_0[i] = _mm256_add_epi32(v_madd_lo[offset][ii], v_madd_lo[offset + 1][ii]);
      }

      // Second round of additions
      __m256i v_add_hi_1[4];
      __m256i v_add_lo_1[4];
      for (int i = 0; i < 4; ++i) {
        const int offset = i * 2;
        v_add_hi_1[i] = _mm256_add_epi32(v_add_hi_0[offset], v_add_hi_0[offset + 1]);
        v_add_lo_1[i] = _mm256_add_epi32(v_add_lo_0[offset], v_add_lo_0[offset + 1]);
      }

      // Third round of addtions
      __m256i v_add_hi_2[2];
      __m256i v_add_lo_2[2];
      for (int i = 0; i < 2; ++i) {
        const int offset = i * 2;
        v_add_hi_2[i] = _mm256_add_epi32(v_add_hi_1[offset], v_add_hi_1[offset + 1]);
        v_add_lo_2[i] = _mm256_add_epi32(v_add_lo_1[offset], v_add_lo_1[offset + 1]);
      }

      // Final round of additions, truncate and store
      __m256i v_trunc_hi = truncate_avx2(_mm256_add_epi32(v_add_hi_2[0], v_add_hi_2[1]), debias, shift_2nd);
      __m256i v_trunc_lo = truncate_avx2(_mm256_add_epi32(v_add_lo_2[0], v_add_lo_2[1]), debias, shift_2nd);
      __m256i v_result = _mm256_packs_epi32(v_trunc_lo, v_trunc_hi);
      _mm256_store_si256((__m256i*)dst, v_result);

      dst += 16;
    }
  }
#undef NUM_PARTS
#undef PART_DIMENSION
#endif

}


static void fast_inverse_tr_16x32_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const int32_t* c_ptr = (int32_t*)coeff; // Handle as 32 bit integer to load two coeffs into vectors at a time
  const __m256i* v_src_raw = (const __m256i*)src;

  __m256i v_tmp16_lo[16];
  __m256i v_tmp16_hi[16];
  for (int d = 0, s = 0; d < 16; ++d, s += 2) {
    v_tmp16_lo[d] = _mm256_unpacklo_epi16(v_src_raw[s + 0], v_src_raw[s + 1]);
    v_tmp16_hi[d] = _mm256_unpackhi_epi16(v_src_raw[s + 0], v_src_raw[s + 1]);
  }
  int  row = 0;
  for (; row < 32 - skip_line2; ++row) {
    __m256i v_res_lo = _mm256_setzero_si256();
    __m256i v_res_hi = _mm256_setzero_si256();
    for (int i = 0; i < 16; ++i) {
      const __m256i v_coeff = _mm256_set1_epi32(*c_ptr);
      __m256i v_madd_lo = _mm256_madd_epi16(v_tmp16_lo[i], v_coeff);
      __m256i v_madd_hi = _mm256_madd_epi16(v_tmp16_hi[i], v_coeff);
      c_ptr++;

      v_res_lo = _mm256_add_epi32(v_res_lo, v_madd_lo);
      v_res_hi = _mm256_add_epi32(v_res_hi, v_madd_hi);
    }

    __m256i v_trunc_lo = truncate_avx2(v_res_lo, debias, shift);
    __m256i v_trunc_hi = truncate_avx2(v_res_hi, debias, shift);

    __m256i packed = _mm256_packs_epi32(v_trunc_lo, v_trunc_hi);
    dst[row] = packed;
  }

  for (; row < 32; ++row) {
    dst[row] = _mm256_setzero_si256();
  }
}

static void fast_inverse_tr_16x32_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  int32_t * src_32 = (int32_t *)src;
  for (int row = 0, d = 0; row < 32; ++row) {
    __m256i v_res_0 = _mm256_setzero_si256();
    __m256i v_res_1 = _mm256_setzero_si256();
    __m256i *coeff_start = (__m256i*) coeff;
    for (int i = 0; i < 8; ++i) {
      __m256i v_src = _mm256_set1_epi32(*src_32);
      src_32++;

      __m256i v_madd_0 = _mm256_madd_epi16(v_src, _mm256_load_si256(coeff_start));
      coeff_start++;
      __m256i v_madd_1 = _mm256_madd_epi16(v_src, _mm256_load_si256(coeff_start));
      coeff_start++;

      v_res_0 = _mm256_add_epi32(v_res_0, v_madd_0);
      v_res_1 = _mm256_add_epi32(v_res_1, v_madd_1);
    }

    __m256i v_trunc_0 = truncate_avx2(v_res_0, debias, shift);
    __m256i v_trunc_1 = truncate_avx2(v_res_1, debias, shift);

    __m256i packed = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
    packed = _mm256_permute4x64_epi64(packed, _MM_SHUFFLE(3, 1, 2, 0));
    _mm256_store_si256((__m256i*) dst + row, packed);
  }
}

static void fast_inverse_tr_16x32_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 16;
  const int height = 32;

  int skip_width = 0;
  int skip_height = (ver != DCT2 && height == 32) ? 16 : (height > 32 ? height - 32 : 0);

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = &uvg_g_dct_32_t[0][0];
  const int16_t* hor_coeff = fi_dct2_16x16_coeff_hor;
  if (hor == DST7) {
    hor_coeff = fi_dst7_16x32_coeff_hor; // TODO: coeffs
  } else if (hor == DCT8) {
    hor_coeff = fi_dct8_16x32_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = &uvg_g_dst7_32_t[0][0];
  } else if (ver == DCT8) {
    ver_coeff = &uvg_g_dct8_32[0][0];
  }

  __m256i v_ver_pass_out[32];
  fast_inverse_tr_16x32_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, 0);
  int16_t* ver_pass_out = (int16_t*)v_ver_pass_out;
  fast_inverse_tr_16x32_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static void fast_forward_DCT2_B32_avx2_hor(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2) {

  const int32_t    add = (shift > 0) ? (1 << (shift - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  const int reduced_line = line - skip_line;

  for(int j = 0; j < reduced_line; ++j) {
    int32_t source[16];
    memcpy(source, src, sizeof(int16_t) * 32);
    src += 32;

    __m256i res_0 = _mm256_setzero_si256();
    __m256i res_1 = _mm256_setzero_si256();
    __m256i res_2 = _mm256_setzero_si256();
    __m256i res_3 = _mm256_setzero_si256();
    const int16_t *coeff_start = coeff;
    for(int i = 0; i < 16; i++) {
      __m256i v_src = _mm256_set1_epi32(source[i]);
      __m256i v_coeff_0 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;
      __m256i v_coeff_1 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;
      __m256i v_coeff_2 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;
      __m256i v_coeff_3 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;

      __m256i madd_0 = _mm256_madd_epi16(v_src, v_coeff_0);
      __m256i madd_1 = _mm256_madd_epi16(v_src, v_coeff_1);
      __m256i madd_2 = _mm256_madd_epi16(v_src, v_coeff_2);
      __m256i madd_3 = _mm256_madd_epi16(v_src, v_coeff_3);

      res_0 = _mm256_add_epi32(res_0, madd_0);
      res_1 = _mm256_add_epi32(res_1, madd_1);
      res_2 = _mm256_add_epi32(res_2, madd_2);
      res_3 = _mm256_add_epi32(res_3, madd_3);
    }
    __m256i v_trunc_0 = truncate_avx2(res_0, debias, shift);
    __m256i v_trunc_1 = truncate_avx2(res_1, debias, shift);
    __m256i v_trunc_2 = truncate_avx2(res_2, debias, shift);
    __m256i v_trunc_3 = truncate_avx2(res_3, debias, shift);

    v_trunc_0 = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
    v_trunc_1 = _mm256_packs_epi32(v_trunc_2, v_trunc_3);

    if(line == 32) {
      v_trunc_0 = _mm256_permute4x64_epi64(v_trunc_0, _MM_SHUFFLE(3, 1, 2, 0));
      v_trunc_1 = _mm256_permute4x64_epi64(v_trunc_1, _MM_SHUFFLE(3, 1, 2, 0));
    }

    _mm256_store_si256(dst, v_trunc_0);
    dst++;
    _mm256_store_si256(dst, v_trunc_1);
    dst++;
  }
}

static void fast_forward_DCT8_B32_avx2_hor(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2) {
  const int32_t    add = (shift > 0) ? (1 << (shift - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  const int  cutoff = 32 - skip_line2;
  const int reduced_line = line - skip_line;

  ALIGNED(32) int16_t temp_source[32 * 32];
  __m256i* v_src_p = (__m256i*) src;
  for (int i = 0; i < reduced_line / 2; ++i) {
    __m256i first_half_lo = _mm256_unpacklo_epi32(v_src_p[i * 4], v_src_p[i * 4 + 2]);
    __m256i first_half_hi = _mm256_unpackhi_epi32(v_src_p[i * 4], v_src_p[i * 4 + 2]);
    __m256i second_half_lo = _mm256_unpacklo_epi32(v_src_p[i * 4 + 1], v_src_p[i * 4 + 3]);
    __m256i second_half_hi = _mm256_unpackhi_epi32(v_src_p[i * 4 + 1], v_src_p[i * 4 + 3]);

    _mm256_store_si256((__m256i*)temp_source + i * 4, first_half_lo);
    _mm256_store_si256((__m256i*)temp_source + i * 4 + 1, first_half_hi);
    _mm256_store_si256((__m256i*)temp_source + i * 4 + 2, second_half_lo);
    _mm256_store_si256((__m256i*)temp_source + i * 4 + 3, second_half_hi);
  }

  for (int j = 0; j < reduced_line / 2; ++j) {

    int32_t source[32];
    memcpy(source, temp_source + 64 * j, sizeof(int16_t) * 64);

    __m256i res_0 = _mm256_setzero_si256();
    __m256i res_1 = _mm256_setzero_si256();
    __m256i res_2 = _mm256_setzero_si256();
    __m256i res_3 = _mm256_setzero_si256();
    const int16_t* coeff_start = coeff;

    for (int i = 0; i < 32; i += 2) {
      __m256i v_src0 = _mm256_set1_epi32(source[i]);      
      __m256i v_src1 = _mm256_set1_epi32(source[i + 1]);

      __m256i v_coeff_0 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;
      __m256i v_coeff_1 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 48;

      __m256i madd_0 = _mm256_madd_epi16(v_src0, v_coeff_0);
      __m256i madd_1 = _mm256_madd_epi16(v_src0, v_coeff_1);
      __m256i madd_2 = _mm256_madd_epi16(v_src1, v_coeff_0);
      __m256i madd_3 = _mm256_madd_epi16(v_src1, v_coeff_1);

      res_0 = _mm256_add_epi32(madd_0, res_0);
      res_1 = _mm256_add_epi32(madd_1, res_1);
      res_2 = _mm256_add_epi32(madd_2, res_2);
      res_3 = _mm256_add_epi32(madd_3, res_3);
    }
    __m256i v_trunc_0 = truncate_avx2(res_0, debias, shift);
    __m256i v_trunc_1 = truncate_avx2(res_1, debias, shift);
    __m256i v_trunc_2 = truncate_avx2(res_2, debias, shift);
    __m256i v_trunc_3 = truncate_avx2(res_3, debias, shift);
    
    v_trunc_0 = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
    v_trunc_2 = _mm256_packs_epi32(v_trunc_2, v_trunc_3);

    if (line == 32) {
      v_trunc_0 = _mm256_permute4x64_epi64(v_trunc_0, _MM_SHUFFLE(3, 1, 2, 0));
      v_trunc_2 = _mm256_permute4x64_epi64(v_trunc_2, _MM_SHUFFLE(3, 1, 2, 0));
    }
    _mm256_store_si256(dst, v_trunc_0);
    dst+=2;
    _mm256_store_si256(dst, v_trunc_2);
    dst+=2;
  }
}


static void fast_forward_DCT2_32x2_avx2_ver(const __m256i* src, int16_t* dst, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = (shift > 0) ? (1 << (shift - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_src_ptr = src;

  // Prepare coeffs
  // TODO: either rename these old coeff tables to be consistent with other new avx2 functions
  // or construct them here in place. Should be ease to accomplish with set1_epi32, just use a int32_t combined from two int16_t
  const __m256i v_coeff_0 = _mm256_load_si256((const __m256i*)&fast_forward_dct2_b2_coeff[0]);
  const __m256i v_coeff_1 = _mm256_load_si256((const __m256i*)&fast_forward_dct2_b2_coeff[16]);
  
  // Got data for 4 vectors, 32 lines with 2 samples each
  __m256i v_result_e[4];
  __m256i v_result_o[4];
  for (int j = 0; j < 4; ++j) {
    const __m256i v_src = v_src_ptr[0];

    v_result_e[j] = truncate_avx2(_mm256_madd_epi16(v_src, v_coeff_0), debias, shift);
    v_result_o[j] = truncate_avx2(_mm256_madd_epi16(v_src, v_coeff_1), debias, shift);

    v_src_ptr++;
  }

  __m256i v_tmp[4];
  v_tmp[0] = _mm256_packs_epi32(v_result_e[0], v_result_e[1]);
  v_tmp[1] = _mm256_packs_epi32(v_result_e[2], v_result_e[3]);
  v_tmp[2] = _mm256_packs_epi32(v_result_o[0], v_result_o[1]);
  v_tmp[3] = _mm256_packs_epi32(v_result_o[2], v_result_o[3]);

  v_tmp[0] = _mm256_permute4x64_epi64(v_tmp[0], _MM_SHUFFLE(3, 1, 2, 0));
  v_tmp[1] = _mm256_permute4x64_epi64(v_tmp[1], _MM_SHUFFLE(3, 1, 2, 0));
  v_tmp[2] = _mm256_permute4x64_epi64(v_tmp[2], _MM_SHUFFLE(3, 1, 2, 0));
  v_tmp[3] = _mm256_permute4x64_epi64(v_tmp[3], _MM_SHUFFLE(3, 1, 2, 0));

  for (int i = 0; i < 4; ++i) {
    _mm256_store_si256((__m256i*)&dst[i * 16], v_tmp[i]);
  }
}

static void fast_forward_DCT2_32x4_avx2_ver(const __m256i* src, int16_t* dst, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = (shift > 0) ? (1 << (shift - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  // Got data for 8 vectors, 32 lines with 4 samples each

  // Prepare coeffs
  const int16_t* coeff = &uvg_g_dct_4[0][0];
  const int a = coeff[0];
  const int b = coeff[1 * 4 + 0];
  const int c = coeff[1 * 4 + 1];

  __m256i v_coeff_0 = _mm256_set1_epi16(a);
  __m256i v_coeff_1 = _mm256_setr_epi16(b, c, -c, -b, b, c, -c, -b, b, c, -c, -b, b, c, -c, -b);
  __m256i v_coeff_2 = _mm256_setr_epi16(a, -a, -a, a, a, -a, -a, a, a, -a, -a, a, a, -a, -a, a);
  __m256i v_coeff_3 = _mm256_setr_epi16(c, -b, b, -c, c, -b, b, -c, c, -b, b, -c, c, -b, b, -c);

  const __m256i* v_src_ptr = src;
  __m256i v_trunc_0[8];
  __m256i v_trunc_1[8];
  for (int j = 0; j < 8; ++j) {
    __m256i v_madd_0 = _mm256_madd_epi16(v_src_ptr[0], v_coeff_0);
    __m256i v_madd_1 = _mm256_madd_epi16(v_src_ptr[0], v_coeff_1);
    __m256i v_madd_2 = _mm256_madd_epi16(v_src_ptr[0], v_coeff_2);
    __m256i v_madd_3 = _mm256_madd_epi16(v_src_ptr[0], v_coeff_3);

    v_trunc_0[j] = truncate_avx2(_mm256_hadd_epi32(v_madd_0, v_madd_1), debias, shift);
    v_trunc_1[j] = truncate_avx2(_mm256_hadd_epi32(v_madd_2, v_madd_3), debias, shift);

    v_src_ptr++;
  }

  __m256i v_result[8];
  __m256i v_tmp[8];
  for (int i = 0; i < 8; ++i) {
    v_trunc_0[i] = _mm256_permute4x64_epi64(v_trunc_0[i], _MM_SHUFFLE(3, 1, 2, 0));
    v_trunc_1[i] = _mm256_permute4x64_epi64(v_trunc_1[i], _MM_SHUFFLE(3, 1, 2, 0));
  }
  v_tmp[0] = _mm256_packs_epi32(v_trunc_0[0], v_trunc_0[1]);
  v_tmp[1] = _mm256_packs_epi32(v_trunc_0[2], v_trunc_0[3]);
  v_tmp[2] = _mm256_packs_epi32(v_trunc_0[4], v_trunc_0[5]);
  v_tmp[3] = _mm256_packs_epi32(v_trunc_0[6], v_trunc_0[7]);
  v_tmp[4] = _mm256_packs_epi32(v_trunc_1[0], v_trunc_1[1]);
  v_tmp[5] = _mm256_packs_epi32(v_trunc_1[2], v_trunc_1[3]);
  v_tmp[6] = _mm256_packs_epi32(v_trunc_1[4], v_trunc_1[5]);
  v_tmp[7] = _mm256_packs_epi32(v_trunc_1[6], v_trunc_1[7]);

  v_result[0] = _mm256_permute2x128_si256(v_tmp[0], v_tmp[1], 0x20);
  v_result[1] = _mm256_permute2x128_si256(v_tmp[2], v_tmp[3], 0x20);
  v_result[2] = _mm256_permute2x128_si256(v_tmp[0], v_tmp[1], 0x31);
  v_result[3] = _mm256_permute2x128_si256(v_tmp[2], v_tmp[3], 0x31);

  v_result[4] = _mm256_permute2x128_si256(v_tmp[4], v_tmp[5], 0x20);
  v_result[5] = _mm256_permute2x128_si256(v_tmp[6], v_tmp[7], 0x20);
  v_result[6] = _mm256_permute2x128_si256(v_tmp[4], v_tmp[5], 0x31);
  v_result[7] = _mm256_permute2x128_si256(v_tmp[6], v_tmp[7], 0x31);

  for (int i = 0; i < 8; ++i) {
    _mm256_store_si256((__m256i*)&dst[i * 16], v_result[i]);
  }
}


static void fast_forward_DCT2_32x8_avx2_ver(const __m256i* src, int16_t* dst, int32_t shift, int line, int skip_line, int skip_line2)
{
  int16_t* const p_dst = dst;
  const int32_t    add = (shift > 0) ? (1 << (shift - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  // Re-use coeff table
  const __m256i* v_coeff = (const __m256i*)ff_dct2_16x8_coeff_ver;

  const int reduced_line = line - skip_line;
  const __m256i* v_src_ptr = src;
  __m256i v_tmp_result[16];
  // Handle 2 lines at a time (16 samples, 8 samples per line)
  for (int i = 0; i < 16; ++i) {
    //                    line 1                    line 2
    // src vector:       [s0 s1 s2 s3 s4 s5 s6 s7 | s0 s1 s2 s3 s4 s5 s6 s7]
    // __m256i    v_src = _mm256_load_si256((const __m256i*)src);

    // Rearrange source in a way samples can be added together column-wise using add
    // after first round of madd operations.
    // Need 4 source vectors arranged as follows. High 128 lanes are the same as low:
    // vec_01 = [s0 s1 s0 s1 s0 s1 s0 s1 |...]
    // vec_02 = [s2 s3 s2 s3 s2 s3 s2 s3 |...]
    // vec_03 = [s4 s5 s4 s5 s4 s5 s4 s5 |...]
    // vec_04 = [s6 s7 s6 s7 s6 s7 s6 s7 |...]

    __m256i  v_src_0 = _mm256_shuffle_epi32(v_src_ptr[0], _MM_SHUFFLE(0, 0, 0, 0));
    __m256i  v_src_1 = _mm256_shuffle_epi32(v_src_ptr[0], _MM_SHUFFLE(1, 1, 1, 1));
    __m256i  v_src_2 = _mm256_shuffle_epi32(v_src_ptr[0], _MM_SHUFFLE(2, 2, 2, 2));
    __m256i  v_src_3 = _mm256_shuffle_epi32(v_src_ptr[0], _MM_SHUFFLE(3, 3, 3, 3));

    // Lane 1
    __m256i v_madd_0 = _mm256_madd_epi16(v_src_0, v_coeff[0]);
    __m256i v_madd_1 = _mm256_madd_epi16(v_src_1, v_coeff[1]);
    __m256i v_madd_2 = _mm256_madd_epi16(v_src_2, v_coeff[2]);
    __m256i v_madd_3 = _mm256_madd_epi16(v_src_3, v_coeff[3]);

    __m256i v_add_00 = _mm256_add_epi32(v_madd_0, v_madd_1);
    __m256i v_add_01 = _mm256_add_epi32(v_madd_2, v_madd_3);

    __m256i v_add_10 = _mm256_add_epi32(v_add_00, v_add_01);

    // Lane 2
    __m256i v_madd_4 = _mm256_madd_epi16(v_src_0, v_coeff[4]);
    __m256i v_madd_5 = _mm256_madd_epi16(v_src_1, v_coeff[5]);
    __m256i v_madd_6 = _mm256_madd_epi16(v_src_2, v_coeff[6]);
    __m256i v_madd_7 = _mm256_madd_epi16(v_src_3, v_coeff[7]);

    __m256i v_add_02 = _mm256_add_epi32(v_madd_4, v_madd_5);
    __m256i v_add_03 = _mm256_add_epi32(v_madd_6, v_madd_7);

    __m256i v_add_11 = _mm256_add_epi32(v_add_02, v_add_03);

    // Trunc results from both lanes
    __m256i v_trunc_0 = truncate_avx2(v_add_10, debias, shift);
    __m256i v_trunc_1 = truncate_avx2(v_add_11, debias, shift);

    v_tmp_result[i] = _mm256_packs_epi32(v_trunc_0, v_trunc_1);

    v_src_ptr++;
  }

  __m256i v_result[16];
  transpose_avx2(v_tmp_result, v_result, 8, 32);

  for (int i = 0; i < 16; ++i) {
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }

}


static void fast_forward_tr_32x2_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 32;
  const int height = 2;
  
  int skip_width  = (hor != DCT2 && width == 32) ? 16 : width > 32 ? width - 32 : 0;
  int skip_height = 0;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_32xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_32x2_coeff_ver;

  __m256i v_hor_pass_out[4];
  fast_forward_DCT2_B32_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)ver_coeff;

  // Got samples for 4 source vectors, 2 lines 32 samples each
  __m256i v_src_hi_0 = _mm256_unpackhi_epi16(v_hor_pass_out[0], v_hor_pass_out[2]);
  __m256i v_src_lo_0 = _mm256_unpacklo_epi16(v_hor_pass_out[0], v_hor_pass_out[2]);
  __m256i v_src_hi_1 = _mm256_unpackhi_epi16(v_hor_pass_out[1], v_hor_pass_out[3]);
  __m256i v_src_lo_1 = _mm256_unpacklo_epi16(v_hor_pass_out[1], v_hor_pass_out[3]);
  
  __m256i v_madd_hi_00 = _mm256_madd_epi16(v_src_hi_0, v_coeff[0]);
  __m256i v_madd_hi_01 = _mm256_madd_epi16(v_src_hi_0, v_coeff[1]);
  __m256i v_madd_hi_10 = _mm256_madd_epi16(v_src_hi_1, v_coeff[0]);
  __m256i v_madd_hi_11 = _mm256_madd_epi16(v_src_hi_1, v_coeff[1]);

  __m256i v_madd_lo_00 = _mm256_madd_epi16(v_src_lo_0, v_coeff[0]);
  __m256i v_madd_lo_01 = _mm256_madd_epi16(v_src_lo_0, v_coeff[1]);
  __m256i v_madd_lo_10 = _mm256_madd_epi16(v_src_lo_1, v_coeff[0]);
  __m256i v_madd_lo_11 = _mm256_madd_epi16(v_src_lo_1, v_coeff[1]);

  __m256i v_trunc_hi_00 = truncate_avx2(v_madd_hi_00, debias, shift_2nd);
  __m256i v_trunc_hi_01 = truncate_avx2(v_madd_hi_01, debias, shift_2nd);
  __m256i v_trunc_hi_10 = truncate_avx2(v_madd_hi_10, debias, shift_2nd);
  __m256i v_trunc_hi_11 = truncate_avx2(v_madd_hi_11, debias, shift_2nd);

  __m256i v_trunc_lo_00 = truncate_avx2(v_madd_lo_00, debias, shift_2nd);
  __m256i v_trunc_lo_01 = truncate_avx2(v_madd_lo_01, debias, shift_2nd);
  __m256i v_trunc_lo_10 = truncate_avx2(v_madd_lo_10, debias, shift_2nd);
  __m256i v_trunc_lo_11 = truncate_avx2(v_madd_lo_11, debias, shift_2nd);

  __m256i v_result_0 = _mm256_packs_epi32(v_trunc_lo_00, v_trunc_hi_00);
  __m256i v_result_1 = _mm256_packs_epi32(v_trunc_lo_10, v_trunc_hi_10); // Swap middle hi-lo lanes
  __m256i v_result_2 = _mm256_packs_epi32(v_trunc_lo_01, v_trunc_hi_01);
  __m256i v_result_3 = _mm256_packs_epi32(v_trunc_lo_11, v_trunc_hi_11);

  // Swap middle 64-bit chunks
  v_result_0 = _mm256_permute4x64_epi64(v_result_0, _MM_SHUFFLE(3, 1, 2, 0));
  v_result_1 = _mm256_permute4x64_epi64(v_result_1, _MM_SHUFFLE(3, 1, 2, 0));
  v_result_2 = _mm256_permute4x64_epi64(v_result_2, _MM_SHUFFLE(3, 1, 2, 0));
  v_result_3 = _mm256_permute4x64_epi64(v_result_3, _MM_SHUFFLE(3, 1, 2, 0));

  _mm256_store_si256((__m256i*)dst, v_result_0);
  _mm256_store_si256((__m256i*)(dst + 16), v_result_1);
  _mm256_store_si256((__m256i*)(dst + 32), v_result_2);
  _mm256_store_si256((__m256i*)(dst + 48), v_result_3);
}


static void fast_inverse_tr_32x2_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i v_coeff_0 = _mm256_load_si256((const __m256i*) & coeff[0]);
  const __m256i v_coeff_1 = _mm256_load_si256((const __m256i*) & coeff[16]);

  const __m256i* v_src = (const __m256i*)src;

  __m256i v_src_lo_0 = _mm256_unpacklo_epi16(v_src[0], v_src[2]);
  __m256i v_src_lo_1 = _mm256_unpacklo_epi16(v_src[1], v_src[3]);
  __m256i v_src_hi_0 = _mm256_unpackhi_epi16(v_src[0], v_src[2]);
  __m256i v_src_hi_1 = _mm256_unpackhi_epi16(v_src[1], v_src[3]);

  __m256i v_trunc_lo_00 = truncate_avx2(_mm256_madd_epi16(v_src_lo_0, v_coeff_0), debias, shift);
  __m256i v_trunc_lo_01 = truncate_avx2(_mm256_madd_epi16(v_src_lo_0, v_coeff_1), debias, shift);
  __m256i v_trunc_lo_10 = truncate_avx2(_mm256_madd_epi16(v_src_lo_1, v_coeff_0), debias, shift);
  __m256i v_trunc_lo_11 = truncate_avx2(_mm256_madd_epi16(v_src_lo_1, v_coeff_1), debias, shift);

  __m256i v_trunc_hi_00 = truncate_avx2(_mm256_madd_epi16(v_src_hi_0, v_coeff_0), debias, shift);
  __m256i v_trunc_hi_01 = truncate_avx2(_mm256_madd_epi16(v_src_hi_0, v_coeff_1), debias, shift);
  __m256i v_trunc_hi_10 = truncate_avx2(_mm256_madd_epi16(v_src_hi_1, v_coeff_0), debias, shift);
  __m256i v_trunc_hi_11 = truncate_avx2(_mm256_madd_epi16(v_src_hi_1, v_coeff_1), debias, shift);

  __m256i v_tmp0 = _mm256_packs_epi32(v_trunc_lo_00, v_trunc_lo_01);
  __m256i v_tmp1 = _mm256_packs_epi32(v_trunc_hi_00, v_trunc_hi_01);
  __m256i v_tmp2 = _mm256_packs_epi32(v_trunc_lo_10, v_trunc_lo_11);
  __m256i v_tmp3 = _mm256_packs_epi32(v_trunc_hi_10, v_trunc_hi_11);

  dst[0] = _mm256_permute2x128_si256(v_tmp0, v_tmp1, 0x20);
  dst[1] = _mm256_permute2x128_si256(v_tmp0, v_tmp1, 0x31);
  dst[2] = _mm256_permute2x128_si256(v_tmp2, v_tmp3, 0x20);
  dst[3] = _mm256_permute2x128_si256(v_tmp2, v_tmp3, 0x31);
}

static void fast_inverse_tr_32x2_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = (shift > 0) ? (1 << (shift - 1)) : 0;
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;

  __m256i v_src[4];
  for (int i = 0; i < 4; ++i) {
    v_src[i] = _mm256_permute4x64_epi64(src[i], _MM_SHUFFLE(3, 1, 2, 0));
  }

  __m256i v_add[32];
  for (int i = 0; i < 32; ++i) {
    __m256i v_madd_0 = _mm256_madd_epi16(v_src[0], v_coeff[0]);
    __m256i v_madd_1 = _mm256_madd_epi16(v_src[1], v_coeff[1]);
    __m256i v_madd_2 = _mm256_madd_epi16(v_src[2], v_coeff[2]);
    __m256i v_madd_3 = _mm256_madd_epi16(v_src[3], v_coeff[3]);

    __m256i v_add_00 = _mm256_add_epi32(v_madd_0, v_madd_1);
    __m256i v_add_01 = _mm256_add_epi32(v_madd_2, v_madd_3);

    v_add[i] = _mm256_add_epi32(v_add_00, v_add_01);
    v_coeff += 4;
  }

  __m256i v_hadd_0[16];
  for (int src = 0, dst = 0; dst < 16; ++dst, src += 2) {
    v_hadd_0[dst] = _mm256_hadd_epi32(v_add[src + 0], v_add[src + 1]);
  }

  __m256i v_hadd_1[8];
  for (int src = 0, dst = 0; dst < 8; ++dst, src += 2) {
    v_hadd_1[dst] = _mm256_hadd_epi32(v_hadd_0[src + 0], v_hadd_0[src + 1]);
  }

  __m256i v_trunc[8];
  for (int i = 0; i < 8; ++i) {
    v_trunc[i] = truncate_avx2(v_hadd_1[i], debias, shift);
  }

  __m256i v_result[4];
  __m256i v_tmp0 = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
  __m256i v_tmp1 = _mm256_packs_epi32(v_trunc[2], v_trunc[3]);
  __m256i v_tmp2 = _mm256_packs_epi32(v_trunc[4], v_trunc[5]);
  __m256i v_tmp3 = _mm256_packs_epi32(v_trunc[6], v_trunc[7]);

  v_result[0] = _mm256_permute2x128_si256(v_tmp0, v_tmp1, 0x20);
  v_result[1] = _mm256_permute2x128_si256(v_tmp2, v_tmp3, 0x20);
  v_result[2] = _mm256_permute2x128_si256(v_tmp0, v_tmp1, 0x31);
  v_result[3] = _mm256_permute2x128_si256(v_tmp2, v_tmp3, 0x31);

  for (int i = 0; i < 4; ++i) {
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }

  // TODO: cutoff for DCT8 and DST7
}

static void fast_inverse_tr_32x2_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 32;
  const int height = 2;

  int skip_width = 0; // DST7 and DCT8 are not defined for this size. Therefore no skip width needed.
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = ff_dct2_2xN_coeff_hor; // TODO: rename
  const int16_t* hor_coeff = fi_dct2_2x32_coeff_ver; // rename
  // No DST7 and DCT8 tables needed.

  __m256i v_ver_pass_out[4];
  fast_inverse_tr_32x2_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);
  
  fast_inverse_tr_32x2_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static void fast_forward_tr_32x4_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 32;
  const int height = 4;
  
  int skip_width = (ver != DCT2 && width == 32) ? 16 : width > 32 ? width - 32 : 0;
  int skip_height = 0;

  const int  reduced_line = width - skip_width;
  const int  cutoff = height - skip_height;
  int16_t* p_dst = dst;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_32xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_32x4_coeff_ver;
  if (hor == DST7) {
    hor_coeff = ff_dst7_32xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = ff_dct8_32xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = ff_dst7_32x4_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = ff_dct8_32x4_coeff_ver;
  }

  __m256i v_hor_pass_out[8];
  if(hor == DCT2) {
    fast_forward_DCT2_B32_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);
  }
  else {
    fast_forward_DCT8_B32_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);    
  }

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)ver_coeff;

  // Got samples for 8 vectors. 4 lines with 32 samples each. Need 2 vectors for each line
  // Handle two lines at a time
  __m256i v_madd_lo_even[2][4];
  __m256i  v_madd_lo_odd[2][4];
  __m256i v_madd_hi_even[2][4];
  __m256i  v_madd_hi_odd[2][4];
  __m256i* v_src_ptr = v_hor_pass_out;
  for (int i = 0; i < 2; ++i) {
    __m256i v_src_hi_0 = _mm256_unpackhi_epi16(v_src_ptr[0], v_src_ptr[2]);
    __m256i v_src_lo_0 = _mm256_unpacklo_epi16(v_src_ptr[0], v_src_ptr[2]);
    __m256i v_src_hi_1 = _mm256_unpackhi_epi16(v_src_ptr[1], v_src_ptr[3]);
    __m256i v_src_lo_1 = _mm256_unpacklo_epi16(v_src_ptr[1], v_src_ptr[3]);

    // Apply coeffs
    for (int ii = 0; ii < 4; ++ii) {
      v_madd_lo_even[i][ii] = _mm256_madd_epi16(v_src_lo_0, v_coeff[ii]);
      v_madd_hi_even[i][ii] = _mm256_madd_epi16(v_src_hi_0, v_coeff[ii]);
      v_madd_lo_odd[i][ii]  = _mm256_madd_epi16(v_src_lo_1, v_coeff[ii]);
      v_madd_hi_odd[i][ii]  = _mm256_madd_epi16(v_src_hi_1, v_coeff[ii]);
    }

    v_coeff += 4;
    v_src_ptr += 4;
  }

  // Final add and truncate
  __m256i v_trunc_lo_even[4];
  __m256i v_trunc_hi_even[4];
  __m256i v_trunc_lo_odd[4];
  __m256i v_trunc_hi_odd[4];
  for (int ii = 0; ii < 4; ++ii) {
    v_trunc_lo_even[ii] = truncate_avx2(_mm256_add_epi32(v_madd_lo_even[0][ii], v_madd_lo_even[1][ii]), debias, shift_2nd);
    v_trunc_lo_odd[ii]  = truncate_avx2(_mm256_add_epi32( v_madd_lo_odd[0][ii],  v_madd_lo_odd[1][ii]), debias, shift_2nd);
    v_trunc_hi_even[ii] = truncate_avx2(_mm256_add_epi32(v_madd_hi_even[0][ii], v_madd_hi_even[1][ii]), debias, shift_2nd);
    v_trunc_hi_odd[ii]  = truncate_avx2(_mm256_add_epi32( v_madd_hi_odd[0][ii],  v_madd_hi_odd[1][ii]), debias, shift_2nd);
  }

  // Permute and store
  for (int i = 0; i < 4; ++i) {
    __m256i v_result_even = _mm256_packs_epi32(v_trunc_lo_even[i], v_trunc_hi_even[i]);
    __m256i v_result_odd  = _mm256_packs_epi32(v_trunc_lo_odd[i], v_trunc_hi_odd[i]);
    // Flip the middle 64 bit chunks
    v_result_even = _mm256_permute4x64_epi64(v_result_even, _MM_SHUFFLE(3, 1, 2, 0));
    v_result_odd = _mm256_permute4x64_epi64(v_result_odd, _MM_SHUFFLE(3, 1, 2, 0));
    _mm256_store_si256((__m256i*)dst, v_result_even);
    _mm256_store_si256((__m256i*)(dst + 16), v_result_odd);
    dst += 32;
  }

}


static void fast_inverse_tr_32x4_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;

  const __m256i* v_src_raw = (const __m256i*)src;

  __m256i v_madd_lo_0[2][4];
  __m256i v_madd_lo_1[2][4];
  __m256i v_madd_hi_0[2][4];
  __m256i v_madd_hi_1[2][4];
  const __m256i* v_c_ptr = v_coeff;
  for (int src = 0; src < 2; ++src) {
    __m256i v_src_lo_0 = _mm256_unpacklo_epi16(v_src_raw[0], v_src_raw[2]);
    __m256i v_src_lo_1 = _mm256_unpacklo_epi16(v_src_raw[1], v_src_raw[3]);
    __m256i v_src_hi_0 = _mm256_unpackhi_epi16(v_src_raw[0], v_src_raw[2]);
    __m256i v_src_hi_1 = _mm256_unpackhi_epi16(v_src_raw[1], v_src_raw[3]);

    for (int i = 0; i < 4; i++) {
      v_madd_lo_0[src][i] = _mm256_madd_epi16(v_src_lo_0, v_c_ptr[i]);
      v_madd_lo_1[src][i] = _mm256_madd_epi16(v_src_lo_1, v_c_ptr[i]);
      v_madd_hi_0[src][i] = _mm256_madd_epi16(v_src_hi_0, v_c_ptr[i]);
      v_madd_hi_1[src][i] = _mm256_madd_epi16(v_src_hi_1, v_c_ptr[i]);
    }
    v_c_ptr += 4;
    v_src_raw += 4;
  }

  __m256i v_trunc_lo[8];
  __m256i v_trunc_hi[8];
  for (int dst = 0, src = 0; src < 4; ++src, dst += 2) {
    v_trunc_lo[dst + 0] = truncate_avx2(_mm256_add_epi32(v_madd_lo_0[0][src], v_madd_lo_0[1][src]), debias, shift);
    v_trunc_lo[dst + 1] = truncate_avx2(_mm256_add_epi32(v_madd_lo_1[0][src], v_madd_lo_1[1][src]), debias, shift);
    v_trunc_hi[dst + 0] = truncate_avx2(_mm256_add_epi32(v_madd_hi_0[0][src], v_madd_hi_0[1][src]), debias, shift);
    v_trunc_hi[dst + 1] = truncate_avx2(_mm256_add_epi32(v_madd_hi_1[0][src], v_madd_hi_1[1][src]), debias, shift);
  }

  dst[0] = _mm256_packs_epi32(v_trunc_lo[0], v_trunc_lo[2]);
  dst[2] = _mm256_packs_epi32(v_trunc_hi[0], v_trunc_hi[2]);
  dst[4] = _mm256_packs_epi32(v_trunc_lo[4], v_trunc_lo[6]);
  dst[6] = _mm256_packs_epi32(v_trunc_hi[4], v_trunc_hi[6]);

  if(skip_line == 0) {
    dst[1] = _mm256_packs_epi32(v_trunc_lo[1], v_trunc_lo[3]);
    dst[3] = _mm256_packs_epi32(v_trunc_hi[1], v_trunc_hi[3]);
    dst[5] = _mm256_packs_epi32(v_trunc_lo[5], v_trunc_lo[7]);
    dst[7] = _mm256_packs_epi32(v_trunc_hi[5], v_trunc_hi[7]);
  }
  else {
    dst[1] = _mm256_setzero_si256();
    dst[3] = _mm256_setzero_si256();
    dst[5] = _mm256_setzero_si256();
    dst[7] = _mm256_setzero_si256();
  }

  // TODO: mts cutoff
}
static void fast_inverse_tr_32x4_avx2_mts_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2) {
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;

  const __m256i* v_src_raw = (const __m256i*)src;

  __m256i v_madd_lo_0[2][4];
  __m256i v_madd_hi_0[2][4];
  const __m256i* v_c_ptr = v_coeff;
  for (int src = 0; src < 2; ++src) {
    __m256i v_src_lo_0 = _mm256_unpacklo_epi16(v_src_raw[0], v_src_raw[2]);
    __m256i v_src_hi_0 = _mm256_unpackhi_epi16(v_src_raw[0], v_src_raw[2]);

    for (int i = 0; i < 4; i++) {
      v_madd_lo_0[src][i] = _mm256_madd_epi16(v_src_lo_0, v_c_ptr[i]);
      v_madd_hi_0[src][i] = _mm256_madd_epi16(v_src_hi_0, v_c_ptr[i]);
    }
    v_c_ptr += 4;
    v_src_raw += 4;
  }

  __m256i v_trunc_lo[4];
  __m256i v_trunc_hi[4];
  for (int src = 0; src < 4; ++src) {
    v_trunc_lo[src] = truncate_avx2(_mm256_add_epi32(v_madd_lo_0[0][src], v_madd_lo_0[1][src]), debias, shift);
    v_trunc_hi[src] = truncate_avx2(_mm256_add_epi32(v_madd_hi_0[0][src], v_madd_hi_0[1][src]), debias, shift);
  }

  dst[0] = _mm256_packs_epi32(v_trunc_lo[0], v_trunc_lo[1]);
  dst[2] = _mm256_packs_epi32(v_trunc_hi[0], v_trunc_hi[1]);
  dst[4] = _mm256_packs_epi32(v_trunc_lo[2], v_trunc_lo[3]);
  dst[6] = _mm256_packs_epi32(v_trunc_hi[2], v_trunc_hi[3]);

  dst[1] = _mm256_setzero_si256();
  dst[3] = _mm256_setzero_si256();
  dst[5] = _mm256_setzero_si256();
  dst[7] = _mm256_setzero_si256();
  

  // TODO: mts cutoff
}

static void fast_inverse_tr_32x4_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const int64_t* c_ptr = (const int64_t*)coeff; // Cast to 64 bit integer to read four coeffs at a time
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)fi_tr_4x8_result_shuffle_ver); // Can use existing shuffle vector

  __m256i v_src[8];
  v_src[0] = _mm256_permute2x128_si256(src[0], src[4], 0x20);
  v_src[1] = _mm256_permute2x128_si256(src[2], src[6], 0x20);
  v_src[2] = _mm256_permute2x128_si256(src[0], src[4], 0x31);
  v_src[3] = _mm256_permute2x128_si256(src[2], src[6], 0x31);

  v_src[4] = _mm256_permute2x128_si256(src[1], src[5], 0x20);
  v_src[5] = _mm256_permute2x128_si256(src[3], src[7], 0x20);
  v_src[6] = _mm256_permute2x128_si256(src[1], src[5], 0x31);
  v_src[7] = _mm256_permute2x128_si256(src[3], src[7], 0x31);

  __m256i v_add[32];
  for (int i = 0; i < 32; ++i) {
    __m256i v_coeff_0 = _mm256_set1_epi64x(c_ptr[0]);
    __m256i v_coeff_1 = _mm256_set1_epi64x(c_ptr[1]);
    __m256i v_coeff_2 = _mm256_set1_epi64x(c_ptr[2]);
    __m256i v_coeff_3 = _mm256_set1_epi64x(c_ptr[3]);
    __m256i v_coeff_4 = _mm256_set1_epi64x(c_ptr[4]);
    __m256i v_coeff_5 = _mm256_set1_epi64x(c_ptr[5]);
    __m256i v_coeff_6 = _mm256_set1_epi64x(c_ptr[6]);
    __m256i v_coeff_7 = _mm256_set1_epi64x(c_ptr[7]);

    __m256i v_madd_0 = _mm256_madd_epi16(v_src[0], v_coeff_0);
    __m256i v_madd_1 = _mm256_madd_epi16(v_src[1], v_coeff_1);
    __m256i v_madd_2 = _mm256_madd_epi16(v_src[2], v_coeff_2);
    __m256i v_madd_3 = _mm256_madd_epi16(v_src[3], v_coeff_3);
    __m256i v_madd_4 = _mm256_madd_epi16(v_src[4], v_coeff_4);
    __m256i v_madd_5 = _mm256_madd_epi16(v_src[5], v_coeff_5);
    __m256i v_madd_6 = _mm256_madd_epi16(v_src[6], v_coeff_6);
    __m256i v_madd_7 = _mm256_madd_epi16(v_src[7], v_coeff_7);

    __m256i v_add_00 = _mm256_add_epi32(v_madd_0, v_madd_1);
    __m256i v_add_01 = _mm256_add_epi32(v_madd_2, v_madd_3);
    __m256i v_add_02 = _mm256_add_epi32(v_madd_4, v_madd_5);
    __m256i v_add_03 = _mm256_add_epi32(v_madd_6, v_madd_7);

    __m256i v_add_10 = _mm256_add_epi32(v_add_00, v_add_01);
    __m256i v_add_11 = _mm256_add_epi32(v_add_02, v_add_03);

    v_add[i] = _mm256_add_epi32(v_add_10, v_add_11);
    c_ptr += 8;
  }

  __m256i v_hadd[16];
  for (int dst = 0, src = 0; dst < 16; ++dst, src += 2) {
    v_hadd[dst] = _mm256_hadd_epi32(v_add[src + 0], v_add[src + 1]);
  }

  __m256i v_trunc[16];
  for (int i = 0; i < 16; ++i) {
    v_trunc[i] = truncate_avx2(v_hadd[i], debias, shift);
  }

  __m256i v_result[8];
  __m256i v_tmp0 = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
  __m256i v_tmp1 = _mm256_packs_epi32(v_trunc[2], v_trunc[3]);
  __m256i v_tmp2 = _mm256_packs_epi32(v_trunc[4], v_trunc[5]);
  __m256i v_tmp3 = _mm256_packs_epi32(v_trunc[6], v_trunc[7]);
  __m256i v_tmp4 = _mm256_packs_epi32(v_trunc[8], v_trunc[9]);
  __m256i v_tmp5 = _mm256_packs_epi32(v_trunc[10], v_trunc[11]);
  __m256i v_tmp6 = _mm256_packs_epi32(v_trunc[12], v_trunc[13]);
  __m256i v_tmp7 = _mm256_packs_epi32(v_trunc[14], v_trunc[15]);

  v_tmp0 = _mm256_shuffle_epi8(v_tmp0, v_res_shuffle);
  v_tmp1 = _mm256_shuffle_epi8(v_tmp1, v_res_shuffle);
  v_tmp2 = _mm256_shuffle_epi8(v_tmp2, v_res_shuffle);
  v_tmp3 = _mm256_shuffle_epi8(v_tmp3, v_res_shuffle);
  v_tmp4 = _mm256_shuffle_epi8(v_tmp4, v_res_shuffle);
  v_tmp5 = _mm256_shuffle_epi8(v_tmp5, v_res_shuffle);
  v_tmp6 = _mm256_shuffle_epi8(v_tmp6, v_res_shuffle);
  v_tmp7 = _mm256_shuffle_epi8(v_tmp7, v_res_shuffle);

  __m256i v_tmp_lo_0 = _mm256_unpacklo_epi64(v_tmp0, v_tmp1);
  __m256i v_tmp_lo_1 = _mm256_unpacklo_epi64(v_tmp2, v_tmp3);
  __m256i v_tmp_lo_2 = _mm256_unpacklo_epi64(v_tmp4, v_tmp5);
  __m256i v_tmp_lo_3 = _mm256_unpacklo_epi64(v_tmp6, v_tmp7);
  __m256i v_tmp_hi_0 = _mm256_unpackhi_epi64(v_tmp0, v_tmp1);
  __m256i v_tmp_hi_1 = _mm256_unpackhi_epi64(v_tmp2, v_tmp3);
  __m256i v_tmp_hi_2 = _mm256_unpackhi_epi64(v_tmp4, v_tmp5);
  __m256i v_tmp_hi_3 = _mm256_unpackhi_epi64(v_tmp6, v_tmp7);

  v_result[0] = _mm256_permute2x128_si256(v_tmp_lo_0, v_tmp_lo_1, 0x20);
  v_result[1] = _mm256_permute2x128_si256(v_tmp_lo_2, v_tmp_lo_3, 0x20);
  v_result[2] = _mm256_permute2x128_si256(v_tmp_hi_0, v_tmp_hi_1, 0x20);
  v_result[3] = _mm256_permute2x128_si256(v_tmp_hi_2, v_tmp_hi_3, 0x20);

  v_result[4] = _mm256_permute2x128_si256(v_tmp_lo_0, v_tmp_lo_1, 0x31);
  v_result[5] = _mm256_permute2x128_si256(v_tmp_lo_2, v_tmp_lo_3, 0x31);
  v_result[6] = _mm256_permute2x128_si256(v_tmp_hi_0, v_tmp_hi_1, 0x31);
  v_result[7] = _mm256_permute2x128_si256(v_tmp_hi_2, v_tmp_hi_3, 0x31);

  for (int i = 0; i < 8; ++i) {
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }
  // TODO: cutoff for dct8 and dst7
}
static void fast_inverse_tr_32x4_avx2_mts_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2) {
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const int64_t* c_ptr = (const int64_t*)coeff; // Cast to 64 bit integer to read four coeffs at a time
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)fi_tr_4x8_result_shuffle_ver); // Can use existing shuffle vector

  __m256i v_src[8];
  v_src[0] = _mm256_permute2x128_si256(src[0], src[4], 0x20);
  v_src[1] = _mm256_permute2x128_si256(src[2], src[6], 0x20);
  v_src[2] = _mm256_permute2x128_si256(src[0], src[4], 0x31);
  v_src[3] = _mm256_permute2x128_si256(src[2], src[6], 0x31);
  

  __m256i v_add[32];
  for (int i = 0; i < 32; ++i) {
    __m256i v_coeff_0 = _mm256_set1_epi64x(c_ptr[0]);
    __m256i v_coeff_1 = _mm256_set1_epi64x(c_ptr[1]);
    __m256i v_coeff_2 = _mm256_set1_epi64x(c_ptr[2]);
    __m256i v_coeff_3 = _mm256_set1_epi64x(c_ptr[3]);

    __m256i v_madd_0 = _mm256_madd_epi16(v_src[0], v_coeff_0);
    __m256i v_madd_1 = _mm256_madd_epi16(v_src[1], v_coeff_1);
    __m256i v_madd_2 = _mm256_madd_epi16(v_src[2], v_coeff_2);
    __m256i v_madd_3 = _mm256_madd_epi16(v_src[3], v_coeff_3);

    __m256i v_add_00 = _mm256_add_epi32(v_madd_0, v_madd_1);
    __m256i v_add_01 = _mm256_add_epi32(v_madd_2, v_madd_3);

    __m256i v_add_10 = _mm256_add_epi32(v_add_00, v_add_01);

    v_add[i] = v_add_10;
    c_ptr += 8;
  }

  __m256i v_hadd[16];
  for (int dst = 0, src = 0; dst < 16; ++dst, src += 2) {
    v_hadd[dst] = _mm256_hadd_epi32(v_add[src + 0], v_add[src + 1]);
  }

  __m256i v_trunc[16];
  for (int i = 0; i < 16; ++i) {
    v_trunc[i] = truncate_avx2(v_hadd[i], debias, shift);
  }

  __m256i v_result[8];
  __m256i v_tmp0 = _mm256_packs_epi32(v_trunc[0], v_trunc[1]);
  __m256i v_tmp1 = _mm256_packs_epi32(v_trunc[2], v_trunc[3]);
  __m256i v_tmp2 = _mm256_packs_epi32(v_trunc[4], v_trunc[5]);
  __m256i v_tmp3 = _mm256_packs_epi32(v_trunc[6], v_trunc[7]);
  __m256i v_tmp4 = _mm256_packs_epi32(v_trunc[8], v_trunc[9]);
  __m256i v_tmp5 = _mm256_packs_epi32(v_trunc[10], v_trunc[11]);
  __m256i v_tmp6 = _mm256_packs_epi32(v_trunc[12], v_trunc[13]);
  __m256i v_tmp7 = _mm256_packs_epi32(v_trunc[14], v_trunc[15]);

  v_tmp0 = _mm256_shuffle_epi8(v_tmp0, v_res_shuffle);
  v_tmp1 = _mm256_shuffle_epi8(v_tmp1, v_res_shuffle);
  v_tmp2 = _mm256_shuffle_epi8(v_tmp2, v_res_shuffle);
  v_tmp3 = _mm256_shuffle_epi8(v_tmp3, v_res_shuffle);
  v_tmp4 = _mm256_shuffle_epi8(v_tmp4, v_res_shuffle);
  v_tmp5 = _mm256_shuffle_epi8(v_tmp5, v_res_shuffle);
  v_tmp6 = _mm256_shuffle_epi8(v_tmp6, v_res_shuffle);
  v_tmp7 = _mm256_shuffle_epi8(v_tmp7, v_res_shuffle);

  __m256i v_tmp_lo_0 = _mm256_unpacklo_epi64(v_tmp0, v_tmp1);
  __m256i v_tmp_lo_1 = _mm256_unpacklo_epi64(v_tmp2, v_tmp3);
  __m256i v_tmp_lo_2 = _mm256_unpacklo_epi64(v_tmp4, v_tmp5);
  __m256i v_tmp_lo_3 = _mm256_unpacklo_epi64(v_tmp6, v_tmp7);
  __m256i v_tmp_hi_0 = _mm256_unpackhi_epi64(v_tmp0, v_tmp1);
  __m256i v_tmp_hi_1 = _mm256_unpackhi_epi64(v_tmp2, v_tmp3);
  __m256i v_tmp_hi_2 = _mm256_unpackhi_epi64(v_tmp4, v_tmp5);
  __m256i v_tmp_hi_3 = _mm256_unpackhi_epi64(v_tmp6, v_tmp7);

  v_result[0] = _mm256_permute2x128_si256(v_tmp_lo_0, v_tmp_lo_1, 0x20);
  v_result[1] = _mm256_permute2x128_si256(v_tmp_lo_2, v_tmp_lo_3, 0x20);
  v_result[2] = _mm256_permute2x128_si256(v_tmp_hi_0, v_tmp_hi_1, 0x20);
  v_result[3] = _mm256_permute2x128_si256(v_tmp_hi_2, v_tmp_hi_3, 0x20);

  v_result[4] = _mm256_permute2x128_si256(v_tmp_lo_0, v_tmp_lo_1, 0x31);
  v_result[5] = _mm256_permute2x128_si256(v_tmp_lo_2, v_tmp_lo_3, 0x31);
  v_result[6] = _mm256_permute2x128_si256(v_tmp_hi_0, v_tmp_hi_1, 0x31);
  v_result[7] = _mm256_permute2x128_si256(v_tmp_hi_2, v_tmp_hi_3, 0x31);

  for (int i = 0; i < 8; ++i) {
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }
  // TODO: cutoff for dct8 and dst7
}

static void fast_inverse_tr_32x4_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 32;
  const int height = 4;

  int skip_width = (hor != DCT2 && width == 32) ? 16 : width > 32 ? width - 32 : 0; 
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = fi_dct2_4x32_coeff_hor; // TODO: rename
  const int16_t* hor_coeff = &uvg_g_dct_32_t[0][0];
  if (hor == DST7) {
    hor_coeff = &uvg_g_dst7_32_t[0][0];
  } else if (hor == DCT8) {
    hor_coeff = &uvg_g_dct8_32[0][0];
  }
  if (ver == DST7) {
    ver_coeff = fi_dst7_4x32_coeff_hor; // TODO: rename
  } else if (ver == DCT8) {
    ver_coeff = fi_dct8_4x32_coeff_hor; // TODO: rename
  }

  __m256i v_ver_pass_out[8];
  if(ver == DCT2) {
    fast_inverse_tr_32x4_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);
  }
  else {
    fast_inverse_tr_32x4_avx2_mts_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);
  }

  if(hor == DCT2) {
    fast_inverse_tr_32x4_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
  }
  else {
    fast_inverse_tr_32x4_avx2_mts_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
  }
}


static void fast_forward_tr_32x8_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 32;
  const int height = 8;
  
  int skip_width = (hor != DCT2 && width == 32) ? 16 : width > 32 ? width - 32 : 0;
  int skip_height = 0;

  const int  reduced_line = width - skip_width;
  const int  cutoff = height - skip_height;
  int16_t* p_dst = dst;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_32xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_32x8_coeff_ver;
  if (hor == DST7) {
    hor_coeff = ff_dst7_32xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = ff_dct8_32xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = ff_dst7_32x8_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = ff_dct8_32x8_coeff_ver;
  }

  __m256i v_hor_pass_out[16];
  if (hor == DCT2) {
    fast_forward_DCT2_B32_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);
  }
  else {
    fast_forward_DCT8_B32_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);
  }

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);
  // Same as for the other 32 and other dimension 8 or 16
  // However all 1,2,4 seem to be producing similar results as with increasing the value
  // just shifts the pressure from one point to another
#define NUM_PARTS 4
#define PART_DIMENSION (8 / NUM_PARTS)
  for (int part = 0; part < NUM_PARTS; ++part) {
    // Got data for 16 vectors, 8 lines 32 samples each
    // Handle two lines at a time
    __m256i v_madd_lo_even[4][PART_DIMENSION];
    __m256i  v_madd_lo_odd[4][PART_DIMENSION];
    __m256i v_madd_hi_even[4][PART_DIMENSION];
    __m256i  v_madd_hi_odd[4][PART_DIMENSION];
    __m256i* v_src_ptr = v_hor_pass_out;
    const __m256i* v_coeff = (const __m256i*)ver_coeff + part * PART_DIMENSION;
    for (int i = 0; i < 4; ++i) {
      __m256i v_src_hi_0 = _mm256_unpackhi_epi16(v_src_ptr[0], v_src_ptr[2]);
      __m256i v_src_lo_0 = _mm256_unpacklo_epi16(v_src_ptr[0], v_src_ptr[2]);
      __m256i v_src_hi_1 = _mm256_unpackhi_epi16(v_src_ptr[1], v_src_ptr[3]);
      __m256i v_src_lo_1 = _mm256_unpacklo_epi16(v_src_ptr[1], v_src_ptr[3]);

      // Apply coeffs
      for (int ii = 0; ii < PART_DIMENSION; ++ii) {
        v_madd_lo_even[i][ii] = _mm256_madd_epi16(v_src_lo_0, v_coeff[ii]);
        v_madd_hi_even[i][ii] = _mm256_madd_epi16(v_src_hi_0, v_coeff[ii]);
        v_madd_lo_odd[i][ii]  = _mm256_madd_epi16(v_src_lo_1, v_coeff[ii]);
        v_madd_hi_odd[i][ii]  = _mm256_madd_epi16(v_src_hi_1, v_coeff[ii]);
      }

      v_coeff += 8;
      v_src_ptr += 4;
    }

    // First round of additions
    __m256i v_add_lo_even[2][PART_DIMENSION];
    __m256i v_add_hi_even[2][PART_DIMENSION];
    __m256i  v_add_lo_odd[2][PART_DIMENSION];
    __m256i  v_add_hi_odd[2][PART_DIMENSION];
    for (int i = 0; i < 2; ++i) {
      const int offset = 2 * i;
      for (int ii = 0; ii < PART_DIMENSION; ++ii) {
        v_add_lo_even[i][ii] = _mm256_add_epi32(v_madd_lo_even[offset][ii], v_madd_lo_even[offset + 1][ii]);
        v_add_hi_even[i][ii] = _mm256_add_epi32(v_madd_hi_even[offset][ii], v_madd_hi_even[offset + 1][ii]);
        v_add_lo_odd[i][ii]  = _mm256_add_epi32(v_madd_lo_odd[offset][ii], v_madd_lo_odd[offset + 1][ii]);
        v_add_hi_odd[i][ii]  = _mm256_add_epi32(v_madd_hi_odd[offset][ii], v_madd_hi_odd[offset + 1][ii]);
      }
    }

    // Final add and truncate
    __m256i v_trunc_lo_even[PART_DIMENSION];
    __m256i v_trunc_hi_even[PART_DIMENSION];
    __m256i v_trunc_lo_odd[PART_DIMENSION];
    __m256i v_trunc_hi_odd[PART_DIMENSION];
    for (int ii = 0; ii < PART_DIMENSION; ++ii) {
      v_trunc_lo_even[ii] = truncate_avx2(_mm256_add_epi32(v_add_lo_even[0][ii], v_add_lo_even[1][ii]), debias, shift_2nd);
      v_trunc_hi_even[ii] = truncate_avx2(_mm256_add_epi32(v_add_hi_even[0][ii], v_add_hi_even[1][ii]), debias, shift_2nd);
      v_trunc_lo_odd[ii]  = truncate_avx2(_mm256_add_epi32(v_add_lo_odd[0][ii], v_add_lo_odd[1][ii]), debias, shift_2nd);
      v_trunc_hi_odd[ii]  = truncate_avx2(_mm256_add_epi32(v_add_hi_odd[0][ii], v_add_hi_odd[1][ii]), debias, shift_2nd);
    }

    // Permute and store
    for (int i = 0; i < PART_DIMENSION; ++i) {
      __m256i v_result_even = _mm256_packs_epi32(v_trunc_lo_even[i], v_trunc_hi_even[i]);
      __m256i v_result_odd = _mm256_packs_epi32(v_trunc_lo_odd[i], v_trunc_hi_odd[i]);
      // Flip the middle 64 bit chunks
      v_result_even = _mm256_permute4x64_epi64(v_result_even, _MM_SHUFFLE(3, 1, 2, 0));
      v_result_odd = _mm256_permute4x64_epi64(v_result_odd, _MM_SHUFFLE(3, 1, 2, 0));
      _mm256_store_si256((__m256i*)dst, v_result_even);
      _mm256_store_si256((__m256i*)(dst + 16), v_result_odd);
      dst += 32;
    }
  }
#undef NUM_PARTS
#undef PART_DIMENSION

}


static void fast_inverse_tr_32x8_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i* v_src_raw = (const __m256i*)src;

  __m256i v_src_lo[8];
  __m256i v_src_hi[8];
  for (int d = 0, s = 0; d < 8; d += 2, s += 4) {
    v_src_lo[d + 0] = _mm256_unpacklo_epi16(v_src_raw[s + 0], v_src_raw[s + 2]);
    v_src_lo[d + 1] = _mm256_unpacklo_epi16(v_src_raw[s + 1], v_src_raw[s + 3]);

    v_src_hi[d + 0] = _mm256_unpackhi_epi16(v_src_raw[s + 0], v_src_raw[s + 2]);
    v_src_hi[d + 1] = _mm256_unpackhi_epi16(v_src_raw[s + 1], v_src_raw[s + 3]);
  }

  for (int c = 0; c < 8; ++c) {
    __m256i v_madd_lo_0[4];
    __m256i v_madd_lo_1[4];
    __m256i v_madd_hi_0[4];
    __m256i v_madd_hi_1[4];
    for (int d = 0, s = 0; d < 4; ++d, s += 2) {
      v_madd_lo_0[d] = _mm256_madd_epi16(v_src_lo[s + 0], v_coeff[d]);
      v_madd_lo_1[d] = _mm256_madd_epi16(v_src_lo[s + 1], v_coeff[d]);
      v_madd_hi_0[d] = _mm256_madd_epi16(v_src_hi[s + 0], v_coeff[d]);
      v_madd_hi_1[d] = _mm256_madd_epi16(v_src_hi[s + 1], v_coeff[d]);
    }
    v_coeff += 4;

    __m256i v_add_lo_00 = _mm256_add_epi32(v_madd_lo_0[0], v_madd_lo_0[1]);
    __m256i v_add_lo_01 = _mm256_add_epi32(v_madd_lo_0[2], v_madd_lo_0[3]);
    __m256i v_add_lo_10 = _mm256_add_epi32(v_madd_lo_1[0], v_madd_lo_1[1]);
    __m256i v_add_lo_11 = _mm256_add_epi32(v_madd_lo_1[2], v_madd_lo_1[3]);

    __m256i v_add_hi_00 = _mm256_add_epi32(v_madd_hi_0[0], v_madd_hi_0[1]);
    __m256i v_add_hi_01 = _mm256_add_epi32(v_madd_hi_0[2], v_madd_hi_0[3]);
    __m256i v_add_hi_10 = _mm256_add_epi32(v_madd_hi_1[0], v_madd_hi_1[1]);
    __m256i v_add_hi_11 = _mm256_add_epi32(v_madd_hi_1[2], v_madd_hi_1[3]);

    __m256i v_trunc_lo_0 = truncate_avx2(_mm256_add_epi32(v_add_lo_00, v_add_lo_01), debias, shift);
    __m256i v_trunc_lo_1 = truncate_avx2(_mm256_add_epi32(v_add_lo_10, v_add_lo_11), debias, shift);

    __m256i v_trunc_hi_0 = truncate_avx2(_mm256_add_epi32(v_add_hi_00, v_add_hi_01), debias, shift);
    __m256i v_trunc_hi_1 = truncate_avx2(_mm256_add_epi32(v_add_hi_10, v_add_hi_11), debias, shift);

    dst[0] = _mm256_packs_epi32(v_trunc_lo_0, v_trunc_hi_0);
    dst[1] = _mm256_packs_epi32(v_trunc_lo_1, v_trunc_hi_1);
    dst += 2;
  }

  // TODO: mts cutoff
}

static void fast_inverse_tr_32x8_avx2_mts_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_coeff = (const __m256i*)coeff;
  const __m256i* v_src_raw = (const __m256i*)src;

  __m256i v_src_lo[4];
  __m256i v_src_hi[4];
  for (int d = 0, s = 0; d < 4; d += 1, s += 4) {
    v_src_lo[d + 0] = _mm256_unpacklo_epi16(v_src_raw[s + 0], v_src_raw[s + 2]);
    v_src_hi[d + 0] = _mm256_unpackhi_epi16(v_src_raw[s + 0], v_src_raw[s + 2]);
  }

  for (int c = 0; c < 8; ++c) {
    __m256i v_madd_lo_0[4];
    __m256i v_madd_hi_0[4];
    for (int d = 0, s = 0; d < 4; ++d, s += 1) {
      v_madd_lo_0[d] = _mm256_madd_epi16(v_src_lo[s + 0], v_coeff[d]);
      v_madd_hi_0[d] = _mm256_madd_epi16(v_src_hi[s + 0], v_coeff[d]);
    }
    v_coeff += 4;

    __m256i v_add_lo_00 = _mm256_add_epi32(v_madd_lo_0[0], v_madd_lo_0[1]);
    __m256i v_add_lo_01 = _mm256_add_epi32(v_madd_lo_0[2], v_madd_lo_0[3]);

    __m256i v_add_hi_00 = _mm256_add_epi32(v_madd_hi_0[0], v_madd_hi_0[1]);
    __m256i v_add_hi_01 = _mm256_add_epi32(v_madd_hi_0[2], v_madd_hi_0[3]);

    __m256i v_trunc_lo_0 = truncate_avx2(_mm256_add_epi32(v_add_lo_00, v_add_lo_01), debias, shift);

    __m256i v_trunc_hi_0 = truncate_avx2(_mm256_add_epi32(v_add_hi_00, v_add_hi_01), debias, shift);

    dst[0] = _mm256_packs_epi32(v_trunc_lo_0, v_trunc_hi_0);
    dst[1] = _mm256_setzero_si256();
    dst += 2;
  }

  // TODO: mts cutoff
}

static void fast_inverse_tr_32x8_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const int limit = skip_line2 == 16 ? 8 : 16;

  int32_t *src_32 = (int32_t*)src;
  for (int j = 0; j < line; ++j) {
    __m256i res_0 = _mm256_setzero_si256();
    __m256i res_1 = _mm256_setzero_si256();
    __m256i res_2 = _mm256_setzero_si256();
    __m256i res_3 = _mm256_setzero_si256();

    __m256i *coeff_start = (__m256i*)coeff;
    for (int i = 0; i < limit; ++i) {
      __m256i v_src = _mm256_set1_epi32(*src_32);
      src_32++;

      __m256i v_coeff0 = _mm256_loadu_si256(coeff_start);
      coeff_start++;
      __m256i v_coeff1 = _mm256_loadu_si256(coeff_start);
      coeff_start++;
      __m256i v_coeff2 = _mm256_loadu_si256(coeff_start);
      coeff_start++;
      __m256i v_coeff3 = _mm256_loadu_si256(coeff_start);
      coeff_start++;

      __m256i madd0 = _mm256_madd_epi16(v_src, v_coeff0);
      __m256i madd1 = _mm256_madd_epi16(v_src, v_coeff1);
      __m256i madd2 = _mm256_madd_epi16(v_src, v_coeff2);
      __m256i madd3 = _mm256_madd_epi16(v_src, v_coeff3);

      res_0 = _mm256_add_epi32(res_0, madd0);
      res_1 = _mm256_add_epi32(res_1, madd1);
      res_2 = _mm256_add_epi32(res_2, madd2);
      res_3 = _mm256_add_epi32(res_3, madd3);
    }
    src_32 += limit == 8 ? 8 : 0;

    __m256i v_trunk0 = truncate_avx2(res_0, debias, shift);
    __m256i v_trunk1 = truncate_avx2(res_1, debias, shift);
    __m256i v_trunk2 = truncate_avx2(res_2, debias, shift);
    __m256i v_trunk3 = truncate_avx2(res_3, debias, shift);

    __m256i packed0 =  _mm256_packs_epi32(v_trunk0, v_trunk1);
    __m256i packed1 =  _mm256_packs_epi32(v_trunk2, v_trunk3);

    packed0 = _mm256_permute4x64_epi64(packed0, _MM_SHUFFLE(3, 1, 2, 0));
    packed1 = _mm256_permute4x64_epi64(packed1, _MM_SHUFFLE(3, 1, 2, 0));

    _mm256_store_si256((__m256i*)dst, packed0);
    _mm256_store_si256((__m256i*)dst + 1, packed1);
    dst += 32;
  }
  
  // TODO: cutoff for dct8 and dst7
}

static void fast_inverse_tr_32x8_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 32;
  const int height = 8;

  int skip_width = (hor != DCT2 && width == 32) ? 16 : width > 32 ? width - 32 : 0;
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = fi_dct2_8x32_coeff_hor; // TODO: rename this table
  const int16_t* hor_coeff = fi_dct2_32xN_coeff_hor;
  if (hor == DST7) {
    hor_coeff = fi_dst7_32xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = fi_dct8_32xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = fi_dst7_8x32_coeff_hor; // TODO: rename
  } else if (ver == DCT8) {
    ver_coeff = fi_dct8_8x32_coeff_hor; // TODO: rename
  }

  __m256i v_ver_pass_out[16];
  if(ver == DCT2 || hor == DCT2) {
    fast_inverse_tr_32x8_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, height, 0, skip_width);
  }
  else {
    fast_inverse_tr_32x8_avx2_mts_ver(src, v_ver_pass_out, ver_coeff, shift_1st, height, 0, skip_width);    
  }

  fast_inverse_tr_32x8_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static void fast_forward_tr_32x16_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 32;
  const int height = 16;
  
  int skip_width = (ver != DCT2 && width == 32) ? 16 : width > 32 ? width - 32 : 0;
  int skip_height = 0;

  const int  reduced_line = width - skip_width;
  const int  cutoff = height - skip_height;
  int16_t* p_dst = dst;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_32xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_32x16_coeff_ver;
  if (hor == DST7) {
    hor_coeff = ff_dst7_32xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = ff_dct8_32xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = ff_dst7_32x16_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = ff_dct8_32x16_coeff_ver;
  }

  __m256i v_hor_pass_out[32];
  if (hor == DCT2) {
    fast_forward_DCT2_B32_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);
  }
  else {
    fast_forward_DCT8_B32_avx2_hor(src, v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);
  }

  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);
  // Same as for 8x32 and 16x32, 4 parts is optimal
#define NUM_PARTS 4
#define PART_DIMENSION (16 / NUM_PARTS)
  for (int part = 0; part < NUM_PARTS; ++part) {
    // Got samples for 32 vectors, 16 lines with 32 samples each
    // Handle two lines at a time
    __m256i v_madd_lo_even[8][PART_DIMENSION];
    __m256i  v_madd_lo_odd[8][PART_DIMENSION];
    __m256i v_madd_hi_even[8][PART_DIMENSION];
    __m256i  v_madd_hi_odd[8][PART_DIMENSION];
    __m256i* v_src_ptr = v_hor_pass_out;
    const int32_t* line_coeff = (const int32_t*)ver_coeff + part * PART_DIMENSION;
    for (int i = 0; i < 8; ++i) {
      __m256i v_src_hi_0 = _mm256_unpackhi_epi16(v_src_ptr[0], v_src_ptr[2]);
      __m256i v_src_lo_0 = _mm256_unpacklo_epi16(v_src_ptr[0], v_src_ptr[2]);
      __m256i v_src_hi_1 = _mm256_unpackhi_epi16(v_src_ptr[1], v_src_ptr[3]);
      __m256i v_src_lo_1 = _mm256_unpacklo_epi16(v_src_ptr[1], v_src_ptr[3]);

      // Apply coeffs
      for (int ii = 0; ii < PART_DIMENSION; ++ii) {
        const int32_t coeff = line_coeff[ii];
        const __m256i v_coeff = _mm256_set1_epi32(coeff);
        v_madd_lo_even[i][ii] = _mm256_madd_epi16(v_src_lo_0, v_coeff);
        v_madd_hi_even[i][ii] = _mm256_madd_epi16(v_src_hi_0, v_coeff);
        v_madd_lo_odd[i][ii] = _mm256_madd_epi16(v_src_lo_1, v_coeff);
        v_madd_hi_odd[i][ii] = _mm256_madd_epi16(v_src_hi_1, v_coeff);
      }

      line_coeff += 16;
      v_src_ptr += 4;
    }

    for (int ii = 0; ii < PART_DIMENSION; ++ii) {
      // First round of additions
      __m256i v_add_lo_even_0[4];
      __m256i v_add_hi_even_0[4];
      __m256i v_add_lo_odd_0[4];
      __m256i v_add_hi_odd_0[4];
      for (int i = 0; i < 4; ++i) {
        const int offset = i * 2;
        v_add_lo_even_0[i] = _mm256_add_epi32(v_madd_lo_even[offset][ii], v_madd_lo_even[offset + 1][ii]);
        v_add_hi_even_0[i] = _mm256_add_epi32(v_madd_hi_even[offset][ii], v_madd_hi_even[offset + 1][ii]);
        v_add_lo_odd_0[i] = _mm256_add_epi32(v_madd_lo_odd[offset][ii], v_madd_lo_odd[offset + 1][ii]);
        v_add_hi_odd_0[i] = _mm256_add_epi32(v_madd_hi_odd[offset][ii], v_madd_hi_odd[offset + 1][ii]);
      }

      // Second round of additions
      __m256i v_add_lo_even_1[2];
      __m256i v_add_hi_even_1[2];
      __m256i v_add_lo_odd_1[2];
      __m256i v_add_hi_odd_1[2];
      for (int i = 0; i < 2; ++i) {
        const int offset = 2 * i;
        v_add_lo_even_1[i] = _mm256_add_epi32(v_add_lo_even_0[offset], v_add_lo_even_0[offset + 1]);
        v_add_hi_even_1[i] = _mm256_add_epi32(v_add_hi_even_0[offset], v_add_hi_even_0[offset + 1]);
        v_add_lo_odd_1[i] = _mm256_add_epi32(v_add_lo_odd_0[offset], v_add_lo_odd_0[offset + 1]);
        v_add_hi_odd_1[i] = _mm256_add_epi32(v_add_hi_odd_0[offset], v_add_hi_odd_0[offset + 1]);
      }

      // Final add and truncate
      __m256i v_trunc_lo_even;
      __m256i v_trunc_hi_even;
      __m256i v_trunc_lo_odd;
      __m256i v_trunc_hi_odd;
      v_trunc_lo_even = truncate_avx2(_mm256_add_epi32(v_add_lo_even_1[0], v_add_lo_even_1[1]), debias, shift_2nd);
      v_trunc_hi_even = truncate_avx2(_mm256_add_epi32(v_add_hi_even_1[0], v_add_hi_even_1[1]), debias, shift_2nd);
      v_trunc_lo_odd = truncate_avx2(_mm256_add_epi32(v_add_lo_odd_1[0], v_add_lo_odd_1[1]), debias, shift_2nd);
      v_trunc_hi_odd = truncate_avx2(_mm256_add_epi32(v_add_hi_odd_1[0], v_add_hi_odd_1[1]), debias, shift_2nd);


      // Permute and store
      __m256i v_result_even = _mm256_packs_epi32(v_trunc_lo_even, v_trunc_hi_even);
      __m256i v_result_odd = _mm256_packs_epi32(v_trunc_lo_odd, v_trunc_hi_odd);
      // Flip the middle 64 bit chunks
      v_result_even = _mm256_permute4x64_epi64(v_result_even, _MM_SHUFFLE(3, 1, 2, 0));
      v_result_odd = _mm256_permute4x64_epi64(v_result_odd, _MM_SHUFFLE(3, 1, 2, 0));
      _mm256_store_si256((__m256i*)dst, v_result_even);
      _mm256_store_si256((__m256i*)(dst + 16), v_result_odd);
      dst += 32;
    }
  }
#undef NUM_PARTS
#undef PART_DIMENSION
}


static void fast_inverse_tr_32x16_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const int limit = 32 - skip_line;
  __m256i temp[32];
  for (int j = 0; j < limit; ++j) {
    __m256i res_0 = _mm256_setzero_si256();
    __m256i res_1 = _mm256_setzero_si256();

    __m256i* coeff_start = (__m256i*)coeff;
    for (int i = 0; i < 8; ++i) {
      int16_t source[2];
      source[0] = src[j + i * 64];
      source[1] = src[j + i * 64 + 32];
      int32_t paired_source;
      memcpy(&paired_source, source, sizeof(int32_t));

      __m256i v_src = _mm256_set1_epi32(paired_source);

      __m256i v_coeff0 = _mm256_load_si256(coeff_start);
      coeff_start++;
      __m256i v_coeff1 = _mm256_load_si256(coeff_start);
      coeff_start++;

      __m256i v_madd0 = _mm256_madd_epi16(v_src, v_coeff0);
      __m256i v_madd1 = _mm256_madd_epi16(v_src, v_coeff1);

      res_0 = _mm256_add_epi32(res_0, v_madd0);
      res_1 = _mm256_add_epi32(res_1, v_madd1);
    }

    __m256i v_trunc0 = truncate_avx2(res_0, debias, shift);
    __m256i v_trunc1 = truncate_avx2(res_1, debias, shift);

    __m256i packed = _mm256_packs_epi32(v_trunc0, v_trunc1);
    packed = _mm256_permute4x64_epi64(packed, _MM_SHUFFLE(3, 1, 2, 0));
    temp[j] = packed;
  }
  for (int j = limit; j < 32; ++j) {
    temp[j] = _mm256_setzero_si256();
  }
  transpose_avx2(temp, dst, 16, 32);
}

static void fast_inverse_tr_32x16_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i* v_src_raw = src;
  const __m256i v_res_shuffle = _mm256_load_si256((const __m256i*)shuffle_16b_0415);

  // Do a 32-bit transpose to arrange result from previous pass
  __m256i v_tmp32_lo_e[8];
  __m256i v_tmp32_hi_e[8];
  __m256i v_tmp32_lo_o[8];
  __m256i v_tmp32_hi_o[8];
  for (int d = 0, s = 0; d < 8; ++d, s += 4) {
    v_tmp32_lo_e[d] = _mm256_unpacklo_epi32(v_src_raw[s + 0], v_src_raw[s + 2]);
    v_tmp32_hi_e[d] = _mm256_unpackhi_epi32(v_src_raw[s + 0], v_src_raw[s + 2]);
    v_tmp32_lo_o[d] = _mm256_unpacklo_epi32(v_src_raw[s + 1], v_src_raw[s + 3]);
    v_tmp32_hi_o[d] = _mm256_unpackhi_epi32(v_src_raw[s + 1], v_src_raw[s + 3]);
  }

  __m256i v_tmp64_lo_e[8];
  __m256i v_tmp64_hi_e[8];
  __m256i v_tmp64_lo_o[8];
  __m256i v_tmp64_hi_o[8];
  for (int d = 0, s = 0; d < 4; ++d, s += 2) {
    v_tmp64_lo_e[0 + d] = _mm256_unpacklo_epi64(v_tmp32_lo_e[s + 0], v_tmp32_lo_e[s + 1]);
    v_tmp64_lo_e[4 + d] = _mm256_unpacklo_epi64(v_tmp32_hi_e[s + 0], v_tmp32_hi_e[s + 1]);

    v_tmp64_hi_e[0 + d] = _mm256_unpackhi_epi64(v_tmp32_lo_e[s + 0], v_tmp32_lo_e[s + 1]);
    v_tmp64_hi_e[4 + d] = _mm256_unpackhi_epi64(v_tmp32_hi_e[s + 0], v_tmp32_hi_e[s + 1]);

    v_tmp64_lo_o[0 + d] = _mm256_unpacklo_epi64(v_tmp32_lo_o[s + 0], v_tmp32_lo_o[s + 1]);
    v_tmp64_lo_o[4 + d] = _mm256_unpacklo_epi64(v_tmp32_hi_o[s + 0], v_tmp32_hi_o[s + 1]);

    v_tmp64_hi_o[0 + d] = _mm256_unpackhi_epi64(v_tmp32_lo_o[s + 0], v_tmp32_lo_o[s + 1]);
    v_tmp64_hi_o[4 + d] = _mm256_unpackhi_epi64(v_tmp32_hi_o[s + 0], v_tmp32_hi_o[s + 1]);
  }

  __m256i v_src[32];
  v_src[0] = _mm256_permute2x128_si256(v_tmp64_lo_e[0], v_tmp64_lo_e[1], 0x20);
  v_src[1] = _mm256_permute2x128_si256(v_tmp64_hi_e[0], v_tmp64_hi_e[1], 0x20);
  v_src[2] = _mm256_permute2x128_si256(v_tmp64_lo_e[4], v_tmp64_lo_e[5], 0x20);
  v_src[3] = _mm256_permute2x128_si256(v_tmp64_hi_e[4], v_tmp64_hi_e[5], 0x20);

  v_src[4] = _mm256_permute2x128_si256(v_tmp64_lo_e[0], v_tmp64_lo_e[1], 0x31);
  v_src[5] = _mm256_permute2x128_si256(v_tmp64_hi_e[0], v_tmp64_hi_e[1], 0x31);
  v_src[6] = _mm256_permute2x128_si256(v_tmp64_lo_e[4], v_tmp64_lo_e[5], 0x31);
  v_src[7] = _mm256_permute2x128_si256(v_tmp64_hi_e[4], v_tmp64_hi_e[5], 0x31);

  v_src[8] = _mm256_permute2x128_si256(v_tmp64_lo_o[0], v_tmp64_lo_o[1], 0x20);
  v_src[9] = _mm256_permute2x128_si256(v_tmp64_hi_o[0], v_tmp64_hi_o[1], 0x20);
  v_src[10] = _mm256_permute2x128_si256(v_tmp64_lo_o[4], v_tmp64_lo_o[5], 0x20);
  v_src[11] = _mm256_permute2x128_si256(v_tmp64_hi_o[4], v_tmp64_hi_o[5], 0x20);

  v_src[12] = _mm256_permute2x128_si256(v_tmp64_lo_o[0], v_tmp64_lo_o[1], 0x31);
  v_src[13] = _mm256_permute2x128_si256(v_tmp64_hi_o[0], v_tmp64_hi_o[1], 0x31);
  v_src[14] = _mm256_permute2x128_si256(v_tmp64_lo_o[4], v_tmp64_lo_o[5], 0x31);
  v_src[15] = _mm256_permute2x128_si256(v_tmp64_hi_o[4], v_tmp64_hi_o[5], 0x31);

  v_src[16] = _mm256_permute2x128_si256(v_tmp64_lo_e[2], v_tmp64_lo_e[3], 0x20);
  v_src[17] = _mm256_permute2x128_si256(v_tmp64_hi_e[2], v_tmp64_hi_e[3], 0x20);
  v_src[18] = _mm256_permute2x128_si256(v_tmp64_lo_e[6], v_tmp64_lo_e[7], 0x20);
  v_src[19] = _mm256_permute2x128_si256(v_tmp64_hi_e[6], v_tmp64_hi_e[7], 0x20);

  v_src[20] = _mm256_permute2x128_si256(v_tmp64_lo_e[2], v_tmp64_lo_e[3], 0x31);
  v_src[21] = _mm256_permute2x128_si256(v_tmp64_hi_e[2], v_tmp64_hi_e[3], 0x31);
  v_src[22] = _mm256_permute2x128_si256(v_tmp64_lo_e[6], v_tmp64_lo_e[7], 0x31);
  v_src[23] = _mm256_permute2x128_si256(v_tmp64_hi_e[6], v_tmp64_hi_e[7], 0x31);

  v_src[24] = _mm256_permute2x128_si256(v_tmp64_lo_o[2], v_tmp64_lo_o[3], 0x20);
  v_src[25] = _mm256_permute2x128_si256(v_tmp64_hi_o[2], v_tmp64_hi_o[3], 0x20);
  v_src[26] = _mm256_permute2x128_si256(v_tmp64_lo_o[6], v_tmp64_lo_o[7], 0x20);
  v_src[27] = _mm256_permute2x128_si256(v_tmp64_hi_o[6], v_tmp64_hi_o[7], 0x20);

  v_src[28] = _mm256_permute2x128_si256(v_tmp64_lo_o[2], v_tmp64_lo_o[3], 0x31);
  v_src[29] = _mm256_permute2x128_si256(v_tmp64_hi_o[2], v_tmp64_hi_o[3], 0x31);
  v_src[30] = _mm256_permute2x128_si256(v_tmp64_lo_o[6], v_tmp64_lo_o[7], 0x31);
  v_src[31] = _mm256_permute2x128_si256(v_tmp64_hi_o[6], v_tmp64_hi_o[7], 0x31);

  __m256i v_trunc[64];
  __m256i* v_src_ptr = v_src;
  __m256i* v_tr_ptr = v_trunc;


  for (int chunk = 0; chunk < 2; ++chunk) {
    const int32_t* c_ptr = (const int32_t*)coeff; // Handle as 32 bit integer to load two coeffs into vector at the same time
    for (int c = 0; c < 32; ++c) {
      __m256i v_madd[16];
      for (int i = 0; i < 16; ++i) {
        const __m256i v_coeff = _mm256_set1_epi32(*c_ptr);
        v_madd[i] = _mm256_madd_epi16(v_src_ptr[i], v_coeff);
        c_ptr++;
      }

      __m256i v_add_0[8];
      for (int d = 0, s = 0; d < 8; ++d, s += 2) {
        v_add_0[d] = _mm256_add_epi32(v_madd[s + 0], v_madd[s + 1]);
      }

      __m256i v_add_1[4];
      for (int d = 0, s = 0; d < 4; ++d, s += 2) {
        v_add_1[d] = _mm256_add_epi32(v_add_0[s + 0], v_add_0[s + 1]);
      }

      __m256i v_add_2[2];
      for (int d = 0, s = 0; d < 2; ++d, s += 2) {
        v_add_2[d] = _mm256_add_epi32(v_add_1[s + 0], v_add_1[s + 1]);
      }

      v_tr_ptr[c] = truncate_avx2(_mm256_add_epi32(v_add_2[0], v_add_2[1]), debias, shift);
    }
    v_tr_ptr += 32;
    v_src_ptr += 16;
  }

  __m256i v_tmp[32];
  __m256i v_result[32];
  for (int i = 0, s = 0; i < 32; ++i, s += 2) {
    v_tmp[i] = _mm256_packs_epi32(v_trunc[s + 0], v_trunc[s + 1]);
    v_tmp[i] = _mm256_shuffle_epi8(v_tmp[i], v_res_shuffle);
  }

  __m256i v_rtmp32_lo[16];
  __m256i v_rtmp32_hi[16];
  for (int d = 0, s = 0; d < 16; ++d, s += 2) {
    v_rtmp32_lo[d] = _mm256_unpacklo_epi32(v_tmp[s + 0], v_tmp[s + 1]);
    v_rtmp32_hi[d] = _mm256_unpackhi_epi32(v_tmp[s + 0], v_tmp[s + 1]);
  }

  __m256i v_rtmp64_lo[16];
  __m256i v_rtmp64_hi[16];
  for (int d = 0, s = 0; d < 8; ++d, s += 2) {
    v_rtmp64_lo[0 + d] = _mm256_unpacklo_epi64(v_rtmp32_lo[s + 0], v_rtmp32_lo[s + 1]);
    v_rtmp64_lo[8 + d] = _mm256_unpacklo_epi64(v_rtmp32_hi[s + 0], v_rtmp32_hi[s + 1]);

    v_rtmp64_hi[0 + d] = _mm256_unpackhi_epi64(v_rtmp32_lo[s + 0], v_rtmp32_lo[s + 1]);
    v_rtmp64_hi[8 + d] = _mm256_unpackhi_epi64(v_rtmp32_hi[s + 0], v_rtmp32_hi[s + 1]);
  }

  v_result[0] = _mm256_permute2x128_si256(v_rtmp64_lo[0], v_rtmp64_lo[1], 0x20);
  v_result[1] = _mm256_permute2x128_si256(v_rtmp64_lo[2], v_rtmp64_lo[3], 0x20);
  v_result[2] = _mm256_permute2x128_si256(v_rtmp64_hi[0], v_rtmp64_hi[1], 0x20);
  v_result[3] = _mm256_permute2x128_si256(v_rtmp64_hi[2], v_rtmp64_hi[3], 0x20);

  v_result[4] = _mm256_permute2x128_si256(v_rtmp64_lo[8], v_rtmp64_lo[9], 0x20);
  v_result[5] = _mm256_permute2x128_si256(v_rtmp64_lo[10], v_rtmp64_lo[11], 0x20);
  v_result[6] = _mm256_permute2x128_si256(v_rtmp64_hi[8], v_rtmp64_hi[9], 0x20);
  v_result[7] = _mm256_permute2x128_si256(v_rtmp64_hi[10], v_rtmp64_hi[11], 0x20);

  v_result[8] = _mm256_permute2x128_si256(v_rtmp64_lo[0], v_rtmp64_lo[1], 0x31);
  v_result[9] = _mm256_permute2x128_si256(v_rtmp64_lo[2], v_rtmp64_lo[3], 0x31);
  v_result[10] = _mm256_permute2x128_si256(v_rtmp64_hi[0], v_rtmp64_hi[1], 0x31);
  v_result[11] = _mm256_permute2x128_si256(v_rtmp64_hi[2], v_rtmp64_hi[3], 0x31);

  v_result[12] = _mm256_permute2x128_si256(v_rtmp64_lo[8], v_rtmp64_lo[9], 0x31);
  v_result[13] = _mm256_permute2x128_si256(v_rtmp64_lo[10], v_rtmp64_lo[11], 0x31);
  v_result[14] = _mm256_permute2x128_si256(v_rtmp64_hi[8], v_rtmp64_hi[9], 0x31);
  v_result[15] = _mm256_permute2x128_si256(v_rtmp64_hi[10], v_rtmp64_hi[11], 0x31);

  v_result[16] = _mm256_permute2x128_si256(v_rtmp64_lo[4], v_rtmp64_lo[5], 0x20);
  v_result[17] = _mm256_permute2x128_si256(v_rtmp64_lo[6], v_rtmp64_lo[7], 0x20);
  v_result[18] = _mm256_permute2x128_si256(v_rtmp64_hi[4], v_rtmp64_hi[5], 0x20);
  v_result[19] = _mm256_permute2x128_si256(v_rtmp64_hi[6], v_rtmp64_hi[7], 0x20);

  v_result[20] = _mm256_permute2x128_si256(v_rtmp64_lo[12], v_rtmp64_lo[13], 0x20);
  v_result[21] = _mm256_permute2x128_si256(v_rtmp64_lo[14], v_rtmp64_lo[15], 0x20);
  v_result[22] = _mm256_permute2x128_si256(v_rtmp64_hi[12], v_rtmp64_hi[13], 0x20);
  v_result[23] = _mm256_permute2x128_si256(v_rtmp64_hi[14], v_rtmp64_hi[15], 0x20);

  v_result[24] = _mm256_permute2x128_si256(v_rtmp64_lo[4], v_rtmp64_lo[5], 0x31);
  v_result[25] = _mm256_permute2x128_si256(v_rtmp64_lo[6], v_rtmp64_lo[7], 0x31);
  v_result[26] = _mm256_permute2x128_si256(v_rtmp64_hi[4], v_rtmp64_hi[5], 0x31);
  v_result[27] = _mm256_permute2x128_si256(v_rtmp64_hi[6], v_rtmp64_hi[7], 0x31);

  v_result[28] = _mm256_permute2x128_si256(v_rtmp64_lo[12], v_rtmp64_lo[13], 0x31);
  v_result[29] = _mm256_permute2x128_si256(v_rtmp64_lo[14], v_rtmp64_lo[15], 0x31);
  v_result[30] = _mm256_permute2x128_si256(v_rtmp64_hi[12], v_rtmp64_hi[13], 0x31);
  v_result[31] = _mm256_permute2x128_si256(v_rtmp64_hi[14], v_rtmp64_hi[15], 0x31);

  for (int i = 0; i < 32; ++i) {
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }

  // TODO: MTS cutoff
}

static void fast_inverse_tr_32x16_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 32;
  const int height = 16;

  int skip_width = (hor != DCT2 && width == 32) ? 16 : width > 32 ? width - 32 : 0;
  int skip_height = 0;

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = fi_dct2_32x16_coeff_ver;
  const int16_t* hor_coeff = fi_dct2_32xN_coeff_hor;
  if (hor == DST7) {
    hor_coeff = fi_dst7_32xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = fi_dct8_32xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = fi_dst7_32x16_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = fi_dct8_32x16_coeff_ver;
  }

  __m256i v_ver_pass_out[32];
  fast_inverse_tr_32x16_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);

  fast_inverse_tr_32x8_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static void fast_forward_tr_32x32_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 32;
  const int height = 32;

  int skip_width = (hor != DCT2 && width == 32) ? 16 : width > 32 ? width - 32 : 0;
  int skip_height = (ver != DCT2 && height == 32) ? 16 : (height > 32 ? height - 32 : 0);

  const int  reduced_line = width - skip_width;
  const int  cutoff = height - skip_height;
  int16_t* p_dst = dst;

  const int log2_width_minus1 = uvg_g_convert_to_log2[width] - 1;
  const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
  const int32_t shift_1st = log2_width_minus1 + UVG_BIT_DEPTH - 8;
  const int32_t shift_2nd = log2_height_minus1 + 7;

  const int16_t* hor_coeff = ff_dct2_32xN_coeff_hor;
  const int16_t* ver_coeff = ff_dct2_32x32_coeff_ver;
  if (hor == DST7) {
    hor_coeff = ff_dst7_32xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = ff_dct8_32xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = ff_dst7_32x32_coeff_ver;
  } else if (ver == DCT8) {
    ver_coeff = ff_dct8_32x32_coeff_ver;
  }

  ALIGNED(32) int16_t v_hor_pass_out[32 * 32];
  if(hor == DCT2) {
    fast_forward_DCT2_B32_avx2_hor(src, (__m256i*)v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);
  }
  else {
    fast_forward_DCT8_B32_avx2_hor(src, (__m256i*)v_hor_pass_out, hor_coeff, shift_1st, height, 0, skip_width);    
  }

  __m256i temp_out[32 * 2];
  // Vertical pass
  const int32_t    add = (shift_2nd > 0) ? (1 << (shift_2nd - 1)) : 0; // ISP_TODO: optimize (shift > 0) check out if shift is always gt 0
  const __m256i debias = _mm256_set1_epi32(add);
  for (int j = 0; j < reduced_line; ++j) {
    __m256i res_0 = _mm256_setzero_si256();
    __m256i res_1 = _mm256_setzero_si256();
    __m256i res_2 = _mm256_setzero_si256();
    __m256i res_3 = _mm256_setzero_si256();
    const int16_t* coeff_start = ver_coeff;
    for (int i = 0; i < 16; ++i) {
      int16_t source[2];
      source[0] = v_hor_pass_out[j + i * 64];
      source[1] = v_hor_pass_out[j + i * 64 + 32];
      int32_t paired_source;
      memcpy(&paired_source, source, sizeof(int32_t));

      __m256i v_src = _mm256_set1_epi32(paired_source);
      __m256i v_coeff_0 = _mm256_load_si256((__m256i*) coeff_start);
      coeff_start += 16;
      __m256i v_coeff_1 = _mm256_load_si256((__m256i*) coeff_start);
      __m256i v_coeff_2;
      __m256i v_coeff_3;
      if(skip_height == 0) {
        coeff_start += 16;
        v_coeff_2 = _mm256_load_si256((__m256i*) coeff_start);
        coeff_start += 16;
        v_coeff_3 = _mm256_load_si256((__m256i*) coeff_start);
        coeff_start += 16;
      }
      else {
        coeff_start += 48;
      }

      __m256i madd_0 = _mm256_madd_epi16(v_src, v_coeff_0);
      __m256i madd_1 = _mm256_madd_epi16(v_src, v_coeff_1);
      __m256i madd_2;
      __m256i madd_3;
      if(skip_height == 0) {
        madd_2 = _mm256_madd_epi16(v_src, v_coeff_2);
        madd_3 = _mm256_madd_epi16(v_src, v_coeff_3);
      }

      res_0 = _mm256_add_epi32(res_0, madd_0);
      res_1 = _mm256_add_epi32(res_1, madd_1);
      if(skip_height == 0) {
        res_2 = _mm256_add_epi32(res_2, madd_2);
        res_3 = _mm256_add_epi32(res_3, madd_3);
      }
    }
    __m256i v_trunc_0 = truncate_avx2(res_0, debias, shift_2nd);
    __m256i v_trunc_1 = truncate_avx2(res_1, debias, shift_2nd);
    __m256i v_trunc_2;
    __m256i v_trunc_3;
    if(skip_height == 0) {
      v_trunc_2 = truncate_avx2(res_2, debias, shift_2nd);
      v_trunc_3 = truncate_avx2(res_3, debias, shift_2nd);
    }

    v_trunc_0 = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
    v_trunc_0 = _mm256_permute4x64_epi64(v_trunc_0, _MM_SHUFFLE(3, 1, 2, 0));
    _mm256_store_si256(temp_out + 2 * j, v_trunc_0);
    if(skip_height == 0) {
      v_trunc_2 = _mm256_packs_epi32(v_trunc_2, v_trunc_3);
      v_trunc_2 = _mm256_permute4x64_epi64(v_trunc_2, _MM_SHUFFLE(3, 1, 2, 0));
      _mm256_store_si256(temp_out + 2 * j + 1, v_trunc_2);
    }
  }
  transpose_avx2(temp_out, (__m256i*) dst, 32, 32);
#if 0
  // 8 is probably best, though difference to 16 is not that large
#define NUM_PARTS 8
#define PART_DIMENSION (32 / NUM_PARTS)
  for (int part = 0; part < NUM_PARTS; ++part) {
    const int32_t* coeff_ptr = (const int32_t*)ver_coeff + part * PART_DIMENSION; // Cast to 32 bit integer to read 2 coeffs at a time
    const __m256i* v_src_ptr = v_hor_pass_out;

    __m256i v_madd_lo_e[16][PART_DIMENSION];
    __m256i v_madd_lo_o[16][PART_DIMENSION];
    __m256i v_madd_hi_e[16][PART_DIMENSION];
    __m256i v_madd_hi_o[16][PART_DIMENSION];
    for (int i = 0; i < 16; ++i) {
      __m256i v_src_lo_e = _mm256_unpacklo_epi16(v_src_ptr[0], v_src_ptr[2]);
      __m256i v_src_lo_o = _mm256_unpacklo_epi16(v_src_ptr[1], v_src_ptr[3]);
      __m256i v_src_hi_e = _mm256_unpackhi_epi16(v_src_ptr[0], v_src_ptr[2]);
      __m256i v_src_hi_o = _mm256_unpackhi_epi16(v_src_ptr[1], v_src_ptr[3]);


      for (int c = 0; c < PART_DIMENSION; ++c) {
        const __m256i v_coeff = _mm256_set1_epi32(coeff_ptr[c]);
        v_madd_lo_e[i][c] = _mm256_madd_epi16(v_src_lo_e, v_coeff);
        v_madd_lo_o[i][c] = _mm256_madd_epi16(v_src_lo_o, v_coeff);
        v_madd_hi_e[i][c] = _mm256_madd_epi16(v_src_hi_e, v_coeff);
        v_madd_hi_o[i][c] = _mm256_madd_epi16(v_src_hi_o, v_coeff);
      }
      coeff_ptr += 32;
      v_src_ptr += 4;
    }

    for (int c = 0; c < PART_DIMENSION; ++c) {
      __m256i v_add_lo_e0[8];
      __m256i v_add_lo_o0[8];
      __m256i v_add_hi_e0[8];
      __m256i v_add_hi_o0[8];
      for (int dst = 0, src = 0; dst < 8; ++dst, src += 2) {
        v_add_lo_e0[dst] = _mm256_add_epi32(v_madd_lo_e[src + 0][c], v_madd_lo_e[src + 1][c]);
        v_add_lo_o0[dst] = _mm256_add_epi32(v_madd_lo_o[src + 0][c], v_madd_lo_o[src + 1][c]);
        v_add_hi_e0[dst] = _mm256_add_epi32(v_madd_hi_e[src + 0][c], v_madd_hi_e[src + 1][c]);
        v_add_hi_o0[dst] = _mm256_add_epi32(v_madd_hi_o[src + 0][c], v_madd_hi_o[src + 1][c]);      
      }

      __m256i v_add_lo_e1[4];
      __m256i v_add_lo_o1[4];
      __m256i v_add_hi_e1[4];
      __m256i v_add_hi_o1[4];
      for (int dst = 0, src = 0; dst < 4; ++dst, src += 2) {
        v_add_lo_e1[dst] = _mm256_add_epi32(v_add_lo_e0[src + 0], v_add_lo_e0[src + 1]);
        v_add_lo_o1[dst] = _mm256_add_epi32(v_add_lo_o0[src + 0], v_add_lo_o0[src + 1]);
        v_add_hi_e1[dst] = _mm256_add_epi32(v_add_hi_e0[src + 0], v_add_hi_e0[src + 1]);
        v_add_hi_o1[dst] = _mm256_add_epi32(v_add_hi_o0[src + 0], v_add_hi_o0[src + 1]);      
      }

      __m256i v_add_lo_e2[2];
      __m256i v_add_lo_o2[2];
      __m256i v_add_hi_e2[2];
      __m256i v_add_hi_o2[2];
      for (int dst = 0, src = 0; dst < 2; ++dst, src += 2) {
        v_add_lo_e2[dst] = _mm256_add_epi32(v_add_lo_e1[src + 0], v_add_lo_e1[src + 1]);
        v_add_lo_o2[dst] = _mm256_add_epi32(v_add_lo_o1[src + 0], v_add_lo_o1[src + 1]);
        v_add_hi_e2[dst] = _mm256_add_epi32(v_add_hi_e1[src + 0], v_add_hi_e1[src + 1]);
        v_add_hi_o2[dst] = _mm256_add_epi32(v_add_hi_o1[src + 0], v_add_hi_o1[src + 1]);
      }

      __m256i v_trunc_lo_e = truncate_avx2(_mm256_add_epi32(v_add_lo_e2[0], v_add_lo_e2[1]), debias, shift_2nd);
      __m256i v_trunc_lo_o = truncate_avx2(_mm256_add_epi32(v_add_lo_o2[0], v_add_lo_o2[1]), debias, shift_2nd);
      __m256i v_trunc_hi_e = truncate_avx2(_mm256_add_epi32(v_add_hi_e2[0], v_add_hi_e2[1]), debias, shift_2nd);
      __m256i v_trunc_hi_o = truncate_avx2(_mm256_add_epi32(v_add_hi_o2[0], v_add_hi_o2[1]), debias, shift_2nd);

      __m256i v_result_e = _mm256_packs_epi32(v_trunc_lo_e, v_trunc_hi_e);
      __m256i v_result_o = _mm256_packs_epi32(v_trunc_lo_o, v_trunc_hi_o);

      v_result_e = _mm256_permute4x64_epi64(v_result_e, _MM_SHUFFLE(3, 1, 2, 0));
      v_result_o = _mm256_permute4x64_epi64(v_result_o, _MM_SHUFFLE(3, 1, 2, 0));

      _mm256_store_si256((__m256i*)dst, v_result_e);
      dst += 16;
      _mm256_store_si256((__m256i*)dst, v_result_o);
      dst += 16;
    }
  }
#undef NUM_PARTS
#undef PART_DIMENSION
#endif

}


static void fast_inverse_tr_32x32_avx2_ver(const int16_t* src, __m256i* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const int32_t* c_ptr = (const int32_t*)coeff; // Handle as 32 bit integer to load two coeffs into vector at the same time
  const __m256i* v_src_raw = (const __m256i*)src;

  __m256i v_src[16][4];
  for (int d = 0, s = 0; d < 16; ++d, s += 4) {
    v_src[d][0] = _mm256_unpacklo_epi16(v_src_raw[s + 0], v_src_raw[s + 2]);
    v_src[d][1] = _mm256_unpackhi_epi16(v_src_raw[s + 0], v_src_raw[s + 2]);
    v_src[d][2] = _mm256_unpacklo_epi16(v_src_raw[s + 1], v_src_raw[s + 3]);
    v_src[d][3] = _mm256_unpackhi_epi16(v_src_raw[s + 1], v_src_raw[s + 3]);
  }

  for (int row = 0, d = 0; row < 32; ++row, d += 2) {
    __m256i v_res_0 = _mm256_setzero_si256();
    __m256i v_res_1 = _mm256_setzero_si256();
    __m256i v_res_2 = _mm256_setzero_si256();
    __m256i v_res_3 = _mm256_setzero_si256();
    if(skip_line == 0) {
      for (int i = 0; i < 16; ++i) {
        const __m256i v_coeff = _mm256_set1_epi32(*c_ptr);
        __m256i v_madd_0 = _mm256_madd_epi16(v_src[i][0], v_coeff);
        __m256i v_madd_1 = _mm256_madd_epi16(v_src[i][1], v_coeff);
        __m256i v_madd_2 = _mm256_madd_epi16(v_src[i][2], v_coeff);
        __m256i v_madd_3 = _mm256_madd_epi16(v_src[i][3], v_coeff);
        v_res_0 = _mm256_add_epi32(v_res_0, v_madd_0);
        v_res_1 = _mm256_add_epi32(v_res_1, v_madd_1);
        v_res_2 = _mm256_add_epi32(v_res_2, v_madd_2);
        v_res_3 = _mm256_add_epi32(v_res_3, v_madd_3);
        c_ptr++;
      }

      __m256i v_trunc_0 = truncate_avx2(v_res_0, debias, shift);
      __m256i v_trunc_1 = truncate_avx2(v_res_1, debias, shift);
      __m256i v_trunc_2 = truncate_avx2(v_res_2, debias, shift);
      __m256i v_trunc_3 = truncate_avx2(v_res_3, debias, shift);

      dst[d + 0] = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
      dst[d + 1] = _mm256_packs_epi32(v_trunc_2, v_trunc_3);
    }
    else {
      for (int i = 0; i < 16; ++i) {
        const __m256i v_coeff = _mm256_set1_epi32(*c_ptr);
        __m256i v_madd_0 = _mm256_madd_epi16(v_src[i][0], v_coeff);
        __m256i v_madd_1 = _mm256_madd_epi16(v_src[i][1], v_coeff);
        v_res_0 = _mm256_add_epi32(v_res_0, v_madd_0);
        v_res_1 = _mm256_add_epi32(v_res_1, v_madd_1);
        c_ptr++;
      }

      __m256i v_trunc_0 = truncate_avx2(v_res_0, debias, shift);
      __m256i v_trunc_1 = truncate_avx2(v_res_1, debias, shift);

      dst[d + 0] = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
      dst[d + 1] = _mm256_setzero_si256();
    }
  }
}

static void fast_inverse_tr_32x32_avx2_hor(const __m256i* src, int16_t* dst, const int16_t* coeff, int32_t shift, int line, int skip_line, int skip_line2)
{
  const int32_t    add = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const int32_t* c_ptr = (const int32_t*)coeff; // Handle as 32 bit integer to load two coeffs into vector at the same time

  // Do a 32 bit transpose on input
  __m256i v_tmp32_lo[32];
  __m256i v_tmp32_hi[32];
  for (int d = 0, s = 0; d < 32; d += 2, s += 4) {
    v_tmp32_lo[d + 0] = _mm256_unpacklo_epi32(src[s + 0], src[s + 2]);
    v_tmp32_lo[d + 1] = _mm256_unpacklo_epi32(src[s + 1], src[s + 3]);
    v_tmp32_hi[d + 0] = _mm256_unpackhi_epi32(src[s + 0], src[s + 2]);
    v_tmp32_hi[d + 1] = _mm256_unpackhi_epi32(src[s + 1], src[s + 3]);
  }

  __m256i v_tmp64_lo[32];
  __m256i v_tmp64_hi[32];
  for (int i = 0; i < 32; i += 4) {
    v_tmp64_lo[i + 0] = _mm256_unpacklo_epi64(v_tmp32_lo[i + 0], v_tmp32_lo[i + 2]);
    v_tmp64_lo[i + 1] = _mm256_unpacklo_epi64(v_tmp32_lo[i + 1], v_tmp32_lo[i + 3]);
    v_tmp64_lo[i + 2] = _mm256_unpacklo_epi64(v_tmp32_hi[i + 0], v_tmp32_hi[i + 2]);
    v_tmp64_lo[i + 3] = _mm256_unpacklo_epi64(v_tmp32_hi[i + 1], v_tmp32_hi[i + 3]);

    v_tmp64_hi[i + 0] = _mm256_unpackhi_epi64(v_tmp32_lo[i + 0], v_tmp32_lo[i + 2]);
    v_tmp64_hi[i + 1] = _mm256_unpackhi_epi64(v_tmp32_lo[i + 1], v_tmp32_lo[i + 3]);
    v_tmp64_hi[i + 2] = _mm256_unpackhi_epi64(v_tmp32_hi[i + 0], v_tmp32_hi[i + 2]);
    v_tmp64_hi[i + 3] = _mm256_unpackhi_epi64(v_tmp32_hi[i + 1], v_tmp32_hi[i + 3]);
  }

  __m256i v_src[64];
  for (int d = 0, s = 0; d < 64; d += 16, s += 8) {
    v_src[d + 0]  = _mm256_permute2x128_si256(v_tmp64_lo[s + 0], v_tmp64_lo[s + 4], 0x20);
    v_src[d + 1]  = _mm256_permute2x128_si256(v_tmp64_hi[s + 0], v_tmp64_hi[s + 4], 0x20);
    v_src[d + 2]  = _mm256_permute2x128_si256(v_tmp64_lo[s + 2], v_tmp64_lo[s + 6], 0x20);
    v_src[d + 3]  = _mm256_permute2x128_si256(v_tmp64_hi[s + 2], v_tmp64_hi[s + 6], 0x20);

    v_src[d + 4]  = _mm256_permute2x128_si256(v_tmp64_lo[s + 0], v_tmp64_lo[s + 4], 0x31);
    v_src[d + 5]  = _mm256_permute2x128_si256(v_tmp64_hi[s + 0], v_tmp64_hi[s + 4], 0x31);
    v_src[d + 6]  = _mm256_permute2x128_si256(v_tmp64_lo[s + 2], v_tmp64_lo[s + 6], 0x31);
    v_src[d + 7]  = _mm256_permute2x128_si256(v_tmp64_hi[s + 2], v_tmp64_hi[s + 6], 0x31);

    v_src[d + 8]  = _mm256_permute2x128_si256(v_tmp64_lo[s + 1], v_tmp64_lo[s + 5], 0x20);
    v_src[d + 9]  = _mm256_permute2x128_si256(v_tmp64_hi[s + 1], v_tmp64_hi[s + 5], 0x20);
    v_src[d + 10] = _mm256_permute2x128_si256(v_tmp64_lo[s + 3], v_tmp64_lo[s + 7], 0x20);
    v_src[d + 11] = _mm256_permute2x128_si256(v_tmp64_hi[s + 3], v_tmp64_hi[s + 7], 0x20);

    v_src[d + 12] = _mm256_permute2x128_si256(v_tmp64_lo[s + 1], v_tmp64_lo[s + 5], 0x31);
    v_src[d + 13] = _mm256_permute2x128_si256(v_tmp64_hi[s + 1], v_tmp64_hi[s + 5], 0x31);
    v_src[d + 14] = _mm256_permute2x128_si256(v_tmp64_lo[s + 3], v_tmp64_lo[s + 7], 0x31);
    v_src[d + 15] = _mm256_permute2x128_si256(v_tmp64_hi[s + 3], v_tmp64_hi[s + 7], 0x31);
  }

  __m256i v_tmp[64];
  for (int row = 0, d = 0; row < 32; ++row, d += 2) {
    __m256i v_res_0 = _mm256_setzero_si256();
    __m256i v_res_1 = _mm256_setzero_si256();
    __m256i v_res_2 = _mm256_setzero_si256();
    __m256i v_res_3 = _mm256_setzero_si256();
    for (int i = 0; i < 16; ++i) {
      const __m256i v_coeff = _mm256_set1_epi32(*c_ptr);
      __m256i v_madd_0 = _mm256_madd_epi16(v_src[i + 0], v_coeff);
      __m256i v_madd_1 = _mm256_madd_epi16(v_src[i + 16], v_coeff);
      __m256i v_madd_2 = _mm256_madd_epi16(v_src[i + 32], v_coeff);
      __m256i v_madd_3 = _mm256_madd_epi16(v_src[i + 48], v_coeff);

      v_res_0 = _mm256_add_epi32(v_madd_0, v_res_0);
      v_res_1 = _mm256_add_epi32(v_madd_1, v_res_1);
      v_res_2 = _mm256_add_epi32(v_madd_2, v_res_2);
      v_res_3 = _mm256_add_epi32(v_madd_3, v_res_3);
      c_ptr++;
    }

    __m256i v_trunc_0 = truncate_avx2(v_res_0, debias, shift);
    __m256i v_trunc_1 = truncate_avx2(v_res_1, debias, shift);
    __m256i v_trunc_2 = truncate_avx2(v_res_2, debias, shift);
    __m256i v_trunc_3 = truncate_avx2(v_res_3, debias, shift);

    v_tmp[d + 0] = _mm256_packs_epi32(v_trunc_0, v_trunc_1);
    v_tmp[d + 1] = _mm256_packs_epi32(v_trunc_2, v_trunc_3);
  }

  for (int i = 0; i < 64; ++i) {
    v_tmp[i] = _mm256_permute4x64_epi64(v_tmp[i], _MM_SHUFFLE(3, 1, 2, 0));
  }

  __m256i v_result[64];
  transpose_avx2(v_tmp, v_result, 32, 32);

  for (int i = 0; i < 64; ++i) {
    _mm256_store_si256((__m256i*)dst, v_result[i]);
    dst += 16;
  }
}

static void fast_inverse_tr_32x32_avx2(const int16_t* src, int16_t* dst, tr_type_t hor, tr_type_t ver)
{
  const int width = 32;
  const int height = 32;

  int skip_width = (hor != DCT2 && width == 32) ? 16 : width > 32 ? width - 32 : 0;
  int skip_height = (ver != DCT2 && height == 32) ? 16 : (height > 32 ? height - 32 : 0);

  const int32_t shift_1st = INVERSE_SHIFT_1ST;
  const int32_t shift_2nd = INVERSE_SHIFT_2ND;

  const int16_t* ver_coeff = &uvg_g_dct_32_t[0][0];
  const int16_t* hor_coeff = fi_dct2_32xN_coeff_hor;
  if (hor == DST7) {
    hor_coeff = fi_dst7_32xN_coeff_hor;
  } else if (hor == DCT8) {
    hor_coeff = fi_dct8_32xN_coeff_hor;
  }
  if (ver == DST7) {
    ver_coeff = &uvg_g_dst7_32_t[0][0];
  } else if (ver == DCT8) {
    ver_coeff = &uvg_g_dct8_32[0][0];
  }

  __m256i v_ver_pass_out[64];
  fast_inverse_tr_32x32_avx2_ver(src, v_ver_pass_out, ver_coeff, shift_1st, width, skip_width, skip_height);

  fast_inverse_tr_32x8_avx2_hor(v_ver_pass_out, dst, hor_coeff, shift_2nd, height, 0, skip_width);
}


static dct_full_pass* dct_function_table[6][6] = {
  { NULL,                      NULL,                      fast_forward_tr_2x8_avx2,  fast_forward_tr_2x16_avx2,  fast_forward_tr_2x32_avx2,  NULL },
  { NULL,                      fast_forward_tr_4x4_avx2,  fast_forward_tr_4x8_avx2,  fast_forward_tr_4x16_avx2,  fast_forward_tr_4x32_avx2,  NULL },
  { fast_forward_tr_8x2_avx2,  fast_forward_tr_8x4_avx2,  fast_forward_tr_8x8_avx2,  fast_forward_tr_8x16_avx2,  fast_forward_tr_8x32_avx2,  NULL },
  { fast_forward_tr_16x2_avx2, fast_forward_tr_16x4_avx2, fast_forward_tr_16x8_avx2, fast_forward_tr_16x16_avx2, fast_forward_tr_16x32_avx2, NULL },
  { fast_forward_tr_32x2_avx2, fast_forward_tr_32x4_avx2, fast_forward_tr_32x8_avx2, fast_forward_tr_32x16_avx2, fast_forward_tr_32x32_avx2, NULL },
  { NULL,                      NULL,                      NULL,                      NULL,                       NULL,                       NULL }
};


static dct_full_pass* idct_function_table[6][6] = {
  { NULL,                      NULL,                      fast_inverse_tr_2x8_avx2,  fast_inverse_tr_2x16_avx2,  fast_inverse_tr_2x32_avx2,  NULL },
  { NULL,                      fast_inverse_tr_4x4_avx2,  fast_inverse_tr_4x8_avx2,  fast_inverse_tr_4x16_avx2,  fast_inverse_tr_4x32_avx2,  NULL },
  { fast_inverse_tr_8x2_avx2,  fast_inverse_tr_8x4_avx2,  fast_inverse_tr_8x8_avx2,  fast_inverse_tr_8x16_avx2,  fast_inverse_tr_8x32_avx2,  NULL },
  { fast_inverse_tr_16x2_avx2, fast_inverse_tr_16x4_avx2, fast_inverse_tr_16x8_avx2, fast_inverse_tr_16x16_avx2, fast_inverse_tr_16x32_avx2, NULL },
  { fast_inverse_tr_32x2_avx2, fast_inverse_tr_32x4_avx2, fast_inverse_tr_32x8_avx2, fast_inverse_tr_32x16_avx2, fast_inverse_tr_32x32_avx2, NULL },
  { NULL,                      NULL,                      NULL,                      NULL,                       NULL,                       NULL },
};


extern void uvg_get_tr_type(
  int8_t width,
  int8_t height,
  color_t color,
  const cu_info_t* tu,
  tr_type_t* hor_out,
  tr_type_t* ver_out,
  const int8_t mts_type);

static void mts_dct_avx2(
  const int8_t bitdepth,
  const color_t color,
  const cu_info_t* tu,
  const int8_t width,
  const int8_t height,
  const int16_t* input,
  int16_t* output,
  const int8_t mts_type)
{
  tr_type_t type_hor;
  tr_type_t type_ver;

  uvg_get_tr_type(width, height, color, tu, &type_hor, &type_ver, mts_type);

  if (type_hor == DCT2 && type_ver == DCT2 && !tu->lfnst_idx && width == height)
  {
    dct_func* dct_func = uvg_get_dct_func(width, height, color, tu->type);
    dct_func(bitdepth, input, output);
  }
  else{
    const int log2_width_minus1  = uvg_g_convert_to_log2[width] - 1;
    const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
    // Transforms with 1 lenght dimensions are handled separately since their interface differ from other full pass functions
    if (height == 1) {
      if (width == 16) {
        fast_forward_DCT2_B16_avx2_hor(input, (__m256i*)output, type_hor == DCT2 ? ff_dct2_16xN_coeff_hor : ff_dst7_16xN_coeff_hor, 3, 1, 0, 0);
      } else if (width == 32) {
        fast_forward_DCT2_B32_avx2_hor(input, (__m256i*)output, ff_dct2_32xN_coeff_hor, 4, 1, 0, 0);        
      }
    }
    else if (width == 1){
      if (height == 16) {
        fast_forward_DCT2_B16_avx2_hor(input, (__m256i*)output, type_ver == DCT2 ? ff_dct2_16xN_coeff_hor : ff_dst7_16xN_coeff_hor, 3, 1, 0, 0);
      } else if (height == 32) {
        fast_forward_DCT2_B32_avx2_hor(input, (__m256i*)output, ff_dct2_32xN_coeff_hor, 4, 1, 0, 0);        
      }
    }
    else {
      dct_full_pass* dct_func = dct_function_table[log2_width_minus1][log2_height_minus1];
      dct_func(input, output, type_hor, type_ver);
    }
  }
}


static void mts_idct_avx2(
  const int8_t bitdepth,
  const color_t color,
  const cu_info_t* tu,
  const int8_t width,
  const int8_t height,
  const int16_t* input,
  int16_t* output,
  const int8_t mts_type)
{
  tr_type_t type_hor;
  tr_type_t type_ver;

  uvg_get_tr_type(width, height, color, tu, &type_hor, &type_ver, mts_type);

  if (type_hor == DCT2 && type_ver == DCT2 && width == height)
  {
    dct_func* idct_func = uvg_get_idct_func(width, height, color, tu->type);
    idct_func(bitdepth, input, output);
  }
  else {
    const int log2_width_minus1  = uvg_g_convert_to_log2[width] - 1;
    const int log2_height_minus1 = uvg_g_convert_to_log2[height] - 1;
    // Transforms with 1 lenght dimensions can be transformed with existing forward functions
    if (height == 1) {
      if (width == 16) {
        fast_forward_DCT2_B16_avx2_hor(input, (__m256i*)output, type_hor == DCT2 ? fi_dct2_16x1_coeff_hor : fi_dst7_16x1_coeff_hor, 13, 1, 0, 0);
        _mm256_store_si256((__m256i*)output, _mm256_permute4x64_epi64(_mm256_load_si256((__m256i*)output), _MM_SHUFFLE(3, 1, 2, 0)));
      } else if (width == 32) {
        fast_forward_DCT2_B32_avx2_hor(input, (__m256i*)output, fi_dct2_32xN_coeff_hor, 13, 1, 0, 0);        
      }
    }
    else if (width == 1){
      if (height == 16) {
        fast_forward_DCT2_B16_avx2_hor(input, (__m256i*)output, type_ver == DCT2 ? fi_dct2_16x1_coeff_hor : fi_dst7_16x1_coeff_hor, 13, 1, 0, 0);
        _mm256_store_si256((__m256i*)output, _mm256_permute4x64_epi64(_mm256_load_si256((__m256i*)output), _MM_SHUFFLE(3, 1, 2, 0)));
      } else if (height == 32) {
        fast_forward_DCT2_B32_avx2_hor(input, (__m256i*)output, fi_dct2_32xN_coeff_hor, 13, 1, 0, 0);        
      }
    }
    else {
      dct_full_pass* idct_func = idct_function_table[log2_width_minus1][log2_height_minus1];
      idct_func(input, output, type_hor, type_ver);
    }
  }
}

#endif // UVG_BIT_DEPTH == 8
#endif //COMPILE_INTEL_AVX2

int uvg_strategy_register_dct_avx2(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if COMPILE_INTEL_AVX2
#if UVG_BIT_DEPTH == 8
  if (bitdepth == 8){
    //success &= uvg_strategyselector_register(opaque, "fast_forward_dst_4x4", "avx2", 40, &matrix_dst_4x4_avx2);

    success &= uvg_strategyselector_register(opaque, "dct_4x4", "avx2", 40, &matrix_dct_4x4_avx2);
    success &= uvg_strategyselector_register(opaque, "dct_8x8", "avx2", 40, &matrix_dct_8x8_avx2);
    success &= uvg_strategyselector_register(opaque, "dct_16x16", "avx2", 40, &matrix_dct_16x16_avx2);
    success &= uvg_strategyselector_register(opaque, "dct_32x32", "avx2", 40, &matrix_dct_32x32_avx2);

    // success &= uvg_strategyselector_register(opaque, "fast_inverse_dst_4x4", "avx2", 40, &matrix_idst_4x4_avx2);

    success &= uvg_strategyselector_register(opaque, "idct_4x4", "avx2", 40, &matrix_idct_4x4_avx2);
    success &= uvg_strategyselector_register(opaque, "idct_8x8", "avx2", 40, &matrix_idct_8x8_avx2);
    success &= uvg_strategyselector_register(opaque, "idct_16x16", "avx2", 40, &matrix_idct_16x16_avx2);
    success &= uvg_strategyselector_register(opaque, "idct_32x32", "avx2", 40, &matrix_idct_32x32_avx2);

    success &= uvg_strategyselector_register(opaque, "mts_dct", "avx2", 40, &mts_dct_avx2);
    success &= uvg_strategyselector_register(opaque, "mts_idct", "avx2", 40, &mts_idct_avx2);

  }
#endif // UVG_BIT_DEPTH == 8
#endif //COMPILE_INTEL_AVX2  
  return success;
}
