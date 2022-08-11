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
#define DEFINE_DST7_P4_MATRIX(a,b,c,d) { \
    { a,  b,  c,  d},\
    { c,  c,  0, -c},\
    { d, -a, -c,  b},\
    { b, -d,  c, -a},\
}

#define DEFINE_DST7_P4_MATRIX_T(a,b,c,d) { \
    { a,  c,  d,  b},\
    { b,  c, -a, -d},\
    { c,  0, -c,  c},\
    { d, -c,  b, -a},\
}

#define DEFINE_DST7_P8_MATRIX(a,b,c,d,e,f,g,h) \
{\
   { a,  b,  c,  d,  e,  f,  g,  h},\
   { c,  f,  h,  e,  b, -a, -d, -g},\
   { e,  g,  b, -c, -h, -d,  a,  f},\
   { g,  c, -d, -f,  a,  h,  b, -e},\
   { h, -a, -g,  b,  f, -c, -e,  d},\
   { f, -e, -a,  g, -d, -b,  h, -c},\
   { d, -h,  e, -a, -c,  g, -f,  b},\
   { b, -d,  f, -h,  g, -e,  c, -a},\
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

#define DEFINE_DST7_P16_MATRIX(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) \
{ \
   { a,  b,  c,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p}, \
   { c,  f,  i,  l,  o,  o,  l,  i,  f,  c,  0, -c, -f, -i, -l, -o}, \
   { e,  j,  o,  m,  h,  c, -b, -g, -l, -p, -k, -f, -a,  d,  i,  n}, \
   { g,  n,  l,  e, -b, -i, -p, -j, -c,  d,  k,  o,  h,  a, -f, -m}, \
   { i,  o,  f, -c, -l, -l, -c,  f,  o,  i,  0, -i, -o, -f,  c,  l}, \
   { k,  k,  0, -k, -k,  0,  k,  k,  0, -k, -k,  0,  k,  k,  0, -k}, \
   { m,  g, -f, -n, -a,  l,  h, -e, -o, -b,  k,  i, -d, -p, -c,  j}, \
   { o,  c, -l, -f,  i,  i, -f, -l,  c,  o,  0, -o, -c,  l,  f, -i}, \
   { p, -a, -o,  b,  n, -c, -m,  d,  l, -e, -k,  f,  j, -g, -i,  h}, \
   { n, -e, -i,  j,  d, -o,  a,  m, -f, -h,  k,  c, -p,  b,  l, -g}, \
   { l, -i, -c,  o, -f, -f,  o, -c, -i,  l,  0, -l,  i,  c, -o,  f}, \
   { j, -m,  c,  g, -p,  f,  d, -n,  i,  a, -k,  l, -b, -h,  o, -e}, \
   { h, -p,  i, -a, -g,  o, -j,  b,  f, -n,  k, -c, -e,  m, -l,  d}, \
   { f, -l,  o, -i,  c,  c, -i,  o, -l,  f,  0, -f,  l, -o,  i, -c}, \
   { d, -h,  l, -p,  m, -i,  e, -a, -c,  g, -k,  o, -n,  j, -f,  b}, \
   { b, -d,  f, -h,  j, -l,  n, -p,  o, -m,  k, -i,  g, -e,  c, -a}, \
}

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



#define DEFINE_DST7_P32_MATRIX(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,A,B,C,D,E,F) \
{ \
    {a,  b,  c,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  t,  u,  v,  w,  x,  y,  z,  A,  B,  C,  D,  E,  F}, \
    {c,  f,  i,  l,  o,  r,  u,  x,  A,  D,  F,  C,  z,  w,  t,  q,  n,  k,  h,  e,  b, -a, -d, -g, -j, -m, -p, -s, -v, -y, -B, -E}, \
    {e,  j,  o,  t,  y,  D,  D,  y,  t,  o,  j,  e,  0, -e, -j, -o, -t, -y, -D, -D, -y, -t, -o, -j, -e,  0,  e,  j,  o,  t,  y,  D}, \
    {g,  n,  u,  B,  D,  w,  p,  i,  b, -e, -l, -s, -z, -F, -y, -r, -k, -d,  c,  j,  q,  x,  E,  A,  t,  m,  f, -a, -h, -o, -v, -C}, \
    {i,  r,  A,  C,  t,  k,  b, -g, -p, -y, -E, -v, -m, -d,  e,  n,  w,  F,  x,  o,  f, -c, -l, -u, -D, -z, -q, -h,  a,  j,  s,  B}, \
    {k,  v,  F,  u,  j, -a, -l, -w, -E, -t, -i,  b,  m,  x,  D,  s,  h, -c, -n, -y, -C, -r, -g,  d,  o,  z,  B,  q,  f, -e, -p, -A}, \
    {m,  z,  z,  m,  0, -m, -z, -z, -m,  0,  m,  z,  z,  m,  0, -m, -z, -z, -m,  0,  m,  z,  z,  m,  0, -m, -z, -z, -m,  0,  m,  z}, \
    {o,  D,  t,  e, -j, -y, -y, -j,  e,  t,  D,  o,  0, -o, -D, -t, -e,  j,  y,  y,  j, -e, -t, -D, -o,  0,  o,  D,  t,  e, -j, -y}, \
    {q,  E,  n, -c, -t, -B, -k,  f,  w,  y,  h, -i, -z, -v, -e,  l,  C,  s,  b, -o, -F, -p,  a,  r,  D,  m, -d, -u, -A, -j,  g,  x}, \
    {s,  A,  h, -k, -D, -p,  c,  v,  x,  e, -n, -F, -m,  f,  y,  u,  b, -q, -C, -j,  i,  B,  r, -a, -t, -z, -g,  l,  E,  o, -d, -w}, \
    {u,  w,  b, -s, -y, -d,  q,  A,  f, -o, -C, -h,  m,  E,  j, -k, -F, -l,  i,  D,  n, -g, -B, -p,  e,  z,  r, -c, -x, -t,  a,  v}, \
    {w,  s, -d, -A, -o,  h,  E,  k, -l, -D, -g,  p,  z,  c, -t, -v,  a,  x,  r, -e, -B, -n,  i,  F,  j, -m, -C, -f,  q,  y,  b, -u}, \
    {y,  o, -j, -D, -e,  t,  t, -e, -D, -j,  o,  y,  0, -y, -o,  j,  D,  e, -t, -t,  e,  D,  j, -o, -y,  0,  y,  o, -j, -D, -e,  t}, \
    {A,  k, -p, -v,  e,  F,  f, -u, -q,  j,  B,  a, -z, -l,  o,  w, -d, -E, -g,  t,  r, -i, -C, -b,  y,  m, -n, -x,  c,  D,  h, -s}, \
    {C,  g, -v, -n,  o,  u, -h, -B,  a,  D,  f, -w, -m,  p,  t, -i, -A,  b,  E,  e, -x, -l,  q,  s, -j, -z,  c,  F,  d, -y, -k,  r}, \
    {E,  c, -B, -f,  y,  i, -v, -l,  s,  o, -p, -r,  m,  u, -j, -x,  g,  A, -d, -D,  a,  F,  b, -C, -e,  z,  h, -w, -k,  t,  n, -q}, \
    {F, -a, -E,  b,  D, -c, -C,  d,  B, -e, -A,  f,  z, -g, -y,  h,  x, -i, -w,  j,  v, -k, -u,  l,  t, -m, -s,  n,  r, -o, -q,  p}, \
    {D, -e, -y,  j,  t, -o, -o,  t,  j, -y, -e,  D,  0, -D,  e,  y, -j, -t,  o,  o, -t, -j,  y,  e, -D,  0,  D, -e, -y,  j,  t, -o}, \
    {B, -i, -s,  r,  j, -A, -a,  C, -h, -t,  q,  k, -z, -b,  D, -g, -u,  p,  l, -y, -c,  E, -f, -v,  o,  m, -x, -d,  F, -e, -w,  n}, \
    {z, -m, -m,  z,  0, -z,  m,  m, -z,  0,  z, -m, -m,  z,  0, -z,  m,  m, -z,  0,  z, -m, -m,  z,  0, -z,  m,  m, -z,  0,  z, -m}, \
    {x, -q, -g,  E, -j, -n,  A, -c, -u,  t,  d, -B,  m,  k, -D,  f,  r, -w, -a,  y, -p, -h,  F, -i, -o,  z, -b, -v,  s,  e, -C,  l}, \
    {v, -u, -a,  w, -t, -b,  x, -s, -c,  y, -r, -d,  z, -q, -e,  A, -p, -f,  B, -o, -g,  C, -n, -h,  D, -m, -i,  E, -l, -j,  F, -k}, \
    {t, -y,  e,  o, -D,  j,  j, -D,  o,  e, -y,  t,  0, -t,  y, -e, -o,  D, -j, -j,  D, -o, -e,  y, -t,  0,  t, -y,  e,  o, -D,  j}, \
    {r, -C,  k,  g, -y,  v, -d, -n,  F, -o, -c,  u, -z,  h,  j, -B,  s, -a, -q,  D, -l, -f,  x, -w,  e,  m, -E,  p,  b, -t,  A, -i}, \
    {p, -F,  q, -a, -o,  E, -r,  b,  n, -D,  s, -c, -m,  C, -t,  d,  l, -B,  u, -e, -k,  A, -v,  f,  j, -z,  w, -g, -i,  y, -x,  h}, \
    {n, -B,  w, -i, -e,  s, -F,  r, -d, -j,  x, -A,  m,  a, -o,  C, -v,  h,  f, -t,  E, -q,  c,  k, -y,  z, -l, -b,  p, -D,  u, -g}, \
    {l, -x,  C, -q,  e,  g, -s,  E, -v,  j,  b, -n,  z, -A,  o, -c, -i,  u, -F,  t, -h, -d,  p, -B,  y, -m,  a,  k, -w,  D, -r,  f}, \
    {j, -t,  D, -y,  o, -e, -e,  o, -y,  D, -t,  j,  0, -j,  t, -D,  y, -o,  e,  e, -o,  y, -D,  t, -j,  0,  j, -t,  D, -y,  o, -e}, \
    {h, -p,  x, -F,  y, -q,  i, -a, -g,  o, -w,  E, -z,  r, -j,  b,  f, -n,  v, -D,  A, -s,  k, -c, -e,  m, -u,  C, -B,  t, -l,  d}, \
    {f, -l,  r, -x,  D, -C,  w, -q,  k, -e, -a,  g, -m,  s, -y,  E, -B,  v, -p,  j, -d, -b,  h, -n,  t, -z,  F, -A,  u, -o,  i, -c}, \
    {d, -h,  l, -p,  t, -x,  B, -F,  C, -y,  u, -q,  m, -i,  e, -a, -c,  g, -k,  o, -s,  w, -A,  E, -D,  z, -v,  r, -n,  j, -f,  b}, \
    {b, -d,  f, -h,  j, -l,  n, -p,  r, -t,  v, -x,  z, -B,  D, -F,  E, -C,  A, -y,  w, -u,  s, -q,  o, -m,  k, -i,  g, -e,  c, -a}, \
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

// DCT-8
#define DEFINE_DCT8_P4_MATRIX(a,b,c,d) \
{ \
    {a,  b,  c,  d}, \
    {b,  0, -b, -b}, \
    {c, -b, -d,  a}, \
    {d, -b,  a, -c}, \
}

#define DEFINE_DCT8_P8_MATRIX(a,b,c,d,e,f,g,h) \
{ \
    {a,  b,  c,  d,  e,  f,  g,  h}, \
    {b,  e,  h, -g, -d, -a, -c, -f}, \
    {c,  h, -e, -a, -f,  g,  b,  d}, \
    {d, -g, -a, -h,  c,  e, -f, -b}, \
    {e, -d, -f,  c,  g, -b, -h,  a}, \
    {f, -a,  g,  e, -b,  h,  d, -c}, \
    {g, -c,  b, -f, -h,  d, -a,  e}, \
    {h, -f,  d, -b,  a, -c,  e, -g}, \
}

#define DEFINE_DCT8_P16_MATRIX(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) \
{ \
    {a,  b,  c,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p}, \
    {b,  e,  h,  k,  n,  0, -n, -k, -h, -e, -b, -b, -e, -h, -k, -n}, \
    {c,  h,  m, -p, -k, -f, -a, -e, -j, -o,  n,  i,  d,  b,  g,  l}, \
    {d,  k, -p, -i, -b, -f, -m,  n,  g,  a,  h,  o, -l, -e, -c, -j}, \
    {e,  n, -k, -b, -h,  0,  h,  b,  k, -n, -e, -e, -n,  k,  b,  h}, \
    {f,  0, -f, -f,  0,  f,  f,  0, -f, -f,  0,  f,  f,  0, -f, -f}, \
    {g, -n, -a, -m,  h,  f, -o, -b, -l,  i,  e, -p, -c, -k,  j,  d}, \
    {h, -k, -e,  n,  b,  0, -b, -n,  e,  k, -h, -h,  k,  e, -n, -b}, \
    {i, -h, -j,  g,  k, -f, -l,  e,  m, -d, -n,  c,  o, -b, -p,  a}, \
    {j, -e, -o,  a, -n, -f,  i,  k, -d, -p,  b, -m, -g,  h,  l, -c}, \
    {k, -b,  n,  h, -e,  0,  e, -h, -n,  b, -k, -k,  b, -n, -h,  e}, \
    {l, -b,  i,  o, -e,  f, -p, -h,  c, -m, -k,  a, -j, -n,  d, -g}, \
    {m, -e,  d, -l, -n,  f, -c,  k,  o, -g,  b, -j, -p,  h, -a,  i}, \
    {n, -h,  b, -e,  k,  0, -k,  e, -b,  h, -n, -n,  h, -b,  e, -k}, \
    {o, -k,  g, -c,  b, -f,  j, -n, -p,  l, -h,  d, -a,  e, -i,  m}, \
    {p, -n,  l, -j,  h, -f,  d, -b,  a, -c,  e, -g,  i, -k,  m, -o}, \
}


#define DEFINE_DCT8_P32_MATRIX(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,A,B,C,D,E,F) \
{ \
    {a,  b,  c,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  t,  u,  v,  w,  x,  y,  z,  A,  B,  C,  D,  E,  F}, \
    {b,  e,  h,  k,  n,  q,  t,  w,  z,  C,  F, -E, -B, -y, -v, -s, -p, -m, -j, -g, -d, -a, -c, -f, -i, -l, -o, -r, -u, -x, -A, -D}, \
    {c,  h,  m,  r,  w,  B,  0, -B, -w, -r, -m, -h, -c, -c, -h, -m, -r, -w, -B,  0,  B,  w,  r,  m,  h,  c,  c,  h,  m,  r,  w,  B}, \
    {d,  k,  r,  y,  F, -A, -t, -m, -f, -b, -i, -p, -w, -D,  C,  v,  o,  h,  a,  g,  n,  u,  B, -E, -x, -q, -j, -c, -e, -l, -s, -z}, \
    {e,  n,  w,  F, -y, -p, -g, -c, -l, -u, -D,  A,  r,  i,  a,  j,  s,  B, -C, -t, -k, -b, -h, -q, -z,  E,  v,  m,  d,  f,  o,  x}, \
    {f,  q,  B, -A, -p, -e, -g, -r, -C,  z,  o,  d,  h,  s,  D, -y, -n, -c, -i, -t, -E,  x,  m,  b,  j,  u,  F, -w, -l, -a, -k, -v}, \
    {g,  t,  0, -t, -g, -g, -t,  0,  t,  g,  g,  t,  0, -t, -g, -g, -t,  0,  t,  g,  g,  t,  0, -t, -g, -g, -t,  0,  t,  g,  g,  t}, \
    {h,  w, -B, -m, -c, -r,  0,  r,  c,  m,  B, -w, -h, -h, -w,  B,  m,  c,  r,  0, -r, -c, -m, -B,  w,  h,  h,  w, -B, -m, -c, -r}, \
    {i,  z, -w, -f, -l, -C,  t,  c,  o,  F, -q, -a, -r,  E,  n,  d,  u, -B, -k, -g, -x,  y,  h,  j,  A, -v, -e, -m, -D,  s,  b,  p}, \
    {j,  C, -r, -b, -u,  z,  g,  m,  F, -o, -e, -x,  w,  d,  p, -E, -l, -h, -A,  t,  a,  s, -B, -i, -k, -D,  q,  c,  v, -y, -f, -n}, \
    {k,  F, -m, -i, -D,  o,  g,  B, -q, -e, -z,  s,  c,  x, -u, -a, -v,  w,  b,  t, -y, -d, -r,  A,  f,  p, -C, -h, -n,  E,  j,  l}, \
    {l, -E, -h, -p,  A,  d,  t, -w, -a, -x,  s,  e,  B, -o, -i, -F,  k,  m, -D, -g, -q,  z,  c,  u, -v, -b, -y,  r,  f,  C, -n, -j}, \
    {m, -B, -c, -w,  r,  h,  0, -h, -r,  w,  c,  B, -m, -m,  B,  c,  w, -r, -h,  0,  h,  r, -w, -c, -B,  m,  m, -B, -c, -w,  r,  h}, \
    {n, -y, -c, -D,  i,  s, -t, -h,  E,  d,  x, -o, -m,  z,  b,  C, -j, -r,  u,  g, -F, -e, -w,  p,  l, -A, -a, -B,  k,  q, -v, -f}, \
    {o, -v, -h,  C,  a,  D, -g, -w,  n,  p, -u, -i,  B,  b,  E, -f, -x,  m,  q, -t, -j,  A,  c,  F, -e, -y,  l,  r, -s, -k,  z,  d}, \
    {p, -s, -m,  v,  j, -y, -g,  B,  d, -E, -a, -F,  c,  C, -f, -z,  i,  w, -l, -t,  o,  q, -r, -n,  u,  k, -x, -h,  A,  e, -D, -b}, \
    {q, -p, -r,  o,  s, -n, -t,  m,  u, -l, -v,  k,  w, -j, -x,  i,  y, -h, -z,  g,  A, -f, -B,  e,  C, -d, -D,  c,  E, -b, -F,  a}, \
    {r, -m, -w,  h,  B, -c,  0,  c, -B, -h,  w,  m, -r, -r,  m,  w, -h, -B,  c,  0, -c,  B,  h, -w, -m,  r,  r, -m, -w,  h,  B, -c}, \
    {s, -j, -B,  a, -C, -i,  t,  r, -k, -A,  b, -D, -h,  u,  q, -l, -z,  c, -E, -g,  v,  p, -m, -y,  d, -F, -f,  w,  o, -n, -x,  e}, \
    {t, -g,  0,  g, -t, -t,  g,  0, -g,  t,  t, -g,  0,  g, -t, -t,  g,  0, -g,  t,  t, -g,  0,  g, -t, -t,  g,  0, -g,  t,  t, -g}, \
    {u, -d,  B,  n, -k, -E,  g, -r, -x,  a, -y, -q,  h, -F, -j,  o,  A, -c,  v,  t, -e,  C,  m, -l, -D,  f, -s, -w,  b, -z, -p,  i}, \
    {v, -a,  w,  u, -b,  x,  t, -c,  y,  s, -d,  z,  r, -e,  A,  q, -f,  B,  p, -g,  C,  o, -h,  D,  n, -i,  E,  m, -j,  F,  l, -k}, \
    {w, -c,  r,  B, -h,  m,  0, -m,  h, -B, -r,  c, -w, -w,  c, -r, -B,  h, -m,  0,  m, -h,  B,  r, -c,  w,  w, -c,  r,  B, -h,  m}, \
    {x, -f,  m, -E, -q,  b, -t, -B,  j, -i,  A,  u, -c,  p,  F, -n,  e, -w, -y,  g, -l,  D,  r, -a,  s,  C, -k,  h, -z, -v,  d, -o}, \
    {y, -i,  h, -x, -z,  j, -g,  w,  A, -k,  f, -v, -B,  l, -e,  u,  C, -m,  d, -t, -D,  n, -c,  s,  E, -o,  b, -r, -F,  p, -a,  q}, \
    {z, -l,  c, -q,  E,  u, -g,  h, -v, -D,  p, -b,  m, -A, -y,  k, -d,  r, -F, -t,  f, -i,  w,  C, -o,  a, -n,  B,  x, -j,  e, -s}, \
    {A, -o,  c, -j,  v,  F, -t,  h, -e,  q, -C, -y,  m, -a,  l, -x, -D,  r, -f,  g, -s,  E,  w, -k,  b, -n,  z,  B, -p,  d, -i,  u}, \
    {B, -r,  h, -c,  m, -w,  0,  w, -m,  c, -h,  r, -B, -B,  r, -h,  c, -m,  w,  0, -w,  m, -c,  h, -r,  B,  B, -r,  h, -c,  m, -w}, \
    {C, -u,  m, -e,  d, -l,  t, -B, -D,  v, -n,  f, -c,  k, -s,  A,  E, -w,  o, -g,  b, -j,  r, -z, -F,  x, -p,  h, -a,  i, -q,  y}, \
    {D, -x,  r, -l,  f, -a,  g, -m,  s, -y,  E,  C, -w,  q, -k,  e, -b,  h, -n,  t, -z,  F,  B, -v,  p, -j,  d, -c,  i, -o,  u, -A}, \
    {E, -A,  w, -s,  o, -k,  g, -c,  b, -f,  j, -n,  r, -v,  z, -D, -F,  B, -x,  t, -p,  l, -h,  d, -a,  e, -i,  m, -q,  u, -y,  C}, \
    {F, -D,  B, -z,  x, -v,  t, -r,  p, -n,  l, -j,  h, -f,  d, -b,  a, -c,  e, -g,  i, -k,  m, -o,  q, -s,  u, -w,  y, -A,  C, -E}, \
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
  if (skip_line)
  {
    const int reduced_line = 8, cutoff = 8;
    int16_t* dst2 = output + reduced_line;
    for (int j = 0; j < cutoff; j++)
    {
      memset(dst2, 0, sizeof(int16_t) * skip_line);
      dst2 += 16;
    }
  }

  if (skip_line2)
  {
    int16_t* dst2 = output + 16 * 8;
    memset(dst2, 0, sizeof(int16_t) * 16 * skip_line2);
  }
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
  
  if (skip_line)
  {
    int16_t* dst2 = dst + reduced_line;
    for (j = 0; j < cutoff; j++)
    {
      memset(dst2, 0, sizeof(int16_t) * skip_line);
      dst2 += 32;
    }
  }

  if (skip_line2)
  {
    int16_t* dst2 = dst + 32 * cutoff;
    memset(dst2, 0, sizeof(int16_t) * 32 * skip_line2);
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


extern void uvg_get_tr_type(
  int8_t width,
  color_t color,
  const cu_info_t* tu,
  tr_type_t* hor_out,
  tr_type_t* ver_out,
  const int8_t mts_idx);

static void mts_dct_avx2(
  const int8_t bitdepth,
  const color_t color,
  const cu_info_t* tu,
  const int8_t width,
  const int8_t height,
  const int16_t* input,
  int16_t* output,
  const int8_t mts_idx)
{
  tr_type_t type_hor;
  tr_type_t type_ver;
  // ISP_TODO: height passed but not used

  uvg_get_tr_type(width, color, tu, &type_hor, &type_ver, mts_idx);

  if (type_hor == DCT2 && type_ver == DCT2 && !tu->lfnst_idx)
  {
    dct_func* dct_func = uvg_get_dct_func(width, height, color, tu->type);
    dct_func(bitdepth, input, output);
  }
  else
  {
    const int log2_width_minus2 = uvg_g_convert_to_bit[width];

    tr_func* dct = dct_table[log2_width_minus2];

    dct(input, output, type_hor, type_ver, bitdepth, tu->lfnst_idx);
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
  const int8_t mts_idx)
{
  tr_type_t type_hor;
  tr_type_t type_ver;

  uvg_get_tr_type(width, color, tu, &type_hor, &type_ver, mts_idx);

  if (type_hor == DCT2 && type_ver == DCT2)
  {
    dct_func* idct_func = uvg_get_idct_func(width, color, tu->type);
    idct_func(bitdepth, input, output);
  }
  else
  {
    const int log2_width_minus2 = uvg_g_convert_to_bit[width];

    tr_func* idct = idct_table[log2_width_minus2];

    idct(input, output, type_hor, type_ver, bitdepth, tu->lfnst_idx);
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
    success &= uvg_strategyselector_register(opaque, "fast_forward_dst_4x4", "avx2", 40, &matrix_dst_4x4_avx2);

    success &= uvg_strategyselector_register(opaque, "dct_4x4", "avx2", 40, &matrix_dct_4x4_avx2);
    success &= uvg_strategyselector_register(opaque, "dct_8x8", "avx2", 40, &matrix_dct_8x8_avx2);
    success &= uvg_strategyselector_register(opaque, "dct_16x16", "avx2", 40, &matrix_dct_16x16_avx2);
    success &= uvg_strategyselector_register(opaque, "dct_32x32", "avx2", 40, &matrix_dct_32x32_avx2);

    success &= uvg_strategyselector_register(opaque, "fast_inverse_dst_4x4", "avx2", 40, &matrix_idst_4x4_avx2);

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
