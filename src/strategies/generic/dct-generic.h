#ifndef STRATEGIES_DCT_GENERIC_H_
#define STRATEGIES_DCT_GENERIC_H_
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

/**
 * \ingroup Optimization
 * \file
 * Generic C implementations of optimized functions.
 */

#include  "global.h" // IWYU pragma: keep
#include "cu.h"

ALIGNED(16) extern const int16_t uvg_g_dct2_2x2[2][2];
ALIGNED(16) extern const int16_t uvg_g_dct2_2x2_t[2][2];

ALIGNED(32) extern const int16_t uvg_g_dct2_4x4[4][4];
ALIGNED(32) extern const int16_t uvg_g_dct2_4x4_t[4][4];

ALIGNED(64) extern const int16_t uvg_g_dct2_8x8[8][8];
ALIGNED(64) extern const int16_t uvg_g_dct2_8x8_t[8][8];

ALIGNED(64) extern const int16_t uvg_g_dct2_16x16[16][16];
ALIGNED(64) extern const int16_t uvg_g_dct2_16x16_t[16][16];

ALIGNED(64) extern const int16_t uvg_g_dct2_32x32[32][32];
ALIGNED(64) extern const int16_t uvg_g_dct2_32x32_t[32][32];

ALIGNED(64) extern const int16_t uvg_g_dct2_64x64[64][64];
ALIGNED(64) extern const int16_t uvg_g_dct2_64x64_t[64][64];

ALIGNED(32) extern const int16_t uvg_g_dct8_4x4[4][4];
ALIGNED(32) extern const int16_t uvg_g_dct8_4x4_t[4][4];

ALIGNED(64) extern const int16_t uvg_g_dct8_8x8[8][8];
ALIGNED(64) extern const int16_t uvg_g_dct8_8x8_t[8][8];

ALIGNED(64) extern const int16_t uvg_g_dct8_16x16[16][16];
ALIGNED(64) extern const int16_t uvg_g_dct8_16x16_t[16][16];

ALIGNED(64) extern const int16_t uvg_g_dct8_32x32[32][32];
ALIGNED(64) extern const int16_t uvg_g_dct8_32x32_t[32][32];

ALIGNED(32) extern const int16_t uvg_g_dst7_4x4[4][4];
ALIGNED(32) extern const int16_t uvg_g_dst7_4x4_t[4][4];

ALIGNED(64) extern const int16_t uvg_g_dst7_8x8[8][8];
ALIGNED(64) extern const int16_t uvg_g_dst7_8x8_t[8][8];

ALIGNED(64) extern const int16_t uvg_g_dst7_16x16[16][16];
ALIGNED(64) extern const int16_t uvg_g_dst7_16x16_t[16][16];

ALIGNED(64) extern const int16_t uvg_g_dst7_32x32[32][32];
ALIGNED(64) extern const int16_t uvg_g_dst7_32x32_t[32][32];

void uvg_get_tr_type(
  int8_t width,
  int8_t height,
  color_t color,
  const cu_info_t* tu,
  tr_type_t* hor_out,
  tr_type_t* ver_out,
  int8_t mts_type
);

void uvg_mts_dct_generic(
  const int8_t bitdepth,
  const color_t color,
  const cu_info_t* tu,
  const int8_t width,
  const int8_t height,
  const int16_t* input,
  int16_t* output,
  const int8_t mts_type);

void uvg_mts_idct_generic(
  const int8_t bitdepth,
  const color_t color,
  const cu_info_t* tu,
  const int8_t width,
  const int8_t height,
  const int16_t* input,
  int16_t* output,
  const int8_t mts_type);

int uvg_strategy_register_dct_generic(void* opaque, uint8_t bitdepth);

// Deprecated.
void uvg_dct_4x4_generic(int8_t bitdepth, const int16_t* input, int16_t* output);
void uvg_dct_8x8_generic(int8_t bitdepth, const int16_t* input, int16_t* output);
void uvg_dct_16x16_generic(int8_t bitdepth, const int16_t* input, int16_t* output);
void uvg_dct_32x32_generic(int8_t bitdepth, const int16_t* input, int16_t* output);
void uvg_idct_4x4_generic(int8_t bitdepth, const int16_t* input, int16_t* output);
void uvg_idct_8x8_generic(int8_t bitdepth, const int16_t* input, int16_t* output);
void uvg_idct_16x16_generic(int8_t bitdepth, const int16_t* input, int16_t* output);
void uvg_idct_32x32_generic(int8_t bitdepth, const int16_t* input, int16_t* output);

#endif //STRATEGIES_DCT_GENERIC_H_
