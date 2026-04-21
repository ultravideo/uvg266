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

#include "uvg266.h"
extern "C" {
#include "transform-neon.h"
#include "transpose-neon.h"
#include "strategies/generic/dct-generic.h"
#include "global.h"
}

#include "debug-neon.hpp"
#include "transform-common-neon.hpp"
#include "utility-neon.hpp"

#include <arm_neon.h>

#if !defined(NEON_FORCE_FAST_DCT2)
  #define NEON_FORCE_FAST_DCT2 0
#endif

#if !defined(NEON_PREFER_FAST_DCT2)
  #define NEON_PREFER_FAST_DCT2 1
#endif

// ---------------------------------------------------------
// Coefficient tables
// ---------------------------------------------------------

static constexpr const int16_t* transform_coeffs[7][3] = {
   // DCT-2                 // DCT-8                 // DST-7
  {nullptr,                 nullptr,                 nullptr                }, // size=1 (padding)
  {&uvg_g_dct2_2x2[0][0],   nullptr,                 nullptr                }, // size=2
  {&uvg_g_dct2_4x4[0][0],   &uvg_g_dct8_4x4[0][0],   &uvg_g_dst7_4x4[0][0]  }, // size=4
  {&uvg_g_dct2_8x8[0][0],   &uvg_g_dct8_8x8[0][0],   &uvg_g_dst7_8x8[0][0]  }, // size=8
  {&uvg_g_dct2_16x16[0][0], &uvg_g_dct8_16x16[0][0], &uvg_g_dst7_16x16[0][0]}, // size=16
  {&uvg_g_dct2_32x32[0][0], &uvg_g_dct8_32x32[0][0], &uvg_g_dst7_32x32[0][0]}, // size=32
  {nullptr,                 nullptr,                 nullptr                }, // size=64 (to be implemented)
};

static constexpr const int16_t* transform_coeffs_t[7][3] = {
   // DCT-2                   // DCT-8                   // DST-7
  {nullptr,                   nullptr,                   nullptr                  }, // size=1 (padding)
  {&uvg_g_dct2_2x2_t[0][0],   nullptr,                   nullptr                  }, // size=2
  {&uvg_g_dct2_4x4_t[0][0],   &uvg_g_dct8_4x4_t[0][0],   &uvg_g_dst7_4x4_t[0][0]  }, // size=4
  {&uvg_g_dct2_8x8_t[0][0],   &uvg_g_dct8_8x8_t[0][0],   &uvg_g_dst7_8x8_t[0][0]  }, // size=8
  {&uvg_g_dct2_16x16_t[0][0], &uvg_g_dct8_16x16_t[0][0], &uvg_g_dst7_16x16_t[0][0]}, // size=16
  {&uvg_g_dct2_32x32_t[0][0], &uvg_g_dct8_32x32_t[0][0], &uvg_g_dst7_32x32_t[0][0]}, // size=32
  {nullptr,                   nullptr,                   nullptr                  }, // size=64 (to be implemented)
};

// ---------------------------------------------------------
// Reference implementation
// ---------------------------------------------------------

// Horizontal/vertical inverse transform. Produces a transposed output.
static void transform_inverse_naive(
  int16_t* __restrict output,
  const int16_t* __restrict input,
  const int16_t* __restrict coeffs,
  int height,
  int width,
  int shift,
  int skip_width,
  int skip_height)
{
  const int32_t round = (shift > 0) ? (1 << (shift - 1)) : 0;

  for (int x = 0; x < (width - skip_width); ++x) {
    for (int y = 0; y < height; ++y) {
      int32_t acc = 0;
      for (int k = 0; k < (height - skip_height); ++k) {
        acc += input[k * width + x] * coeffs[k * height + y];
      }
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[x * height + y] = result;
    }
  }

  for (int x = (width - skip_width); x < width; ++x) {
    for (int y = 0; y < height; ++y) {
      output[x * height + y] = 0;
    }
  }
}

// Vertical inverse transform. Produces a non-transposed output.
static void transform_inverse_vertical_naive(
  const int16_t* __restrict input,
  int16_t* __restrict output,
  tr_type_t type,
  int width,
  int height,
  int shift,
  int skip_width,
  int skip_height)
{
  NEON_ASSERT(shift > 0);
  const int16_t* __restrict coeffs_t = transform_coeffs_t[neon_fast_log2(height)][type];
  const int32_t round = 1 << (shift - 1);

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width - skip_width; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < height - skip_height; ++k) {
        acc += input[k * width + x] * coeffs_t[y * height + k];
      }
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * width + x] = result;
    }
    for (int x = width - skip_width; x < width; ++x) {
      output[y * width + x] = 0;
    }
  }
}

// Horizontal inverse transform. Produces a non-transposed output.
static void transform_inverse_horizontal_naive(
  const int16_t* __restrict input,
  int16_t* __restrict output,
  tr_type_t type,
  int width,
  int height,
  int shift,
  int skip_width)
{
  NEON_ASSERT(shift > 0);
  const int16_t* __restrict coeffs = transform_coeffs[neon_fast_log2(width)][type];
  const int32_t round = 1 << (shift - 1);

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < width - skip_width; ++k) {
        acc += input[y * width + k] * coeffs[k * width + x];
      }
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * width + x] = result;
    }
  }
}

// ---------------------------------------------------------
// Matrix multiply implementation
// ---------------------------------------------------------

template<int width, int height>
static void transform_inverse_vertical(
  const int16_t* __restrict input,
  int16_t* __restrict output,
  tr_type_t type,
  int shift,
  int skip_width,
  int skip_height)
{
  const int16_t* __restrict coeffs_t = transform_coeffs_t[neon_comptime_log2(height)][type];

  if constexpr (!NEON_ENABLE_INTRINSICS || width < 4) {
    NEON_ASSERT(shift > 0);
    const int32_t round = 1 << (shift - 1);
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width - skip_width; ++x) {
        int32_t acc = 0;
        for (int k = 0; k < height - skip_height; ++k) {
          acc += input[k * width + x] * coeffs_t[y * height + k];
        }
        acc += round;
        acc >>= shift;
        acc = MAX(acc, INT16_MIN);
        acc = MIN(acc, INT16_MAX);
        const int16_t result = (int16_t)acc;
        output[y * width + x] = result;
      }
      for (int x = width - skip_width; x < width; ++x) {
        output[y * width + x] = 0;
      }
    }
  }
  else if constexpr (height < 4) {
    NEON_ASSERT(shift > 0);
    NEON_ASSERT(width - skip_width >= 4);
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width - skip_width; x += 4) {
        int32x4_t vacc = vdupq_n_s32(0);
        for (int k = 0; k < height - skip_height; ++k) {
          vacc = vmlal_n_s16(vacc, vld1_s16(&input[k * width + x]), coeffs_t[y * height + k]);
        }
        const int32x4_t vshift = vdupq_n_s32(-shift);
        vacc = vqrshlq_s32(vacc, vshift);
        const int16x4_t vresult = vqmovn_s32(vacc);
        vst1_s16(&output[y * width + x], vresult);
      }
      for (int x = width - skip_width; x < width; ++x) {
        output[y * width + x] = 0;
      }
    }
  }
  else {
    NEON_ASSERT(shift > 0);
    NEON_ASSERT(width - skip_width >= 4);
    NEON_ASSERT(height - skip_height >= 4);
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width - skip_width; x += 4) {
        int32x4_t vacc = vdupq_n_s32(0);
        for (int k = 0; k < height - skip_height; k += 4) {
          vacc = vmlal_lane_s16(vacc, vld1_s16(&input[(k + 0) * width + x]), vld1_s16(&coeffs_t[y * height + k]), 0);
          vacc = vmlal_lane_s16(vacc, vld1_s16(&input[(k + 1) * width + x]), vld1_s16(&coeffs_t[y * height + k]), 1);
          vacc = vmlal_lane_s16(vacc, vld1_s16(&input[(k + 2) * width + x]), vld1_s16(&coeffs_t[y * height + k]), 2);
          vacc = vmlal_lane_s16(vacc, vld1_s16(&input[(k + 3) * width + x]), vld1_s16(&coeffs_t[y * height + k]), 3);
        }
        const int32x4_t vshift = vdupq_n_s32(-shift);
        vacc = vqrshlq_s32(vacc, vshift);
        const int16x4_t vresult = vqmovn_s32(vacc);
        vst1_s16(&output[y * width + x], vresult);
      }
      for (int x = width - skip_width; x < width; ++x) {
        output[y * width + x] = 0;
      }
    }
  }
}

template<int width>
static void transform_inverse_horizontal(
  const int16_t* __restrict input,
  int16_t* __restrict output,
  tr_type_t type,
  int height,
  int shift,
  int skip_width)
{
  const int16_t* __restrict coeffs = transform_coeffs[neon_comptime_log2(width)][type];

  if constexpr (!NEON_ENABLE_INTRINSICS || width < 4) {
    NEON_ASSERT(shift > 0);
    const int32_t round = 1 << (shift - 1);
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        int32_t acc = 0;
        for (int k = 0; k < width - skip_width; ++k) {
          acc += input[y * width + k] * coeffs[k * width + x];
        }
        acc += round;
        acc >>= shift;
        acc = MAX(acc, INT16_MIN);
        acc = MIN(acc, INT16_MAX);
        const int16_t result = (int16_t)acc;
        output[y * width + x] = result;
      }
    }
  }
  else {
    NEON_ASSERT(shift > 0);
    NEON_ASSERT(width - skip_width >= 4);
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; x += 4) {
        int32x4_t vacc = vdupq_n_s32(0);
        for (int k = 0; k < width - skip_width; k += 4) {
          vacc = vmlal_lane_s16(vacc, vld1_s16(&coeffs[(k + 0) * width + x]), vld1_s16(&input[y * width + k]), 0);
          vacc = vmlal_lane_s16(vacc, vld1_s16(&coeffs[(k + 1) * width + x]), vld1_s16(&input[y * width + k]), 1);
          vacc = vmlal_lane_s16(vacc, vld1_s16(&coeffs[(k + 2) * width + x]), vld1_s16(&input[y * width + k]), 2);
          vacc = vmlal_lane_s16(vacc, vld1_s16(&coeffs[(k + 3) * width + x]), vld1_s16(&input[y * width + k]), 3);
        }
        const int32x4_t vshift = vdupq_n_s32(-shift);
        vacc = vqrshlq_s32(vacc, vshift);
        const int16x4_t vresult = vqmovn_s32(vacc);
        vst1_s16(&output[y * width + x], vresult);
      }
    }
  }
}

// ---------------------------------------------------------
// Fast DCT2 implementation
// ---------------------------------------------------------

static void transform_inverse_vertical_Wx2_dct2(const int16_t* __restrict input, int16_t* __restrict output, int W, int shift)
{
  constexpr const auto& T = uvg_g_dct2_2x2;
  const int32_t round = 1 << (shift - 1);

  alignas(16) int32_t E[1][TR_MAX_WIDTH];
  alignas(16) int32_t O[1][TR_MAX_WIDTH];

  for (int x = 0; x < W; ++x) {
    E[0][x] = T[0][0] * (input[0 * W + x] + input[1 * W + x]);
    O[0][x] = T[0][1] * (input[0 * W + x] - input[1 * W + x]);
  }

  {
    constexpr int y = 0;
    for (int x = 0; x < W; ++x) {
      int32_t acc = E[0][x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
  {
    constexpr int y = 1;
    for (int x = 0; x < W; ++x) {
      int32_t acc = O[0][x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
}

static void transform_inverse_horizontal_2xH_dct2(const int16_t* __restrict input, int16_t* __restrict output, int H, int shift)
{
  constexpr const auto& T = uvg_g_dct2_2x2;
  const int32_t round = 1 << (shift - 1);

  alignas(16) int32_t E[TR_MAX_WIDTH][1];
  alignas(16) int32_t O[TR_MAX_WIDTH][1];

  for (int y = 0; y < H; ++y) {
    E[y][0] = T[0][0] * (input[y * 2 + 0] + input[y * 2 + 1]);
    O[y][0] = T[0][1] * (input[y * 2 + 0] - input[y * 2 + 1]);
  }

  for (int y = 0; y < H; ++y) {
    {
      constexpr int x = 0;
      int32_t acc = E[y][0];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * 2 + x] = result;
    }
    {
      constexpr int x = 1;
      int32_t acc = O[y][0];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * 2 + x] = result;
    }
  }
}

static void transform_inverse_vertical_Wx4_dct2(const int16_t* __restrict input, int16_t* __restrict output, int W, int shift)
{
  constexpr const auto& T = uvg_g_dct2_4x4;
  const int32_t round = 1 << (shift - 1);

  alignas(16) int32_t O[2][TR_MAX_WIDTH];
  alignas(16) int32_t E[2][TR_MAX_WIDTH];

  for (int y = 0; y < 2; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 1; k < 4; k += 2) {
        acc += T[k][y] * input[k * W + x];
      }
      O[y][x] = acc;
    }
  }
  for (int y = 0; y < 2; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 4; k += 2) {
        acc += T[k][y] * input[k * W + x];
      }
      E[y][x] = acc;
    }
  }

  for (int y = 0; y < 2; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = E[y][x] + O[y][x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
  for (int y = 2; y < 4; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = E[3 - y][x] - O[3 - y][x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
}

static void transform_inverse_horizontal_4xH_dct2(const int16_t* __restrict input, int16_t* __restrict output, int H, int shift)
{
  constexpr const auto& T = uvg_g_dct2_4x4;
  const int32_t round = 1 << (shift - 1);

  alignas(16) int32_t E[TR_MAX_WIDTH][2];
  alignas(16) int32_t O[TR_MAX_WIDTH][2];

  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < 2; ++x) {
      int32_t acc = 0;
      for (int k = 1; k < 4; k += 2) {
        acc += T[k][x] * input[y * 4 + k];
      }
      O[y][x] = acc;
    }
    for (int x = 0; x < 2; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 4; k += 2) {
        acc += T[k][x] * input[y * 4 + k];
      }
      E[y][x] = acc;
    }
  }

  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < 2; ++x) {
      int32_t acc = E[y][x] + O[y][x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * 4 + x] = result;
    }
    for (int x = 2; x < 4; ++x) {
      int32_t acc = E[y][3 - x] - O[y][3 - x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * 4 + x] = result;
    }
  }
}

static void transform_inverse_vertical_Wx8_dct2(const int16_t* __restrict input, int16_t* __restrict output, int W, int shift)
{
  constexpr const auto& T = uvg_g_dct2_8x8;
  const int32_t round = 1 << (shift - 1);

  alignas(16) int32_t O[4][TR_MAX_WIDTH];
  alignas(16) int32_t EO[2][TR_MAX_WIDTH];
  alignas(16) int32_t EE[2][TR_MAX_WIDTH];
  alignas(16) int32_t E[4][TR_MAX_WIDTH];

  for (int y = 0; y < 4; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 1; k < 8; k += 2) {
        acc += T[k][y] * input[k * W + x];
      }
      O[y][x] = acc;
    } 
  }
  for (int y = 0; y < 2; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 2; k < 8; k += 4) {
        acc += T[k][y] * input[k * W + x];
      }
      EO[y][x] = acc;
    }
  }
  for (int y = 0; y < 2; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 8; k += 4) {
        acc += T[k][y] * input[k * W + x];
      }
      EE[y][x] = acc;
    }
  }
  for (int y = 0; y < 2; ++y) {
    for (int x = 0; x < W; ++x) {
      E[y][x] = EE[y][x] + EO[y][x];
      E[y + 2][x] = EE[1 - y][x] - EO[1 - y][x];
    }
  }

  for (int y = 0; y < 4; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = E[y][x] + O[y][x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
  for (int y = 4; y < 8; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = E[7 - y][x] - O[7 - y][x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
}

static void transform_inverse_horizontal_8xH_dct2(const int16_t* __restrict input, int16_t* __restrict output, int H, int shift)
{
  constexpr const auto& T = uvg_g_dct2_8x8;
  const int32_t round = 1 << (shift - 1);

  alignas(16) int32_t E[TR_MAX_WIDTH][4];
  alignas(16) int32_t O[TR_MAX_WIDTH][4];

  for (int y = 0; y < H; ++y) {
    alignas(16) int32_t EE[2];
    alignas(16) int32_t EO[2];
    for (int x = 0; x < 4; ++x) {
      int32_t acc = 0;
      for (int k = 1; k < 8; k += 2) {
        acc += T[k][x] * input[y * 8 + k];
      }
      O[y][x] = acc;
    }
    for (int x = 0; x < 2; ++x) {
      int32_t acc = 0;
      for (int k = 2; k < 8; k += 4) {
        acc += T[k][x] * input[y * 8 + k];
      }
      EO[x] = acc;
    }
    for (int x = 0; x < 2; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 8; k += 4) {
        acc += T[k][x] * input[y * 8 + k];
      }
      EE[x] = acc;
    }
    for (int x = 0; x < 2; ++x) {
      E[y][x] = EE[x] + EO[x];
      E[y][x + 2] = EE[1 - x] - EO[1 - x];
    }
  }

  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < 4; ++x) {
      int32_t acc = E[y][x] + O[y][x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * 8 + x] = result;
    }
  }
  for (int y = 0; y < H; ++y) {
    for (int x = 4; x < 8; ++x) {
      int32_t acc = E[y][7 - x] - O[y][7 - x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * 8 + x] = result;
    }
  }
}

static void transform_inverse_vertical_Wx16_dct2(const int16_t* __restrict input, int16_t* __restrict output, int W, int shift)
{
  constexpr const auto& T = uvg_g_dct2_16x16;
  const int32_t round = 1 << (shift - 1);

  alignas(16) int32_t O[8][TR_MAX_WIDTH];
  alignas(16) int32_t EO[4][TR_MAX_WIDTH];
  alignas(16) int32_t EEO[2][TR_MAX_WIDTH];
  alignas(16) int32_t EEE[2][TR_MAX_WIDTH];
  alignas(16) int32_t EE[4][TR_MAX_WIDTH];
  alignas(16) int32_t E[8][TR_MAX_WIDTH];

  for (int y = 0; y < 8; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 1; k < 16; k += 2) {
        acc += T[k][y] * input[k * W + x];
      }
      O[y][x] = acc;
    }
  }
  for (int y = 0; y < 4; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 2; k < 16; k += 4) {
        acc += T[k][y] * input[k * W + x];
      }
      EO[y][x] = acc;
    }
  }
  for (int y = 0; y < 2; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 4; k < 16; k += 8) {
        acc += T[k][y] * input[k * W + x];
      }
      EEO[y][x] = acc;
    }
  }
  for (int y = 0; y < 2; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 16; k += 8) {
        acc += T[k][y] * input[k * W + x];
      }
      EEE[y][x] = acc;
    }
  }
  for (int y = 0; y < 2; ++y) {
    for (int x = 0; x < W; ++x) {
      EE[y][x] = EEE[y][x] + EEO[y][x];
      EE[y + 2][x] = EEE[1 - y][x] - EEO[1 - y][x];
    }
  }
  for (int y = 0; y < 4; ++y) {
    for (int x = 0; x < W; ++x) {
      E[y][x] = EE[y][x] + EO[y][x];
      E[y + 4][x] = EE[3 - y][x] - EO[3 - y][x];
    }
  }

  for (int y = 0; y < 8; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = E[y][x] + O[y][x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
  for (int y = 8; y < 16; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = E[15 - y][x] - O[15 - y][x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
}

static void transform_inverse_horizontal_16xH_dct2(const int16_t* __restrict input, int16_t* __restrict output, int H, int shift)
{
  constexpr const auto& T = uvg_g_dct2_16x16;
  const int32_t round = 1 << (shift - 1);

  alignas(16) int32_t E[TR_MAX_WIDTH][8];
  alignas(16) int32_t O[TR_MAX_WIDTH][8];

  for (int y = 0; y < H; ++y) {
    alignas(16) int32_t EE[4];
    alignas(16) int32_t EO[4];
    alignas(16) int32_t EEE[2];
    alignas(16) int32_t EEO[2];
    for (int x = 0; x < 8; ++x) {
      int32_t acc = 0;
      for (int k = 1; k < 16; k += 2) {
        acc += T[k][x] * input[y * 16 + k];
      }
      O[y][x] = acc;
    }
    for (int x = 0; x < 4; ++x) {
      int32_t acc = 0;
      for (int k = 2; k < 16; k += 4) {
        acc += T[k][x] * input[y * 16 + k];
      }
      EO[x] = acc;
    }
    for (int x = 0; x < 2; ++x) {
      int32_t acc = 0;
      for (int k = 4; k < 16; k += 8) {
        acc += T[k][x] * input[y * 16 + k];
      }
      EEO[x] = acc;
    }
    for (int x = 0; x < 2; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 16; k += 8) {
        acc += T[k][x] * input[y * 16 + k];
      }
      EEE[x] = acc;
    }
    for (int x = 0; x < 2; ++x) {
      EE[x] = EEE[x] + EEO[x];
      EE[x + 2] = EEE[1 - x] - EEO[1 - x];
    }
    for (int x = 0; x < 4; ++x) {
      E[y][x] = EE[x] + EO[x];
      E[y][x + 4] = EE[3 - x] - EO[3 - x];
    }
  }

  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < 8; ++x) {
      int32_t acc = E[y][x] + O[y][x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * 16 + x] = result;
    }
  }
  for (int y = 0; y < H; ++y) {
    for (int x = 8; x < 16; ++x) {
      int32_t acc = E[y][15 - x] - O[y][15 - x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * 16 + x] = result;
    }
  }
}

static void transform_inverse_vertical_Wx32_dct2(const int16_t* __restrict input, int16_t* __restrict output, int W, int shift)
{
  constexpr const auto& T = uvg_g_dct2_32x32;
  const int32_t round = 1 << (shift - 1);

  alignas(16) int32_t O[16][TR_MAX_WIDTH];
  alignas(16) int32_t EO[8][TR_MAX_WIDTH];
  alignas(16) int32_t EEO[4][TR_MAX_WIDTH];
  alignas(16) int32_t EEEO[2][TR_MAX_WIDTH];
  alignas(16) int32_t EEEE[2][TR_MAX_WIDTH];
  alignas(16) int32_t EEE[4][TR_MAX_WIDTH];
  alignas(16) int32_t EE[8][TR_MAX_WIDTH];
  alignas(16) int32_t E[16][TR_MAX_WIDTH];

  for (int y = 0; y < 16; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 1; k < 32; k += 2) {
        acc += T[k][y] * input[k * W + x];
      }
      O[y][x] = acc;
    }
  }
  for (int y = 0; y < 8; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 2; k < 32; k += 4) {
        acc += T[k][y] * input[k * W + x];
      }
      EO[y][x] = acc;
    }
  }
  for (int y = 0; y < 4; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 4; k < 32; k += 8) {
        acc += T[k][y] * input[k * W + x];
      }
      EEO[y][x] = acc;
    }
  }
  for (int y = 0; y < 2; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 8; k < 32; k += 16) {
        acc += T[k][y] * input[k * W + x];
      }
      EEEO[y][x] = acc;
    }
  }
  for (int y = 0; y < 2; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 32; k += 16) {
        acc += T[k][y] * input[k * W + x];
      }
      EEEE[y][x] = acc;
    }
  }
  for (int y = 0; y < 2; ++y) {
    for (int x = 0; x < W; ++x) {
      EEE[y][x] = EEEE[y][x] + EEEO[y][x];
      EEE[y + 2][x] = EEEE[1 - y][x] - EEEO[1 - y][x];
    }
  }
  for (int y = 0; y < 4; ++y) {
    for (int x = 0; x < W; ++x) {
      EE[y][x] = EEE[y][x] + EEO[y][x];
      EE[y + 4][x] = EEE[3 - y][x] - EEO[3 - y][x];
    }
  }
  for (int y = 0; y < 8; ++y) {
    for (int x = 0; x < W; ++x) {
      E[y][x] = EE[y][x] + EO[y][x];
      E[y + 8][x] = EE[7 - y][x] - EO[7 - y][x];
    }
  }

  for (int y = 0; y < 16; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = E[y][x] + O[y][x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
  for (int y = 16; y < 32; ++y) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = E[31 - y][x] - O[31 - y][x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
}

static void transform_inverse_horizontal_32xH_dct2(const int16_t* __restrict input, int16_t* __restrict output, int H, int shift)
{
  constexpr const auto& T = uvg_g_dct2_32x32;
  const int32_t round = 1 << (shift - 1);

  alignas(16) int32_t O[TR_MAX_WIDTH][16];
  alignas(16) int32_t E[TR_MAX_WIDTH][16];

  for (int y = 0; y < H; ++y) {
    alignas(16) int32_t EO[8];
    alignas(16) int32_t EEO[4];
    alignas(16) int32_t EEEO[2];
    alignas(16) int32_t EEEE[2];
    alignas(16) int32_t EEE[4];
    alignas(16) int32_t EE[8];
    for (int x = 0; x < 16; ++x) {
      int32_t acc = 0;
      for (int k = 1; k < 32; k += 2) {
        acc += T[k][x] * input[y * 32 + k];
      }
      O[y][x] = acc;
    }
    for (int x = 0; x < 8; ++x) {
      int32_t acc = 0;
      for (int k = 2; k < 32; k += 4) {
        acc += T[k][x] * input[y * 32 + k];
      }
      EO[x] = acc;
    }
    for (int x = 0; x < 4; ++x) {
      int32_t acc = 0;
      for (int k = 4; k < 32; k += 8) {
        acc += T[k][x] * input[y * 32 + k];
      }
      EEO[x] = acc;
    }
    for (int x = 0; x < 2; ++x) {
      int32_t acc = 0;
      for (int k = 8; k < 32; k += 16) {
        acc += T[k][x] * input[y * 32 + k];
      }
      EEEO[x] = acc;
    }
    for (int x = 0; x < 2; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 32; k += 16) {
        acc += T[k][x] * input[y * 32 + k];
      }
      EEEE[x] = acc;
    }
    for (int x = 0; x < 2; ++x) {
      EEE[x] = EEEE[x] + EEEO[x];
      EEE[x + 2] = EEEE[1 - x] - EEEO[1 - x];
    }
    for (int x = 0; x < 4; ++x) {
      EE[x] = EEE[x] + EEO[x];
      EE[x + 4] = EEE[3 - x] - EEO[3 - x];
    }
    for (int x = 0; x < 8; ++x) {
      E[y][x] = EE[x] + EO[x];
      E[y][x + 8] = EE[7 - x] - EO[7 - x];
    }
  }

  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < 16; ++x) {
      int32_t acc = E[y][x] + O[y][x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * 32 + x] = result;
    }
  }
  for (int y = 0; y < H; ++y) {
    for (int x = 16; x < 32; ++x) {
      int32_t acc = E[y][31 - x] - O[y][31 - x];
      acc += round;
      acc >>= shift;
      acc = MAX(acc, INT16_MIN);
      acc = MIN(acc, INT16_MAX);
      const int16_t result = (int16_t)acc;
      output[y * 32 + x] = result;
    }
  }
}

// ---------------------------------------------------------
// Dispatch
// ---------------------------------------------------------

template<int width, int height>
static void do_inverse_vertical_pass(
  const int16_t* __restrict input,
  int16_t* __restrict output,
  tr_type_t type,
  int shift,
  int skip_width,
  int skip_height)
{
  if (NEON_FORCE_FAST_DCT2 && type == DCT2) {
    switch (height) {
      case  1: transform_inverse_vertical<width, height>(input, output, type, shift, skip_width, skip_height); break;
      case  2: transform_inverse_vertical_Wx2_dct2(input, output, width, shift); break;
      case  4: transform_inverse_vertical_Wx4_dct2(input, output, width, shift); break;
      case  8: transform_inverse_vertical_Wx8_dct2(input, output, width, shift); break;
      case 16: transform_inverse_vertical_Wx16_dct2(input, output, width, shift); break;
      case 32: transform_inverse_vertical_Wx32_dct2(input, output, width, shift); break;
      default: NEON_UNREACHABLE("Unexpected height: %i\n", height); break;
    }
  }
  else if (NEON_PREFER_FAST_DCT2 && type == DCT2)
  {
    if constexpr (height == 16 && width >= 1) {
      transform_inverse_vertical_Wx16_dct2(input, output, width, shift);
    }
    else if constexpr (height == 32 && width >= 1) {
      transform_inverse_vertical_Wx32_dct2(input, output, width, shift);
    }
    else {
      transform_inverse_vertical<width, height>(input, output, type, shift, skip_width, skip_height);
    }
  }
  else {
    transform_inverse_vertical<width, height>(input, output, type, shift, skip_width, skip_height);
  }
}

template<int width, int height>
static void do_inverse_horizontal_pass(
  const int16_t* __restrict input,
  int16_t* __restrict output,
  tr_type_t type,
  int shift,
  int skip_width)
{
  if (NEON_FORCE_FAST_DCT2 && type == DCT2) {
    switch (width) {
      case  1: transform_inverse_horizontal<width>(input, output, type, height, shift, skip_width); break;
      case  2: transform_inverse_horizontal_2xH_dct2(input, output, height, shift); break;
      case  4: transform_inverse_horizontal_4xH_dct2(input, output, height, shift); break;
      case  8: transform_inverse_horizontal_8xH_dct2(input, output, height, shift); break;
      case 16: transform_inverse_horizontal_16xH_dct2(input, output, height, shift); break;
      case 32: transform_inverse_horizontal_32xH_dct2(input, output, height, shift); break;
      default: NEON_UNREACHABLE("Unexpected width: %i\n", width); break;
    }
  }
  else if (NEON_PREFER_FAST_DCT2 && type == DCT2) {
    if constexpr (width == 16 && height >= 1) {
      transform_inverse_horizontal_16xH_dct2(input, output, height, shift);
    }
    else if constexpr (width == 32 && height >= 1) {
      transform_inverse_horizontal_32xH_dct2(input, output, height, shift);
    }
    else {
      transform_inverse_horizontal<width>(input, output, type, height, shift, skip_width);
    }
  }
  else {
    transform_inverse_horizontal<width>(input, output, type, height, shift, skip_width);
  }
}

// ---------------------------------------------------------
// MTS
// ---------------------------------------------------------

template<int width, int height>
static void mts_inverse(
  const cu_info_t* __restrict tu,
  const int16_t* __restrict input,
  int16_t* __restrict output,
  color_t color,
  int bitdepth,
  bool is_isp,
  bool is_lfnst,
  bool explicit_mts,
  bool implicit_mts)
{
  if constexpr (width == 1) {
    const tr_type_t vertical_type = get_vertical_transform_type(tu, height, color, is_isp, is_lfnst, explicit_mts, implicit_mts);
    const int horizontal_skip = 0;
    const int vertical_skip = get_vertical_skip(width, height, vertical_type, is_lfnst);
    const int vertical_shift = 21 - bitdepth;
    do_inverse_vertical_pass<width, height>(input, output, vertical_type, vertical_shift, horizontal_skip, vertical_skip);
  }
  else if constexpr (height == 1) {
    const tr_type_t horizontal_type = get_horizontal_transform_type(tu, width, color, is_isp, is_lfnst, explicit_mts, implicit_mts);
    const int horizontal_skip = get_horizontal_skip(width, height, horizontal_type, is_lfnst);
    const int horizontal_shift = 21 - bitdepth;
    do_inverse_horizontal_pass<width, height>(input, output, horizontal_type, horizontal_shift, horizontal_skip);
  }
  else {
    const tr_type_t horizontal_type = width <= 2 ? DCT2 : get_horizontal_transform_type(tu, width, color, is_isp, is_lfnst, explicit_mts, implicit_mts);
    const tr_type_t vertical_type = height <= 2 ? DCT2 : get_vertical_transform_type(tu, height, color, is_isp, is_lfnst, explicit_mts, implicit_mts);
    const int horizontal_skip = width <= 4 ? 0 : get_horizontal_skip(width, height, horizontal_type, is_lfnst);
    const int vertical_skip = height <= 4 ? 0 : get_vertical_skip(width, height, vertical_type, is_lfnst);
    const int vertical_shift = 7;
    const int horizontal_shift = 20 - bitdepth;
    alignas(16) int16_t tmp[TR_MAX_WIDTH * TR_MAX_WIDTH];
    do_inverse_vertical_pass<width, height>(input, tmp, vertical_type, vertical_shift, horizontal_skip, vertical_skip);
    do_inverse_horizontal_pass<width, height>(tmp, output, horizontal_type, horizontal_shift, horizontal_skip);
  }
}

extern "C" void uvg_mts_idct_neon(
  int8_t bitdepth,
  color_t color,
  const cu_info_t* tu,
  int8_t width,
  int8_t height,
  const int16_t* input,
  int16_t* output,
  int8_t mts_type)
{
  const bool is_isp = tu->type == CU_INTRA && tu->intra.isp_mode && color == COLOR_Y ? tu->intra.isp_mode : 0;
  const bool is_lfnst = (color == COLOR_Y ? tu->lfnst_idx : tu->cr_lfnst_idx);
  const bool explicit_mts = mts_type == UVG_MTS_BOTH || (tu->type == CU_INTRA ? mts_type == UVG_MTS_INTRA : (mts_type == UVG_MTS_INTER && tu->type == CU_INTER));
  const bool implicit_mts = tu->type == CU_INTRA && (mts_type == UVG_MTS_IMPLICIT || mts_type == UVG_MTS_INTER);

  switch (width) {
    case 1: {
        switch (height) {
          case 16: mts_inverse<1, 16>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 32: mts_inverse<1, 32>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          default: NEON_UNREACHABLE("Unexpected transform size: %ix%i\n", width, height);
        }
      }
      break;
    case 2: {
        switch (height) {
          case  4: mts_inverse<2,  4>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  8: mts_inverse<2,  8>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 16: mts_inverse<2, 16>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 32: mts_inverse<2, 32>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          default: NEON_UNREACHABLE("Unexpected transform size: %ix%i\n", width, height);
        }
      }
      break;
    case 4: {
        switch (height) {
          case  2: mts_inverse<4,  2>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  4: mts_inverse<4,  4>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  8: mts_inverse<4,  8>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 16: mts_inverse<4, 16>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 32: mts_inverse<4, 32>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          default: NEON_UNREACHABLE("Unexpected transform size: %ix%i\n", width, height);
        }
      }
      break;
    case 8: {
        switch (height) {
          case  2: mts_inverse<8,  2>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  4: mts_inverse<8,  4>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  8: mts_inverse<8,  8>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 16: mts_inverse<8, 16>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 32: mts_inverse<8, 32>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          default: NEON_UNREACHABLE("Unexpected transform size: %ix%i\n", width, height);
        }
      }
      break;
    case 16: {
        switch (height) {
          case  1: mts_inverse<16,  1>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  2: mts_inverse<16,  2>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  4: mts_inverse<16,  4>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  8: mts_inverse<16,  8>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 16: mts_inverse<16, 16>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 32: mts_inverse<16, 32>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          default: NEON_UNREACHABLE("Unexpected transform size: %ix%i\n", width, height);
        }
      }
      break;
    case 32: {
        switch (height) {
          case  1: mts_inverse<32,  1>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  2: mts_inverse<32,  2>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  4: mts_inverse<32,  4>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  8: mts_inverse<32,  8>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 16: mts_inverse<32, 16>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 32: mts_inverse<32, 32>(tu, input, output, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          default: NEON_UNREACHABLE("Unexpected transform size: %ix%i\n", width ,height);
        }
      }
      break;
    default: {
        NEON_UNREACHABLE("Unexpected transform size: %ix%i\n", width, height);
      }
  }
}

// ---------------------------------------------------------
// Legacy API
// ---------------------------------------------------------

template<int N>
static inline void legacy_idct_NxN(const int16_t* __restrict input, int16_t* __restrict output, int bitdepth)
{
  const int vertical_shift = 7;
  const int horizontal_shift = 20 - bitdepth;
  const int vertical_skip = 0;
  const int horizontal_skip = 0;

  alignas(16) int16_t tmp[N * N];
  do_inverse_vertical_pass<N, N>(input, tmp, DCT2, vertical_shift, horizontal_skip, vertical_skip);
  do_inverse_horizontal_pass<N, N>(tmp, output, DCT2, horizontal_shift, horizontal_skip);
}

extern "C" void uvg_idct_4x4_neon(int8_t bitdepth, const int16_t* input, int16_t* output)
{
  legacy_idct_NxN<4>(input, output, bitdepth);
};

extern "C" void uvg_idct_8x8_neon(int8_t bitdepth, const int16_t* input, int16_t* output)
{
  legacy_idct_NxN<8>(input, output, bitdepth);
};

extern "C" void uvg_idct_16x16_neon(int8_t bitdepth, const int16_t* input, int16_t* output)
{
  legacy_idct_NxN<16>(input, output, bitdepth);
};

extern "C" void uvg_idct_32x32_neon(int8_t bitdepth, const int16_t* input, int16_t* output)
{
  legacy_idct_NxN<32>(input, output, bitdepth);
};
