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

// Horizontal/vertical forward transform. Produces a transposed output.
static void transform_forward_naive(
  const int16_t* __restrict input,
  int16_t* __restrict output,
  const int16_t* __restrict coeffs,
  int width,
  int height,
  int shift,
  int skip_height,
  int skip_width)
{
  const int32_t round = (shift > 0) ? (1 << (shift - 1)) : 0;

  for (int x = 0; x < (width - skip_width); ++x) {
    for (int y = 0; y < (height - skip_height); ++y) {
      int32_t acc = 0;
      for (int i = 0; i < width; ++i) {
        acc += input[y * width + i] * coeffs[x * width + i];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[x * height + y] = result;
    }
  }

  for (int x = 0; x < (width - skip_width); ++x) {
    for (int y = (height - skip_height); y < height; ++y) {
      output[x * height + y] = 0;
    }
  }

  for (int x = (width - skip_width); x < width; ++x) {
    for (int y = 0; y < height; ++y) {
      output[x * height + y] = 0;
    }
  }
}

// Horizontal forward transform. Produces a non-transposed output.
static void transform_forward_horizontal_naive(
  const int16_t* __restrict input,
  int16_t* __restrict output,
  tr_type_t type,
  int width,
  int height,
  int shift,
  int skip_width)
{
  NEON_ASSERT(shift >= 0);
  const int16_t* __restrict coeffs_t = transform_coeffs_t[neon_fast_log2(width)][type];
  const int32_t round = (shift > 0) ? (1 << (shift - 1)) : 0;

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width - skip_width; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < width; ++k) {
        acc += coeffs_t[k * width + x] * input[y * width + k];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * width + x] = result;
    }
    for (int x = width - skip_width; x < width; ++x) {
      output[y * width + x] = 0;
    }
  }
}

// Vertical forward transform. Produces a non-transposed output.
static void transform_forward_vertical_naive(
  const int16_t* __restrict input,
  int16_t* __restrict output,
  tr_type_t type,
  int width,
  int height,
  int shift,
  int skip_height)
{
  NEON_ASSERT(shift > 0);
  const int16_t* __restrict coeffs = transform_coeffs[neon_fast_log2(height)][type];
  const int32_t round = 1 << (shift - 1);

  for (int y = 0; y < height - skip_height; ++y) {
    for (int x = 0; x < width; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < height; ++k) {
        acc += input[k * width + x] * coeffs[y * height + k];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * width + x] = result;
    }
  }

  for (int y = height - skip_height; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      output[y * width + x] = 0;
    }
  }
}

// ---------------------------------------------------------
// Matrix multiply implementation
// ---------------------------------------------------------

template<int width>
static void transform_forward_horizontal(
  const int16_t* __restrict input,
  int16_t* __restrict output,
  tr_type_t type,
  int height,
  int shift,
  int skip_width)
{
  const int16_t* __restrict coeffs_t = transform_coeffs_t[neon_comptime_log2(width)][type];

  if constexpr (!NEON_ENABLE_INTRINSICS || width < 4) {
    NEON_ASSERT(shift >= 0);
    const int32_t round = (shift > 0) ? (1 << (shift - 1)) : 0;
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width - skip_width; ++x) {
        int32_t acc = 0;
        for (int k = 0; k < width; ++k) {
          acc += coeffs_t[k * width + x] * input[y * width + k];
        }
        acc += round;
        acc >>= shift;
        const int16_t result = (int16_t)acc;
        output[y * width + x] = result;
      }
      for (int x = width - skip_width; x < width; ++x) {
        output[y * width + x] = 0;
      }
    }
  }
  else {
    NEON_ASSERT(shift > 0);
    const int32x4_t vshift = vdupq_n_s32(-shift);
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width - skip_width; x += 4) {
        int32x4_t vacc = vdupq_n_s32(0);
        for (int k = 0; k < width; k += 4) {
          const int16x4_t vinput = vld1_s16(&input[y * width + k]);
          vacc = vmlal_lane_s16(vacc, vld1_s16(&coeffs_t[(k + 0) * width + x]), vinput, 0);
          vacc = vmlal_lane_s16(vacc, vld1_s16(&coeffs_t[(k + 1) * width + x]), vinput, 1);
          vacc = vmlal_lane_s16(vacc, vld1_s16(&coeffs_t[(k + 2) * width + x]), vinput, 2);
          vacc = vmlal_lane_s16(vacc, vld1_s16(&coeffs_t[(k + 3) * width + x]), vinput, 3);
        }
        vacc = vqrshlq_s32(vacc, vshift);
        const int16x4_t vresult = vmovn_s32(vacc);
        vst1_s16(&output[y * width + x], vresult);
      }
      for (int x = width - skip_width; x < width; ++x) {
        output[y * width + x] = 0;
      }
    }
  }
}

template<int width, int height>
static void transform_forward_vertical(
  const int16_t* __restrict input,
  int16_t* __restrict output,
  tr_type_t type,
  int shift,
  int skip_height)
{
  const int16_t* __restrict coeffs = transform_coeffs[neon_comptime_log2(height)][type];

  if constexpr (!NEON_ENABLE_INTRINSICS || width < 4) {
    NEON_ASSERT(shift > 0);
    const int32_t round = 1 << (shift - 1);
    for (int y = 0; y < height - skip_height; ++y) {
      for (int x = 0; x < width; ++x) {
        int32_t acc = 0;
        for (int k = 0; k < height; ++k) {
          acc += input[k * width + x] * coeffs[y * height + k];
        }
        acc += round;
        acc >>= shift;
        const int16_t result = (int16_t)acc;
        output[y * width + x] = result;
      }
    }
    for (int y = height - skip_height; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        output[y * width + x] = 0;
      }
    }
  }
  else if constexpr (height < 4) {
    NEON_ASSERT(shift > 0);
    for (int y = 0; y < height - skip_height; ++y) {
      for (int x = 0; x < width; x += 4) {
        int32x4_t vacc = vdupq_n_s32(0);
        for (int k = 0; k < height; ++k) {
          vacc = vmlal_n_s16(vacc, vld1_s16(&input[k * width + x]), coeffs[y * height + k]);
        }
        const int32x4_t vshift = vdupq_n_s32(-shift);
        vacc = vqrshlq_s32(vacc, vshift);
        const int16x4_t vresult = vmovn_s32(vacc);
        vst1_s16(&output[y * width + x], vresult);
      }
    }
    for (int y = height - skip_height; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        output[y * width + x] = 0;
      }
    }
  }
  else {
    NEON_ASSERT(shift > 0);
    for (int y = 0; y < height - skip_height; ++y) {
      for (int x = 0; x < width; x += 4) {
        int32x4_t vacc = vdupq_n_s32(0);
        for (int k = 0; k < height; k += 4) {
          const int16x4_t vcoeff = vld1_s16(&coeffs[y * height + k]);
          vacc = vmlal_lane_s16(vacc, vld1_s16(&input[(k + 0) * width + x]), vcoeff, 0);
          vacc = vmlal_lane_s16(vacc, vld1_s16(&input[(k + 1) * width + x]), vcoeff, 1);
          vacc = vmlal_lane_s16(vacc, vld1_s16(&input[(k + 2) * width + x]), vcoeff, 2);
          vacc = vmlal_lane_s16(vacc, vld1_s16(&input[(k + 3) * width + x]), vcoeff, 3);
        }
        const int32x4_t vshift = vdupq_n_s32(-shift);
        vacc = vqrshlq_s32(vacc, vshift);
        const int16x4_t vresult = vmovn_s32(vacc);
        vst1_s16(&output[y * width + x], vresult);
      }
    }
    for (int y = height - skip_height; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        output[y * width + x] = 0;
      }
    }
  }
}

// ---------------------------------------------------------
// Fast DCT2 implementation
// ---------------------------------------------------------

static void transform_forward_horizontal_2xH_dct2(const int16_t* __restrict input, int16_t* __restrict output, int H, int shift)
{
  constexpr const auto& T = uvg_g_dct2_2x2;
  const int32_t round = (shift > 0) ? (1 << (shift - 1)) : 0;

  alignas(16) int32_t E[1][TR_MAX_WIDTH];
  alignas(16) int32_t O[1][TR_MAX_WIDTH];

  for (int y = 0; y < H; ++y) {
    E[0][y] = input[y * 2 + 0] + input[y * 2 + 1];
    O[0][y] = input[y * 2 + 0] - input[y * 2 + 1];
  }

  // Avoid a separate transpose step by writing the output in an already transposed order.

  for (int y = 0; y < H; ++y) {
    {
      constexpr int x = 0;
      int32_t acc = T[x][0] * E[0][y];
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * 2 + x] = result;
    }
    {
      constexpr int x = 1;
      int32_t acc = T[x][0] * O[0][y];
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * 2 + x] = result;
    }
  }
}

static void transform_forward_vertical_Wx2_dct2(const int16_t* __restrict input, int16_t* __restrict output, int W, int shift)
{
  constexpr const auto& T = uvg_g_dct2_2x2;
  const int32_t round = (shift > 0) ? (1 << (shift - 1)) : 0;

  alignas(16) int32_t E[1][TR_MAX_WIDTH];
  alignas(16) int32_t O[1][TR_MAX_WIDTH];

  for (int x = 0; x < W; ++x) {
    E[0][x] = input[0 * W + x] + input[1 * W + x];
    O[0][x] = input[0 * W + x] - input[1 * W + x];
  }

  for (int x = 0; x < W; ++x) {
    {
      constexpr int y = 0;
      int32_t acc = T[y][0] * E[0][x];
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
    {
      constexpr int y = 1;
      int32_t acc = T[y][0] * O[0][x];
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
}

static void transform_forward_horizontal_4xH_dct2(const int16_t* __restrict input, int16_t* __restrict output, int H, int shift)
{
  constexpr const auto& T = uvg_g_dct2_4x4;
  const int32_t round = (shift > 0) ? (1 << (shift - 1)) : 0;

  alignas(16) int32_t E[2][TR_MAX_WIDTH];
  alignas(16) int32_t O[2][TR_MAX_WIDTH];

  for (int y = 0; y < H; ++y) {
    for (int k = 0; k < 2; ++k) {
      E[k][y] = input[y * 4 + k] + input[y * 4 + (3 - k)];
      O[k][y] = input[y * 4 + k] - input[y * 4 + (3 - k)];
    }
  }

  alignas(16) int16_t tmp[4 * TR_MAX_WIDTH];

  for (int x = 0; x < 4; x += 2) {
    for (int y = 0; y < H; ++y) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[x][k] * E[k][y];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      tmp[x * H + y] = result;
    }
  }
  for (int x = 1; x < 4; x += 2) {
    for (int y = 0; y < H; ++y) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[x][k] * O[k][y];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      tmp[x * H + y] = result;
    }
  }

  uvg_transpose_16bit_neon(tmp, output, H, 4);
}

static void transform_forward_vertical_Wx4_dct2(const int16_t* __restrict input, int16_t* __restrict output, int W, int shift)
{
  constexpr const auto& T = uvg_g_dct2_4x4;
  const int32_t round = (shift > 0) ? (1 << (shift - 1)) : 0;

  alignas(16) int32_t E[2][TR_MAX_WIDTH];
  alignas(16) int32_t O[2][TR_MAX_WIDTH];

  for (int x = 0; x < W; ++x) {
    for (int k = 0; k < 2; ++k) {
      E[k][x] = input[k * W + x] + input[(3 - k) * W + x];
      O[k][x] = input[k * W + x] - input[(3 - k) * W + x];
    }
  }

  for (int y = 0; y < 4; y += 2) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[y][k] * E[k][x];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
  for (int y = 1; y < 4; y += 2) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[y][k] * O[k][x];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
}

static void transform_forward_horizontal_8xH_dct2(const int16_t* __restrict input, int16_t* __restrict output, int H, int shift)
{
  constexpr const auto& T = uvg_g_dct2_8x8;
  const int32_t round = (shift > 0) ? (1 << (shift - 1)) : 0;

  alignas(16) int32_t E[4][TR_MAX_WIDTH];
  alignas(16) int32_t O[4][TR_MAX_WIDTH];
  alignas(16) int32_t EE[2][TR_MAX_WIDTH];
  alignas(16) int32_t EO[2][TR_MAX_WIDTH];

  for (int y = 0; y < H; ++y) {
    for (int k = 0; k < 4; ++k) {
      E[k][y] = input[y * 8 + k] + input[y * 8 + (7 - k)];
      O[k][y] = input[y * 8 + k] - input[y * 8 + (7 - k)];
    }
    for (int k = 0; k < 2; ++k) {
      EE[k][y] = E[k][y] + E[3 - k][y];
      EO[k][y] = E[k][y] - E[3 - k][y];
    }
  }

  alignas(16) int16_t tmp[8 * TR_MAX_WIDTH];

  for (int x = 0; x < 8; x += 4) {
    for (int y = 0; y < H; ++y) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[x][k] * EE[k][y];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      tmp[x * H + y] = result;
    }
  }
  for (int x = 2; x < 8; x += 4) {
    for (int y = 0; y < H; ++y) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[x][k] * EO[k][y];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      tmp[x * H + y] = result;
    }
  }
  for (int x = 1; x < 8; x += 2) {
    for (int y = 0; y < H; ++y) {
      int32_t acc = 0;
      for (int k = 0; k < 4; ++k) {
        acc += T[x][k] * O[k][y];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      tmp[x * H + y] = result;
    }
  }

  uvg_transpose_16bit_neon(tmp, output, H, 8);
}

static void transform_forward_vertical_Wx8_dct2(const int16_t* __restrict input, int16_t* __restrict output, int W, int shift)
{
  constexpr const auto& T = uvg_g_dct2_8x8;
  const int32_t round = (shift > 0) ? (1 << (shift - 1)) : 0;

  alignas(16) int32_t E[4][TR_MAX_WIDTH];
  alignas(16) int32_t O[4][TR_MAX_WIDTH];
  alignas(16) int32_t EE[2][TR_MAX_WIDTH];
  alignas(16) int32_t EO[2][TR_MAX_WIDTH];

  for (int x = 0; x < W; ++x) {
    for (int k = 0; k < 4; ++k) {
      E[k][x] = input[k * W + x] + input[(7 - k) * W + x];
      O[k][x] = input[k * W + x] - input[(7 - k) * W + x];
    }
    for (int k = 0; k < 2; ++k) {
      EE[k][x] = E[k][x] + E[3 - k][x];
      EO[k][x] = E[k][x] - E[3 - k][x];
    }
  }

  for (int y = 0; y < 8; y += 4) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[y][k] * EE[k][x];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
  for (int y = 2; y < 8; y += 4) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[y][k] * EO[k][x];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
  for (int y = 1; y < 8; y += 2) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 4; ++k) {
        acc += T[y][k] * O[k][x];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
}

static void transform_forward_horizontal_16xH_dct2(const int16_t* __restrict input, int16_t* __restrict output, int H, int shift)
{
  constexpr const auto& T = uvg_g_dct2_16x16;
  const int32_t round = (shift > 0) ? (1 << (shift - 1)) : 0;

  alignas(16) int32_t E[8][TR_MAX_WIDTH];
  alignas(16) int32_t O[8][TR_MAX_WIDTH];
  alignas(16) int32_t EE[4][TR_MAX_WIDTH];
  alignas(16) int32_t EO[4][TR_MAX_WIDTH];
  alignas(16) int32_t EEE[2][TR_MAX_WIDTH];
  alignas(16) int32_t EEO[2][TR_MAX_WIDTH];

  for (int y = 0; y < H; ++y) {
    for (int k = 0; k < 8; ++k) {
      E[k][y] = input[y * 16 + k] + input[y * 16 + (15 - k)];
      O[k][y] = input[y * 16 + k] - input[y * 16 + (15 - k)];
    }
    for (int k = 0; k < 4; ++k) {
      EE[k][y] = E[k][y] + E[7 - k][y];
      EO[k][y] = E[k][y] - E[7 - k][y];
    }
    for (int k = 0; k < 2; ++k) {
      EEE[k][y] = EE[k][y] + EE[3 - k][y];
      EEO[k][y] = EE[k][y] - EE[3 - k][y];
    }
  }

  alignas(16) int16_t tmp[16 * TR_MAX_WIDTH];

  for (int x = 0; x < 16; x += 8) {
    for (int y = 0; y < H; ++y) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[x][k] * EEE[k][y];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      tmp[x * H + y] = result;
    }
  }
  for (int x = 4; x < 16; x += 8) {
    for (int y = 0; y < H; ++y) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[x][k] * EEO[k][y];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      tmp[x * H + y] = result;
    }
  }
  for (int x = 2; x < 16; x += 4) {
    for (int y = 0; y < H; ++y) {
      int32_t acc = 0;
      for (int k = 0; k < 4; ++k) {
        acc += T[x][k] * EO[k][y];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      tmp[x * H + y] = result;
    }
  }
  for (int x = 1; x < 16; x += 2) {
    for (int y = 0; y < H; ++y) {
      int32_t acc = 0;
      for (int k = 0; k < 8; ++k) {
        acc += T[x][k] * O[k][y];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      tmp[x * H + y] = result;
    }
  }

  uvg_transpose_16bit_neon(tmp, output, H, 16);
}

static void transform_forward_vertical_Wx16_dct2(const int16_t* __restrict input, int16_t* __restrict output, int W, int shift)
{
  constexpr const auto& T = uvg_g_dct2_16x16;
  const int32_t round = (shift > 0) ? (1 << (shift - 1)) : 0;

  alignas(16) int32_t E[8][TR_MAX_WIDTH];
  alignas(16) int32_t O[8][TR_MAX_WIDTH];
  alignas(16) int32_t EE[4][TR_MAX_WIDTH];
  alignas(16) int32_t EO[4][TR_MAX_WIDTH];
  alignas(16) int32_t EEE[2][TR_MAX_WIDTH];
  alignas(16) int32_t EEO[2][TR_MAX_WIDTH];

  for (int x = 0; x < W; ++x) {
    for (int k = 0; k < 8; ++k) {
      E[k][x] = input[k * W + x] + input[(15 - k) * W + x];
      O[k][x] = input[k * W + x] - input[(15 - k) * W + x];
    }
    for (int k = 0; k < 4; ++k) {
      EE[k][x] = E[k][x] + E[7 - k][x];
      EO[k][x] = E[k][x] - E[7 - k][x];
    }
    for (int k = 0; k < 2; ++k) {
      EEE[k][x] = EE[k][x] + EE[3 - k][x];
      EEO[k][x] = EE[k][x] - EE[3 - k][x];
    }
  }

  for (int y = 0; y < 16; y += 8) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[y][k] * EEE[k][x];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
  for (int y = 4; y < 16; y += 8) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[y][k] * EEO[k][x];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
  for (int y = 2; y < 16; y += 4) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 4; ++k) {
        acc += T[y][k] * EO[k][x];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
  for (int y = 1; y < 16; y += 2) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 8; ++k) {
        acc += T[y][k] * O[k][x];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
}

static void transform_forward_horizontal_32xH_dct2(const int16_t* __restrict input, int16_t* __restrict output, int H, int shift)
{
  constexpr const auto& T = uvg_g_dct2_32x32;
  const int32_t round = (shift > 0) ? (1 << (shift - 1)) : 0;

  alignas(16) int32_t E[16][TR_MAX_WIDTH];
  alignas(16) int32_t O[16][TR_MAX_WIDTH];
  alignas(16) int32_t EE[8][TR_MAX_WIDTH];
  alignas(16) int32_t EO[8][TR_MAX_WIDTH];
  alignas(16) int32_t EEE[4][TR_MAX_WIDTH];
  alignas(16) int32_t EEO[4][TR_MAX_WIDTH];
  alignas(16) int32_t EEEE[2][TR_MAX_WIDTH];
  alignas(16) int32_t EEEO[2][TR_MAX_WIDTH];

  for (int y = 0; y < H; ++y) {
    for (int k = 0; k < 16; ++k) {
      E[k][y] = input[y * 32 + k] + input[y * 32 + (31 - k)];
      O[k][y] = input[y * 32 + k] - input[y * 32 + (31 - k)];
    }
    for (int k = 0; k < 8; ++k) {
      EE[k][y] = E[k][y] + E[15 - k][y];
      EO[k][y] = E[k][y] - E[15 - k][y];
    }
    for (int k = 0; k < 4; ++k) {
      EEE[k][y] = EE[k][y] + EE[7 - k][y];
      EEO[k][y] = EE[k][y] - EE[7 - k][y];
    }
    for (int k = 0; k < 2; ++k) {
      EEEE[k][y] = EEE[k][y] + EEE[3 - k][y];
      EEEO[k][y] = EEE[k][y] - EEE[3 - k][y];
    }
  }

  alignas(16) int16_t tmp[32 * TR_MAX_WIDTH];

  for (int x = 0; x < 32; x += 16) {
    for (int y = 0; y < H; ++y) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[x][k] * EEEE[k][y];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      tmp[x * H + y] = result;
    }
  }
  for (int x = 8; x < 32; x += 16) {
    for (int y = 0; y < H; ++y) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[x][k] * EEEO[k][y];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      tmp[x * H + y] = result;
    }
  }
  for (int x = 4; x < 32; x += 8) {
    for (int y = 0; y < H; ++y) {
      int32_t acc = 0;
      for (int k = 0; k < 4; ++k) {
        acc += T[x][k] * EEO[k][y];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      tmp[x * H + y] = result;
    }
  }
  for (int x = 2; x < 32; x += 4) {
    for (int y = 0; y < H; ++y) {
      int32_t acc = 0;
      for (int k = 0; k < 8; ++k) {
        acc += T[x][k] * EO[k][y];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      tmp[x * H + y] = result;
    }
  }
  for (int x = 1; x < 32; x += 2) {
    for (int y = 0; y < H; ++y) {
      int32_t acc = 0;
      for (int k = 0; k < 16; ++k) {
        acc += T[x][k] * O[k][y];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      tmp[x * H + y] = result;
    }
  }

  uvg_transpose_16bit_neon(tmp, output, H, 32);
}

static void transform_forward_vertical_Wx32_dct2(const int16_t* __restrict input, int16_t* __restrict output, int W, int shift)
{
  constexpr const auto& T = uvg_g_dct2_32x32;
  const int32_t round = (shift > 0) ? (1 << (shift - 1)) : 0;

  alignas(16) int32_t E[16][TR_MAX_WIDTH];
  alignas(16) int32_t O[16][TR_MAX_WIDTH];
  alignas(16) int32_t EE[8][TR_MAX_WIDTH];
  alignas(16) int32_t EO[8][TR_MAX_WIDTH];
  alignas(16) int32_t EEE[4][TR_MAX_WIDTH];
  alignas(16) int32_t EEO[4][TR_MAX_WIDTH];
  alignas(16) int32_t EEEE[2][TR_MAX_WIDTH];
  alignas(16) int32_t EEEO[2][TR_MAX_WIDTH];

  for (int x = 0; x < W; ++x) {
    for (int k = 0; k < 16; ++k) {
      E[k][x] = input[k * W + x] + input[(31 - k) * W + x];
      O[k][x] = input[k * W + x] - input[(31 - k) * W + x];
    }
    for (int k = 0; k < 8; ++k) {
      EE[k][x] = E[k][x] + E[15 - k][x];
      EO[k][x] = E[k][x] - E[15 - k][x];
    }
    for (int k = 0; k < 4; ++k) {
      EEE[k][x] = EE[k][x] + EE[7 - k][x];
      EEO[k][x] = EE[k][x] - EE[7 - k][x];
    }
    for (int k = 0; k < 2; ++k) {
      EEEE[k][x] = EEE[k][x] + EEE[3 - k][x];
      EEEO[k][x] = EEE[k][x] - EEE[3 - k][x];
    }
  }

  for (int y = 0; y < 32; y += 16) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[y][k] * EEEE[k][x];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
  for (int y = 8; y < 32; y += 16) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 2; ++k) {
        acc += T[y][k] * EEEO[k][x];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
  for (int y = 4; y < 32; y += 8) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 4; ++k) {
        acc += T[y][k] * EEO[k][x];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
  for (int y = 2; y < 32; y += 4) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 8; ++k) {
        acc += T[y][k] * EO[k][x];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
  for (int y = 1; y < 32; y += 2) {
    for (int x = 0; x < W; ++x) {
      int32_t acc = 0;
      for (int k = 0; k < 16; ++k) {
        acc += T[y][k] * O[k][x];
      }
      acc += round;
      acc >>= shift;
      const int16_t result = (int16_t)acc;
      output[y * W + x] = result;
    }
  }
}

// ---------------------------------------------------------
// Dispatch
// ---------------------------------------------------------

template<int width, int height>
static void do_forward_horizontal_pass(const int16_t* __restrict input, int16_t* __restrict output, tr_type_t type, int shift, int skip_width)
{
  if (NEON_FORCE_FAST_DCT2 && type == DCT2 && width >= 2) {
    switch (width) {
      case  2: transform_forward_horizontal_2xH_dct2(input, output, height, shift); break;
      case  4: transform_forward_horizontal_4xH_dct2(input, output, height, shift); break;
      case  8: transform_forward_horizontal_8xH_dct2(input, output, height, shift); break;
      case 16: transform_forward_horizontal_16xH_dct2(input, output, height, shift); break;
      case 32: transform_forward_horizontal_32xH_dct2(input, output, height, shift); break;
      default: NEON_UNREACHABLE("Unexpected width: %i\n", width);
    }
  }
  else if (NEON_PREFER_FAST_DCT2 && type == DCT2) {
    if constexpr (width == 16 && height >= 2) {
      transform_forward_horizontal_16xH_dct2(input, output, height, shift);
    }
    else if constexpr (width == 32 && height >= 4) {
      transform_forward_horizontal_32xH_dct2(input, output, height, shift);
    }
    else {
      transform_forward_horizontal<width>(input, output, type, height, shift, skip_width);
    }
  }
  else {
    transform_forward_horizontal<width>(input, output, type, height, shift, skip_width);
  }
}

template<int width, int height>
static void do_forward_vertical_pass(const int16_t* __restrict input, int16_t* __restrict output, tr_type_t type, int shift, int skip_height)
{
  if (NEON_FORCE_FAST_DCT2 && type == DCT2) {
    switch (height) {
      case  1: transform_forward_vertical<width, height>(input, output, type, shift, skip_height); break;
      case  2: transform_forward_vertical_Wx2_dct2(input, output, width, shift); break;
      case  4: transform_forward_vertical_Wx4_dct2(input, output, width, shift); break;
      case  8: transform_forward_vertical_Wx8_dct2(input, output, width, shift); break;
      case 16: transform_forward_vertical_Wx16_dct2(input, output, width, shift); break;
      case 32: transform_forward_vertical_Wx32_dct2(input, output, width, shift); break;
      default: NEON_UNREACHABLE("Unexpected height: %i\n", height);
    }
  }
  else if (NEON_PREFER_FAST_DCT2 && type == DCT2) {
    if constexpr (height == 16 && width >= 2) {
      transform_forward_vertical_Wx16_dct2(input, output, width, shift);
    }
    else if constexpr (height == 32 && width >= 4) {
      transform_forward_vertical_Wx32_dct2(input, output, width, shift);
    }
    else {
      transform_forward_vertical<width, height>(input, output, type, shift, skip_height);
    }
  }
  else {
    transform_forward_vertical<width, height>(input, output, type, shift, skip_height);
  }
}

// ---------------------------------------------------------
// MTS
// ---------------------------------------------------------

template<int width, int height>
static void mts_forward(
  const int16_t* __restrict input,
  int16_t* __restrict output,
  const cu_info_t* __restrict tu,
  color_t color,
  int bitdepth,
  bool is_isp,
  bool is_lfnst,
  bool explicit_mts,
  bool implicit_mts)
{
  if constexpr (width == 1) {
    const tr_type_t vertical_type = get_vertical_transform_type(tu, height, color, is_isp, is_lfnst, explicit_mts, implicit_mts);
    const int vertical_skip = get_vertical_skip(width, height, vertical_type, is_lfnst);
    const int vertical_shift = neon_comptime_log2(height) + bitdepth - 9;
    do_forward_vertical_pass<width, height>(input, output, vertical_type, vertical_shift, vertical_skip);
  }
  else if constexpr (height == 1) {
    const tr_type_t horizontal_type = get_horizontal_transform_type(tu, width, color, is_isp, is_lfnst, explicit_mts, implicit_mts);
    const int horizontal_skip = get_horizontal_skip(width, height, horizontal_type, is_lfnst);
    const int horizontal_shift = neon_comptime_log2(width) + bitdepth - 9;
    do_forward_horizontal_pass<width, height>(input, output, horizontal_type, horizontal_shift, horizontal_skip);
  }
  else {
    const tr_type_t horizontal_type = width <= 2 ? DCT2 : get_horizontal_transform_type(tu, width, color, is_isp, is_lfnst, explicit_mts, implicit_mts);
    const tr_type_t vertical_type = height <= 2 ? DCT2 : get_vertical_transform_type(tu, height, color, is_isp, is_lfnst, explicit_mts, implicit_mts);
    const int horizontal_skip = width <= 4 ? 0 : get_horizontal_skip(width, height, horizontal_type, is_lfnst);
    const int vertical_skip = height <= 4 ? 0 : get_vertical_skip(width, height, vertical_type, is_lfnst);
    const int horizontal_shift = neon_comptime_log2(width) + bitdepth - 9;
    const int vertical_shift = neon_comptime_log2(height) + 6;
    alignas(16) int16_t tmp[TR_MAX_WIDTH * TR_MAX_WIDTH];
    do_forward_horizontal_pass<width, height>(input, tmp, horizontal_type, horizontal_shift, horizontal_skip);
    do_forward_vertical_pass<width, height>(tmp, output, vertical_type, vertical_shift, vertical_skip);
  }
}

// ---------------------------------------------------------
// API
// ---------------------------------------------------------

extern "C" void uvg_mts_dct_neon(
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
          case 16: mts_forward<1, 16>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 32: mts_forward<1, 32>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          default: NEON_UNREACHABLE("Unexpected transform size: %ix%i\n", width, height);
        }
      }
      break;
    case 2: {
        switch (height) {
          case  2: mts_forward<2,  2>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  4: mts_forward<2,  4>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  8: mts_forward<2,  8>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 16: mts_forward<2, 16>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 32: mts_forward<2, 32>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          default: NEON_UNREACHABLE("Unexpected transform size: %ix%i\n", width, height);
        }
      }
      break;
    case 4: {
        switch (height) {
          case  2: mts_forward<4,  2>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  4: mts_forward<4,  4>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  8: mts_forward<4,  8>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 16: mts_forward<4, 16>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 32: mts_forward<4, 32>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          default: NEON_UNREACHABLE("Unexpected transform size: %ix%i\n", width, height);
        }
      }
      break;
    case 8: {
        switch (height) {
          case  2: mts_forward<8,  2>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  4: mts_forward<8,  4>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  8: mts_forward<8,  8>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 16: mts_forward<8, 16>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 32: mts_forward<8, 32>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          default: NEON_UNREACHABLE("Unexpected transform size: %ix%i\n", width, height);
        }
      }
      break;
    case 16: {
        switch (height) {
          case  1: mts_forward<16,  1>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  2: mts_forward<16,  2>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  4: mts_forward<16,  4>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  8: mts_forward<16,  8>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 16: mts_forward<16, 16>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 32: mts_forward<16, 32>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          default: NEON_UNREACHABLE("Unexpected transform size: %ix%i\n", width, height);
        }
      }
      break;
    case 32: {
        switch (height) {
          case  1: mts_forward<32,  1>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  2: mts_forward<32,  2>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  4: mts_forward<32,  4>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case  8: mts_forward<32,  8>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 16: mts_forward<32, 16>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          case 32: mts_forward<32, 32>(input, output, tu, color, bitdepth, is_isp, is_lfnst, explicit_mts, implicit_mts); break;
          default: NEON_UNREACHABLE("Unexpected transform size: %ix%i\n", width, height);
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
static inline void legacy_dct_NxN(const int16_t* __restrict input, int16_t* __restrict output, int bitdepth)
{
  const int horizontal_shift = neon_comptime_log2(N) + bitdepth - 9;
  const int vertical_shift = neon_comptime_log2(N) + 6;
  const int horizontal_skip = 0;
  const int vertical_skip = 0;

  alignas(16) int16_t tmp[N * N];
  do_forward_horizontal_pass<N, N>(input, tmp, DCT2, horizontal_shift, horizontal_skip);
  do_forward_vertical_pass<N, N>(tmp, output, DCT2, vertical_shift, vertical_skip);
}

extern "C" void uvg_dct_4x4_neon(int8_t bitdepth, const int16_t* input, int16_t* output)
{
  legacy_dct_NxN<4>(input, output, bitdepth);
};

extern "C" void uvg_dct_8x8_neon(int8_t bitdepth, const int16_t* input, int16_t* output)
{
  legacy_dct_NxN<8>(input, output, bitdepth);
};

extern "C" void uvg_dct_16x16_neon(int8_t bitdepth, const int16_t* input, int16_t* output)
{
  legacy_dct_NxN<16>(input, output, bitdepth);
};

extern "C" void uvg_dct_32x32_neon(int8_t bitdepth, const int16_t* input, int16_t* output)
{
  legacy_dct_NxN<32>(input, output, bitdepth);
};
