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
#include "transpose-neon.h"
}

#include "config-neon.hpp"
#include "debug-neon.hpp"

#include <arm_neon.h>
#include <cassert>
#include <cstdint>

// ---------------------------------------------------------
// Generic implementation
// ---------------------------------------------------------

// Unoptimized strided transpose.
template<typename T>
static inline auto transpose_naive(
  const T* __restrict in,
  T* __restrict out,
  int_fast32_t width,
  int_fast32_t height,
  int_fast32_t in_stride,
  int_fast32_t out_stride
) -> void
{
  for (int_fast32_t y = 0; y < height; ++y) {
    for (int_fast32_t x = 0; x < width; ++x) {
      out[y + x * out_stride] = in[x + y * in_stride];
    }
  }
}

// Unoptimized strided transpose.
template<typename T, int_fast32_t width, int_fast32_t height>
static inline auto transpose_naive(
  const T* __restrict in,
  T* __restrict out,
  int_fast32_t in_stride,
  int_fast32_t out_stride
) -> void
{
  for (int_fast32_t y = 0; y < height; ++y) {
    for (int_fast32_t x = 0; x < width; ++x) {
      out[y + x * out_stride] = in[x + y * in_stride];
    }
  }
}

// Unoptimized non-strided transpose.
template<typename T>
static inline auto transpose_naive(
  const T* __restrict in,
  T* __restrict out,
  int_fast32_t width,
  int_fast32_t height
) -> void
{
  for (int_fast32_t y = 0; y < height; ++y) {
    for (int_fast32_t x = 0; x < width; ++x) {
      out[y + x * height] = in[x + y * width];
    }
  }
}

// Unoptimized non-strided transpose.
template<typename T, int_fast32_t width, int_fast32_t height>
static inline auto transpose_naive(
  const T* __restrict in,
  T* __restrict out
) -> void
{
  for (int_fast32_t y = 0; y < height; ++y) {
    for (int_fast32_t x = 0; x < width; ++x) {
      out[y + x * height] = in[x + y * width];
    }
  }
}

// Default implementation of strided transpose.
// The default behavior is to keep splitting the matrix into smaller sub-matrices.
template<typename T, int W, int H>
static auto transpose(
  const T* __restrict in,
  T* __restrict out,
  int_fast32_t in_stride,
  int_fast32_t out_stride
) -> void
{
  if constexpr (W == 0 || H == 0) {
    static_assert(sizeof(T) != sizeof(T), "Never-ending recursion: one or more specializations are missing");
  }
  else if constexpr (W > H) {
    // Wide matrix -> transpose as two smaller blocks:
    // |A B|' <=> |A'|
    //            |B'|
    constexpr int_fast32_t w = W / 2;
    constexpr int_fast32_t h = H;
    transpose<T, w, h>(in + 0 * in_stride + 0, out + 0 * out_stride + 0, in_stride, out_stride);
    transpose<T, w, h>(in + 0 * in_stride + w, out + w * out_stride + 0, in_stride, out_stride);
  }
  else if constexpr (W < H) {
    // Tall matrix -> transpose as two smaller blocks:
    // |A|' <=> |A' B'|
    // |B|
    constexpr int_fast32_t w = W;
    constexpr int_fast32_t h = H / 2;
    transpose<T, w, h>(in + 0 * in_stride + 0, out + 0 * out_stride + 0, in_stride, out_stride);
    transpose<T, w, h>(in + h * in_stride + 0, out + 0 * out_stride + h, in_stride, out_stride);
  }
  else {
    // Square matrix -> transpose as four smaller blocks:
    // |A B|' <=> |A' C'|
    // |C D|      |B' D'|
    constexpr int_fast32_t w = W / 2;
    constexpr int_fast32_t h = H / 2;
    transpose<T, w, h>(in + 0 * in_stride + 0, out + 0 * out_stride + 0, in_stride, out_stride);
    transpose<T, w, h>(in + h * in_stride + 0, out + 0 * out_stride + w, in_stride, out_stride);
    transpose<T, w, h>(in + 0 * in_stride + w, out + h * out_stride + 0, in_stride, out_stride);
    transpose<T, w, h>(in + h * in_stride + w, out + h * out_stride + w, in_stride, out_stride);
  }
}

// Default implementation of non-strided transpose.
// The default behavior is to delegate to the strided version.
template<typename T, int W, int H>
static inline auto transpose(
  const T* __restrict in,
  T* __restrict out
) -> void
{
  if constexpr (W == 1 || H == 1) {
    // Special case: 1xN and Nx1 matrices that are not strided can simply be copied
    // because the memory layout remains the same.
    memcpy(out, in, W * H * sizeof(T));
  }
  else {
    // Delegate to the strided version.
    transpose<T, W, H>(in, out, W, H);
  }
}

// ---------------------------------------------------------
// 8-bit specializations
// ---------------------------------------------------------

template<>
inline auto transpose<uint8_t, 2, 2>(
  const uint8_t* __restrict in,
  uint8_t* __restrict out,
  int_fast32_t in_stride,
  int_fast32_t out_stride
) -> void
{
  out[0 * out_stride + 0] = in[0 * in_stride + 0];
  out[0 * out_stride + 1] = in[1 * in_stride + 0];
  out[1 * out_stride + 0] = in[0 * in_stride + 1];
  out[1 * out_stride + 1] = in[1 * in_stride + 1];
}

template<>
inline auto transpose<uint8_t, 2, 8>(
  const uint8_t* __restrict in,
  uint8_t* __restrict out,
  int_fast32_t in_stride,
  int_fast32_t out_stride
) -> void
{
  if constexpr (!NEON_ENABLE_INTRINSICS) {
    transpose_naive<uint8_t, 2, 8>(in, out, in_stride, out_stride);
  }
  else {
    alignas(16) uint8_t tmp[2 * 8];
    memcpy(tmp +  0, in + 0 * in_stride, 2 * sizeof(uint8_t));
    memcpy(tmp +  2, in + 1 * in_stride, 2 * sizeof(uint8_t));
    memcpy(tmp +  4, in + 2 * in_stride, 2 * sizeof(uint8_t));
    memcpy(tmp +  6, in + 3 * in_stride, 2 * sizeof(uint8_t));
    memcpy(tmp +  8, in + 4 * in_stride, 2 * sizeof(uint8_t));
    memcpy(tmp + 10, in + 5 * in_stride, 2 * sizeof(uint8_t));
    memcpy(tmp + 12, in + 6 * in_stride, 2 * sizeof(uint8_t));
    memcpy(tmp + 14, in + 7 * in_stride, 2 * sizeof(uint8_t));
    const uint8x16_t vin = vld1q_u8(tmp);
    const uint8x16_t vshuffle = {0, 2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 9, 11, 13, 15};
    const uint8x16_t vout = vqtbl1q_u8(vin, vshuffle);
    vst1q_u8(tmp, vout);
    memcpy(out + 0 * out_stride, tmp + 0, 8 * sizeof(uint8_t));
    memcpy(out + 1 * out_stride, tmp + 8, 8 * sizeof(uint8_t));
  }
}

template<>
inline auto transpose<uint8_t, 4, 4>(
  const uint8_t* __restrict in,
  uint8_t* __restrict out,
  int_fast32_t in_stride,
  int_fast32_t out_stride
) -> void
{
  if constexpr (!NEON_ENABLE_INTRINSICS) {
    transpose_naive<uint8_t, 4, 4>(in, out, in_stride, out_stride);
  }
  else {
    alignas(16) uint8_t tmp[4 * 4];
    memcpy(tmp +  0, in + 0 * in_stride, 4 * sizeof(uint8_t));
    memcpy(tmp +  4, in + 1 * in_stride, 4 * sizeof(uint8_t));
    memcpy(tmp +  8, in + 2 * in_stride, 4 * sizeof(uint8_t));
    memcpy(tmp + 12, in + 3 * in_stride, 4 * sizeof(uint8_t));
    const uint8x16_t vin = vld1q_u8(tmp);
    const uint8x16_t vshuffle = {0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15};
    const uint8x16_t vout = vqtbl1q_u8(vin, vshuffle);
    vst1q_u8(tmp, vout);
    memcpy(out + 0 * out_stride, tmp +  0, 4 * sizeof(uint8_t));
    memcpy(out + 1 * out_stride, tmp +  4, 4 * sizeof(uint8_t));
    memcpy(out + 2 * out_stride, tmp +  8, 4 * sizeof(uint8_t));
    memcpy(out + 3 * out_stride, tmp + 12, 4 * sizeof(uint8_t));
  }
}

template<>
inline auto transpose<uint8_t, 4, 4>(
  const uint8_t* __restrict in,
  uint8_t* __restrict out
) -> void
{
  if constexpr (!NEON_ENABLE_INTRINSICS) {
    transpose_naive<uint8_t, 4, 4>(in, out);
  }
  else {
    const uint8x16_t vshuffle = {0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15};
    const uint8x16_t vin = vld1q_u8(in);
    const uint8x16_t vout = vqtbl1q_u8(vin, vshuffle);
    vst1q_u8(out, vout);
  }
}

template<>
inline auto transpose<uint8_t, 8, 2>(
  const uint8_t* __restrict in,
  uint8_t* __restrict out,
  int_fast32_t in_stride,
  int_fast32_t out_stride
) -> void
{
  if constexpr (!NEON_ENABLE_INTRINSICS) {
    transpose_naive<uint8_t, 8, 2>(in, out, in_stride, out_stride);
  }
  else {
    alignas(16) uint8_t tmp[8 * 2];
    memcpy(tmp + 0, in + 0 * in_stride, 8 * sizeof(uint8_t));
    memcpy(tmp + 8, in + 1 * in_stride, 8 * sizeof(uint8_t));
    const uint8x16_t vin = vld1q_u8(tmp);
    const uint8x16_t vshuffle = {0, 8, 1, 9, 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15};
    const uint8x16_t vout = vqtbl1q_u8(vin, vshuffle);
    vst1q_u8(tmp, vout);
    memcpy(out + 0 * out_stride, tmp +  0, 2 * sizeof(uint8_t));
    memcpy(out + 1 * out_stride, tmp +  2, 2 * sizeof(uint8_t));
    memcpy(out + 2 * out_stride, tmp +  4, 2 * sizeof(uint8_t));
    memcpy(out + 3 * out_stride, tmp +  6, 2 * sizeof(uint8_t));
    memcpy(out + 4 * out_stride, tmp +  8, 2 * sizeof(uint8_t));
    memcpy(out + 5 * out_stride, tmp + 10, 2 * sizeof(uint8_t));
    memcpy(out + 6 * out_stride, tmp + 12, 2 * sizeof(uint8_t));
    memcpy(out + 7 * out_stride, tmp + 14, 2 * sizeof(uint8_t));
  }
}

template<>
inline auto transpose<uint8_t, 8, 8>(
  const uint8_t* __restrict in,
  uint8_t* __restrict out,
  int_fast32_t in_stride,
  int_fast32_t out_stride
) -> void
{
  if constexpr (!NEON_ENABLE_INTRINSICS) {
    transpose_naive<uint8_t, 8, 8>(in, out, in_stride, out_stride);
  }
  else {
    const auto in_00 = vld1_u8(in + 0 * in_stride);
    const auto in_01 = vld1_u8(in + 1 * in_stride);
    const auto in_02 = vld1_u8(in + 2 * in_stride);
    const auto in_03 = vld1_u8(in + 3 * in_stride);
    const auto in_04 = vld1_u8(in + 4 * in_stride);
    const auto in_05 = vld1_u8(in + 5 * in_stride);
    const auto in_06 = vld1_u8(in + 6 * in_stride);
    const auto in_07 = vld1_u8(in + 7 * in_stride);

    const auto [a_00, a_01] = vtrn_u8(in_00, in_01).val;
    const auto [a_02, a_03] = vtrn_u8(in_02, in_03).val;
    const auto [a_04, a_05] = vtrn_u8(in_04, in_05).val;
    const auto [a_06, a_07] = vtrn_u8(in_06, in_07).val;

    const auto [b_00, b_02] = vtrn_u16(vreinterpret_u16_u8(a_00), vreinterpret_u16_u8(a_02)).val;
    const auto [b_01, b_03] = vtrn_u16(vreinterpret_u16_u8(a_01), vreinterpret_u16_u8(a_03)).val;
    const auto [b_04, b_06] = vtrn_u16(vreinterpret_u16_u8(a_04), vreinterpret_u16_u8(a_06)).val;
    const auto [b_05, b_07] = vtrn_u16(vreinterpret_u16_u8(a_05), vreinterpret_u16_u8(a_07)).val;

    const auto [c_00, c_04] = vtrn_u32(vreinterpret_u32_u16(b_00), vreinterpret_u32_u16(b_04)).val;
    const auto [c_01, c_05] = vtrn_u32(vreinterpret_u32_u16(b_01), vreinterpret_u32_u16(b_05)).val;
    const auto [c_02, c_06] = vtrn_u32(vreinterpret_u32_u16(b_02), vreinterpret_u32_u16(b_06)).val;
    const auto [c_03, c_07] = vtrn_u32(vreinterpret_u32_u16(b_03), vreinterpret_u32_u16(b_07)).val;

    vst1_u8(out + 0 * out_stride, vreinterpret_u8_u32(c_00));
    vst1_u8(out + 1 * out_stride, vreinterpret_u8_u32(c_01));
    vst1_u8(out + 2 * out_stride, vreinterpret_u8_u32(c_02));
    vst1_u8(out + 3 * out_stride, vreinterpret_u8_u32(c_03));
    vst1_u8(out + 4 * out_stride, vreinterpret_u8_u32(c_04));
    vst1_u8(out + 5 * out_stride, vreinterpret_u8_u32(c_05));
    vst1_u8(out + 6 * out_stride, vreinterpret_u8_u32(c_06));
    vst1_u8(out + 7 * out_stride, vreinterpret_u8_u32(c_07));
  }
}

template<>
inline auto transpose<uint8_t, 16, 8>(
  const uint8_t* __restrict in,
  uint8_t* __restrict out,
  int_fast32_t in_stride,
  int_fast32_t out_stride
) -> void
{
  if constexpr (!NEON_ENABLE_INTRINSICS) {
    transpose_naive<uint8_t, 16, 8>(in, out, in_stride, out_stride);
  }
  else {
    const auto in_00 = vld1q_u8(in + 0 * in_stride);
    const auto in_01 = vld1q_u8(in + 1 * in_stride);
    const auto in_02 = vld1q_u8(in + 2 * in_stride);
    const auto in_03 = vld1q_u8(in + 3 * in_stride);
    const auto in_04 = vld1q_u8(in + 4 * in_stride);
    const auto in_05 = vld1q_u8(in + 5 * in_stride);
    const auto in_06 = vld1q_u8(in + 6 * in_stride);
    const auto in_07 = vld1q_u8(in + 7 * in_stride);

    const auto [a_00, a_01] = vtrnq_u8(in_00, in_01).val;
    const auto [a_02, a_03] = vtrnq_u8(in_02, in_03).val;
    const auto [a_04, a_05] = vtrnq_u8(in_04, in_05).val;
    const auto [a_06, a_07] = vtrnq_u8(in_06, in_07).val;

    const auto [b_00, b_02] = vtrnq_u16(vreinterpretq_u16_u8(a_00), vreinterpretq_u16_u8(a_02)).val;
    const auto [b_01, b_03] = vtrnq_u16(vreinterpretq_u16_u8(a_01), vreinterpretq_u16_u8(a_03)).val;
    const auto [b_04, b_06] = vtrnq_u16(vreinterpretq_u16_u8(a_04), vreinterpretq_u16_u8(a_06)).val;
    const auto [b_05, b_07] = vtrnq_u16(vreinterpretq_u16_u8(a_05), vreinterpretq_u16_u8(a_07)).val;

    const auto [c_00, c_04] = vtrnq_u32(vreinterpretq_u32_u16(b_00), vreinterpretq_u32_u16((b_04))).val;
    const auto [c_01, c_05] = vtrnq_u32(vreinterpretq_u32_u16(b_01), vreinterpretq_u32_u16((b_05))).val;
    const auto [c_02, c_06] = vtrnq_u32(vreinterpretq_u32_u16(b_02), vreinterpretq_u32_u16((b_06))).val;
    const auto [c_03, c_07] = vtrnq_u32(vreinterpretq_u32_u16(b_03), vreinterpretq_u32_u16((b_07))).val;

    vst1_u8(out + 0  * out_stride, vget_low_u8(vreinterpretq_u8_u32(c_00)));
    vst1_u8(out + 1  * out_stride, vget_low_u8(vreinterpretq_u8_u32(c_01)));
    vst1_u8(out + 2  * out_stride, vget_low_u8(vreinterpretq_u8_u32(c_02)));
    vst1_u8(out + 3  * out_stride, vget_low_u8(vreinterpretq_u8_u32(c_03)));
    vst1_u8(out + 4  * out_stride, vget_low_u8(vreinterpretq_u8_u32(c_04)));
    vst1_u8(out + 5  * out_stride, vget_low_u8(vreinterpretq_u8_u32(c_05)));
    vst1_u8(out + 6  * out_stride, vget_low_u8(vreinterpretq_u8_u32(c_06)));
    vst1_u8(out + 7  * out_stride, vget_low_u8(vreinterpretq_u8_u32(c_07)));
    vst1_u8(out + 8  * out_stride, vget_high_u8(vreinterpretq_u8_u32(c_00)));
    vst1_u8(out + 9  * out_stride, vget_high_u8(vreinterpretq_u8_u32(c_01)));
    vst1_u8(out + 10 * out_stride, vget_high_u8(vreinterpretq_u8_u32(c_02)));
    vst1_u8(out + 11 * out_stride, vget_high_u8(vreinterpretq_u8_u32(c_03)));
    vst1_u8(out + 12 * out_stride, vget_high_u8(vreinterpretq_u8_u32(c_04)));
    vst1_u8(out + 13 * out_stride, vget_high_u8(vreinterpretq_u8_u32(c_05)));
    vst1_u8(out + 14 * out_stride, vget_high_u8(vreinterpretq_u8_u32(c_06)));
    vst1_u8(out + 15 * out_stride, vget_high_u8(vreinterpretq_u8_u32(c_07)));
  }
}

template<>
inline auto transpose<uint8_t, 16, 16>(
  const uint8_t* __restrict in,
  uint8_t* __restrict out,
  int_fast32_t in_stride,
  int_fast32_t out_stride
) -> void
{
  if constexpr (!NEON_ENABLE_INTRINSICS) {
    transpose_naive<uint8_t, 16, 16>(in, out, in_stride, out_stride);
  }
  else {
    const auto in_00 = vld1q_u8(in +  0 * in_stride);
    const auto in_01 = vld1q_u8(in +  1 * in_stride);
    const auto in_02 = vld1q_u8(in +  2 * in_stride);
    const auto in_03 = vld1q_u8(in +  3 * in_stride);
    const auto in_04 = vld1q_u8(in +  4 * in_stride);
    const auto in_05 = vld1q_u8(in +  5 * in_stride);
    const auto in_06 = vld1q_u8(in +  6 * in_stride);
    const auto in_07 = vld1q_u8(in +  7 * in_stride);
    const auto in_08 = vld1q_u8(in +  8 * in_stride);
    const auto in_09 = vld1q_u8(in +  9 * in_stride);
    const auto in_10 = vld1q_u8(in + 10 * in_stride);
    const auto in_11 = vld1q_u8(in + 11 * in_stride);
    const auto in_12 = vld1q_u8(in + 12 * in_stride);
    const auto in_13 = vld1q_u8(in + 13 * in_stride);
    const auto in_14 = vld1q_u8(in + 14 * in_stride);
    const auto in_15 = vld1q_u8(in + 15 * in_stride);

    const auto [a_00, a_01] = vtrnq_u8(in_00, in_01).val;
    const auto [a_02, a_03] = vtrnq_u8(in_02, in_03).val;
    const auto [a_04, a_05] = vtrnq_u8(in_04, in_05).val;
    const auto [a_06, a_07] = vtrnq_u8(in_06, in_07).val;
    const auto [a_08, a_09] = vtrnq_u8(in_08, in_09).val;
    const auto [a_10, a_11] = vtrnq_u8(in_10, in_11).val;
    const auto [a_12, a_13] = vtrnq_u8(in_12, in_13).val;
    const auto [a_14, a_15] = vtrnq_u8(in_14, in_15).val;

    const auto [b_00, b_02] = vtrnq_u16(vreinterpretq_u16_u8(a_00), vreinterpretq_u16_u8(a_02)).val;
    const auto [b_01, b_03] = vtrnq_u16(vreinterpretq_u16_u8(a_01), vreinterpretq_u16_u8(a_03)).val;
    const auto [b_04, b_06] = vtrnq_u16(vreinterpretq_u16_u8(a_04), vreinterpretq_u16_u8(a_06)).val;
    const auto [b_05, b_07] = vtrnq_u16(vreinterpretq_u16_u8(a_05), vreinterpretq_u16_u8(a_07)).val;
    const auto [b_08, b_10] = vtrnq_u16(vreinterpretq_u16_u8(a_08), vreinterpretq_u16_u8(a_10)).val;
    const auto [b_09, b_11] = vtrnq_u16(vreinterpretq_u16_u8(a_09), vreinterpretq_u16_u8(a_11)).val;
    const auto [b_12, b_14] = vtrnq_u16(vreinterpretq_u16_u8(a_12), vreinterpretq_u16_u8(a_14)).val;
    const auto [b_13, b_15] = vtrnq_u16(vreinterpretq_u16_u8(a_13), vreinterpretq_u16_u8(a_15)).val;

    const auto [c_00, c_04] = vtrnq_u32(vreinterpretq_u32_u16(b_00), vreinterpretq_u32_u16(b_04)).val;
    const auto [c_01, c_05] = vtrnq_u32(vreinterpretq_u32_u16(b_01), vreinterpretq_u32_u16(b_05)).val;
    const auto [c_02, c_06] = vtrnq_u32(vreinterpretq_u32_u16(b_02), vreinterpretq_u32_u16(b_06)).val;
    const auto [c_03, c_07] = vtrnq_u32(vreinterpretq_u32_u16(b_03), vreinterpretq_u32_u16(b_07)).val;
    const auto [c_08, c_12] = vtrnq_u32(vreinterpretq_u32_u16(b_08), vreinterpretq_u32_u16(b_12)).val;
    const auto [c_09, c_13] = vtrnq_u32(vreinterpretq_u32_u16(b_09), vreinterpretq_u32_u16(b_13)).val;
    const auto [c_10, c_14] = vtrnq_u32(vreinterpretq_u32_u16(b_10), vreinterpretq_u32_u16(b_14)).val;
    const auto [c_11, c_15] = vtrnq_u32(vreinterpretq_u32_u16(b_11), vreinterpretq_u32_u16(b_15)).val;

    const auto d_00 = vzip1q_u64(vreinterpretq_u64_u32(c_00), vreinterpretq_u64_u32(c_08));
    const auto d_01 = vzip1q_u64(vreinterpretq_u64_u32(c_01), vreinterpretq_u64_u32(c_09));
    const auto d_02 = vzip1q_u64(vreinterpretq_u64_u32(c_02), vreinterpretq_u64_u32(c_10));
    const auto d_03 = vzip1q_u64(vreinterpretq_u64_u32(c_03), vreinterpretq_u64_u32(c_11));
    const auto d_04 = vzip1q_u64(vreinterpretq_u64_u32(c_04), vreinterpretq_u64_u32(c_12));
    const auto d_05 = vzip1q_u64(vreinterpretq_u64_u32(c_05), vreinterpretq_u64_u32(c_13));
    const auto d_06 = vzip1q_u64(vreinterpretq_u64_u32(c_06), vreinterpretq_u64_u32(c_14));
    const auto d_07 = vzip1q_u64(vreinterpretq_u64_u32(c_07), vreinterpretq_u64_u32(c_15));
    const auto d_08 = vzip2q_u64(vreinterpretq_u64_u32(c_00), vreinterpretq_u64_u32(c_08));
    const auto d_09 = vzip2q_u64(vreinterpretq_u64_u32(c_01), vreinterpretq_u64_u32(c_09));
    const auto d_10 = vzip2q_u64(vreinterpretq_u64_u32(c_02), vreinterpretq_u64_u32(c_10));
    const auto d_11 = vzip2q_u64(vreinterpretq_u64_u32(c_03), vreinterpretq_u64_u32(c_11));
    const auto d_12 = vzip2q_u64(vreinterpretq_u64_u32(c_04), vreinterpretq_u64_u32(c_12));
    const auto d_13 = vzip2q_u64(vreinterpretq_u64_u32(c_05), vreinterpretq_u64_u32(c_13));
    const auto d_14 = vzip2q_u64(vreinterpretq_u64_u32(c_06), vreinterpretq_u64_u32(c_14));
    const auto d_16 = vzip2q_u64(vreinterpretq_u64_u32(c_07), vreinterpretq_u64_u32(c_15));

    vst1q_u8(out +  0 * out_stride, vreinterpretq_u8_u64(d_00));
    vst1q_u8(out +  1 * out_stride, vreinterpretq_u8_u64(d_01));
    vst1q_u8(out +  2 * out_stride, vreinterpretq_u8_u64(d_02));
    vst1q_u8(out +  3 * out_stride, vreinterpretq_u8_u64(d_03));
    vst1q_u8(out +  4 * out_stride, vreinterpretq_u8_u64(d_04));
    vst1q_u8(out +  5 * out_stride, vreinterpretq_u8_u64(d_05));
    vst1q_u8(out +  6 * out_stride, vreinterpretq_u8_u64(d_06));
    vst1q_u8(out +  7 * out_stride, vreinterpretq_u8_u64(d_07));
    vst1q_u8(out +  8 * out_stride, vreinterpretq_u8_u64(d_08));
    vst1q_u8(out +  9 * out_stride, vreinterpretq_u8_u64(d_09));
    vst1q_u8(out + 10 * out_stride, vreinterpretq_u8_u64(d_10));
    vst1q_u8(out + 11 * out_stride, vreinterpretq_u8_u64(d_11));
    vst1q_u8(out + 12 * out_stride, vreinterpretq_u8_u64(d_12));
    vst1q_u8(out + 13 * out_stride, vreinterpretq_u8_u64(d_13));
    vst1q_u8(out + 14 * out_stride, vreinterpretq_u8_u64(d_14));
    vst1q_u8(out + 15 * out_stride, vreinterpretq_u8_u64(d_16));
  }
}

// ---------------------------------------------------------
// 16-bit specializations
// ---------------------------------------------------------

template<>
inline auto transpose<uint16_t, 2, 2>(
  const uint16_t* __restrict in,
  uint16_t* __restrict out,
  int_fast32_t in_stride,
  int_fast32_t out_stride
) -> void
{
  out[0 * out_stride + 0] = in[0 * in_stride + 0];
  out[0 * out_stride + 1] = in[1 * in_stride + 0];
  out[1 * out_stride + 0] = in[0 * in_stride + 1];
  out[1 * out_stride + 1] = in[1 * in_stride + 1];
}

template<>
inline auto transpose<uint16_t, 4, 4>(
  const uint16_t* __restrict in,
  uint16_t* __restrict out,
  int_fast32_t in_stride,
  int_fast32_t out_stride
) -> void
{
  if constexpr (!NEON_ENABLE_INTRINSICS) {
    transpose_naive<uint16_t, 4, 4>(in, out, in_stride, out_stride);
  }
  else {
    const auto in_00 = vld1_u16(in + 0 * in_stride);
    const auto in_01 = vld1_u16(in + 1 * in_stride);
    const auto in_02 = vld1_u16(in + 2 * in_stride);
    const auto in_03 = vld1_u16(in + 3 * in_stride);

    const auto [a_00, a_01] = vtrn_u16(in_00, in_01).val;
    const auto [a_02, a_03] = vtrn_u16(in_02, in_03).val;

    const auto [b_00, b_02] = vtrn_u32(vreinterpret_u32_u16(a_00), vreinterpret_u32_u16(a_02)).val;
    const auto [b_01, b_03] = vtrn_u32(vreinterpret_u32_u16(a_01), vreinterpret_u32_u16(a_03)).val;

    vst1_u16(out + 0 * out_stride, vreinterpret_u16_u32(b_00));
    vst1_u16(out + 1 * out_stride, vreinterpret_u16_u32(b_01));
    vst1_u16(out + 2 * out_stride, vreinterpret_u16_u32(b_02));
    vst1_u16(out + 3 * out_stride, vreinterpret_u16_u32(b_03));
  }
}

template<>
inline auto transpose<uint16_t, 4, 4>(
  const uint16_t* __restrict in,
  uint16_t* __restrict out
) -> void
{
  if constexpr (!NEON_ENABLE_INTRINSICS) {
    transpose_naive<uint16_t, 4, 4>(in, out);
  }
  else {
    vst1_u16_x4(out, vld4_u16(in));
  }
}

template<>
inline auto transpose<uint16_t, 8, 8>(
  const uint16_t* __restrict in,
  uint16_t* __restrict out,
  int_fast32_t in_stride,
  int_fast32_t out_stride
) -> void
{
  if constexpr (!NEON_ENABLE_INTRINSICS) {
    transpose_naive<uint16_t, 8, 8>(in, out, in_stride, out_stride);
  }
  else {
    const auto in_00 = vld1q_u16(in + 0 * in_stride);
    const auto in_01 = vld1q_u16(in + 1 * in_stride);
    const auto in_02 = vld1q_u16(in + 2 * in_stride);
    const auto in_03 = vld1q_u16(in + 3 * in_stride);
    const auto in_04 = vld1q_u16(in + 4 * in_stride);
    const auto in_05 = vld1q_u16(in + 5 * in_stride);
    const auto in_06 = vld1q_u16(in + 6 * in_stride);
    const auto in_07 = vld1q_u16(in + 7 * in_stride);

    const auto [a_00, a_01] = vtrnq_u16(in_00, in_01).val;
    const auto [a_02, a_03] = vtrnq_u16(in_02, in_03).val;
    const auto [a_04, a_05] = vtrnq_u16(in_04, in_05).val;
    const auto [a_06, a_07] = vtrnq_u16(in_06, in_07).val;

    const auto [b_00, b_02] = vtrnq_u32(vreinterpretq_u32_u16(a_00), vreinterpretq_u32_u16(a_02)).val;
    const auto [b_01, b_03] = vtrnq_u32(vreinterpretq_u32_u16(a_01), vreinterpretq_u32_u16(a_03)).val;
    const auto [b_04, b_06] = vtrnq_u32(vreinterpretq_u32_u16(a_04), vreinterpretq_u32_u16(a_06)).val;
    const auto [b_05, b_07] = vtrnq_u32(vreinterpretq_u32_u16(a_05), vreinterpretq_u32_u16(a_07)).val;

    const auto c_00 = vzip1q_u64(vreinterpretq_u64_u32(b_00), vreinterpretq_u64_u32(b_04));
    const auto c_01 = vzip1q_u64(vreinterpretq_u64_u32(b_01), vreinterpretq_u64_u32(b_05));
    const auto c_02 = vzip1q_u64(vreinterpretq_u64_u32(b_02), vreinterpretq_u64_u32(b_06));
    const auto c_03 = vzip1q_u64(vreinterpretq_u64_u32(b_03), vreinterpretq_u64_u32(b_07));
    const auto c_04 = vzip2q_u64(vreinterpretq_u64_u32(b_00), vreinterpretq_u64_u32(b_04));
    const auto c_05 = vzip2q_u64(vreinterpretq_u64_u32(b_01), vreinterpretq_u64_u32(b_05));
    const auto c_06 = vzip2q_u64(vreinterpretq_u64_u32(b_02), vreinterpretq_u64_u32(b_06));
    const auto c_07 = vzip2q_u64(vreinterpretq_u64_u32(b_03), vreinterpretq_u64_u32(b_07));

    vst1q_u16(out + 0 * out_stride, vreinterpretq_u16_u64(c_00));
    vst1q_u16(out + 1 * out_stride, vreinterpretq_u16_u64(c_01));
    vst1q_u16(out + 2 * out_stride, vreinterpretq_u16_u64(c_02));
    vst1q_u16(out + 3 * out_stride, vreinterpretq_u16_u64(c_03));
    vst1q_u16(out + 4 * out_stride, vreinterpretq_u16_u64(c_04));
    vst1q_u16(out + 5 * out_stride, vreinterpretq_u16_u64(c_05));
    vst1q_u16(out + 6 * out_stride, vreinterpretq_u16_u64(c_06));
    vst1q_u16(out + 7 * out_stride, vreinterpretq_u16_u64(c_07));
  }
}

// ---------------------------------------------------------
// Dispatch
// ---------------------------------------------------------

template<typename T>
static auto transpose(
  const void* __restrict src,
  void* __restrict dst,
  int width,
  int height
) -> void
{
  const T* in = (const T*)src;
  T* out = (T*)dst;

  switch (width) {
    case 1:
      switch (height) {
        case  1: transpose<T, 1,  1>(in, out); return;
        case  2: transpose<T, 1,  2>(in, out); return;
        case  4: transpose<T, 1,  4>(in, out); return;
        case  8: transpose<T, 1,  8>(in, out); return;
        case 16: transpose<T, 1, 16>(in, out); return;
        case 32: transpose<T, 1, 32>(in, out); return;
        case 64: transpose<T, 1, 64>(in, out); return;
      }
      break;
    case 2:
      switch (height) {
        case  1: transpose<T, 2,  1>(in, out); return;
        case  2: transpose<T, 2,  2>(in, out); return;
        case  4: transpose<T, 2,  4>(in, out); return;
        case  8: transpose<T, 2,  8>(in, out); return;
        case 16: transpose<T, 2, 16>(in, out); return;
        case 32: transpose<T, 2, 32>(in, out); return;
        case 64: transpose<T, 2, 64>(in, out); return;
      }
      break;
    case 4:
      switch (height) {
        case  1: transpose<T, 4,  1>(in, out); return;
        case  2: transpose<T, 4,  2>(in, out); return;
        case  4: transpose<T, 4,  4>(in, out); return;
        case  8: transpose<T, 4,  8>(in, out); return;
        case 16: transpose<T, 4, 16>(in, out); return;
        case 32: transpose<T, 4, 32>(in, out); return;
        case 64: transpose<T, 4, 64>(in, out); return;
      }
      break;
    case 8:
      switch (height) {
        case  1: transpose<T, 8,  1>(in, out); return;
        case  2: transpose<T, 8,  2>(in, out); return;
        case  4: transpose<T, 8,  4>(in, out); return;
        case  8: transpose<T, 8,  8>(in, out); return;
        case 16: transpose<T, 8, 16>(in, out); return;
        case 32: transpose<T, 8, 32>(in, out); return;
        case 64: transpose<T, 8, 64>(in, out); return;
      }
      break;
    case 16:
      switch (height) {
        case  1: transpose<T, 16,  1>(in, out); return;
        case  2: transpose<T, 16,  2>(in, out); return;
        case  4: transpose<T, 16,  4>(in, out); return;
        case  8: transpose<T, 16,  8>(in, out); return;
        case 16: transpose<T, 16, 16>(in, out); return;
        case 32: transpose<T, 16, 32>(in, out); return;
        case 64: transpose<T, 16, 64>(in, out); return;
      }
      break;
    case 32:
      switch (height) {
        case  1: transpose<T, 32,  1>(in, out); return;
        case  2: transpose<T, 32,  2>(in, out); return;
        case  4: transpose<T, 32,  4>(in, out); return;
        case  8: transpose<T, 32,  8>(in, out); return;
        case 16: transpose<T, 32, 16>(in, out); return;
        case 32: transpose<T, 32, 32>(in, out); return;
        case 64: transpose<T, 32, 64>(in, out); return;
      }
      break;
    case 64:
      switch (height) {
        case  1: transpose<T, 64,  1>(in, out); return;
        case  2: transpose<T, 64,  2>(in, out); return;
        case  4: transpose<T, 64,  4>(in, out); return;
        case  8: transpose<T, 64,  8>(in, out); return;
        case 16: transpose<T, 64, 16>(in, out); return;
        case 32: transpose<T, 64, 32>(in, out); return;
        case 64: transpose<T, 64, 64>(in, out); return;
      }
      break;
  }

  transpose_naive<T>(in, out, width, height);
}

// ---------------------------------------------------------
// API
// ---------------------------------------------------------

extern "C" void uvg_transpose_8bit_neon(
  const void* __restrict src,
  void* __restrict dst,
  int width,
  int height
)
{
  transpose<uint8_t>(src, dst, width, height);
}

extern "C" void uvg_transpose_16bit_neon(
  const void* __restrict src,
  void* __restrict dst,
  int width,
  int height
)
{
  transpose<uint16_t>(src, dst, width, height);
}

extern "C" void uvg_transpose_px_neon(
  const uvg_pixel* __restrict src,
  uvg_pixel* __restrict dst,
  int width,
  int height
)
{
  if constexpr (sizeof(uvg_pixel) == sizeof(uint8_t)) {
    transpose<uint8_t>(src, dst, width, height);
  }
  else if constexpr (sizeof(uvg_pixel) == sizeof(uint16_t)) {
    transpose<uint16_t>(src, dst, width, height);
  }
  else {
    NEON_UNREACHABLE("Unexpected sizeof(uvg_pixel): %lu\n", sizeof(uvg_pixel));
  }
}
