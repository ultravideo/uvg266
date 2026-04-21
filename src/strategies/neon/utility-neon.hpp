#ifndef STRATEGIES_UTILITY_NEON_HPP_
#define STRATEGIES_UTILITY_NEON_HPP_
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
 * Arm Neon optimizations.
 */

#include <arm_neon.h>
#include <cassert>

static constexpr int8_t neon_log2_table[65] = {
  -1,                             // 0
   0,  1, -1,  2, -1, -1, -1,     // 1-7
   3, -1, -1, -1, -1, -1, -1, -1, // 8-15
   4, -1, -1, -1, -1, -1, -1, -1, // 16-23
  -1, -1, -1, -1, -1, -1, -1, -1, // 24-31
   5, -1, -1, -1, -1, -1, -1, -1, // 32-39
  -1, -1, -1, -1, -1, -1, -1, -1, // 40-47
  -1, -1, -1, -1, -1, -1, -1, -1, // 48-55
  -1, -1, -1, -1, -1, -1, -1, -1, // 56-63
   6,                             // 64
};

static inline auto neon_fast_log2(int64_t value) -> int64_t
{
#if defined(__has_builtin) && __has_builtin(__builtin_clzll)
  constexpr int64_t total_bits = sizeof(unsigned long long) * 8;
  const int64_t leading_zeros = __builtin_clzll((unsigned long long)value);
  return (total_bits - 1) - leading_zeros;
#else
  return neon_log2_table[value];
#endif
}

static constexpr auto neon_comptime_log2(int64_t value) -> int64_t
{
  return neon_log2_table[value];
}

#endif // STRATEGIES_UTILITY_NEON_HPP_
