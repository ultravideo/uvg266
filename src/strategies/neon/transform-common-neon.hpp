#ifndef STRATEGIES_TRANSFORM_COMMON_NEON_HPP_
#define STRATEGIES_TRANSFORM_COMMON_NEON_HPP_
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

extern "C" {
#include "cu.h"
#include "global.h"
}

#include <arm_neon.h>

#include <cassert>
#include <climits>

inline int get_horizontal_skip(int width, int height, tr_type_t transform_type, bool is_lfnst)
{
  int skip_width = (transform_type != DCT2 && width == 32) ? 16 : (width > 32 ? width - 32 : 0);

  if (is_lfnst) {
    if ((width == 4 && height > 4) || (width > 4 && height == 4)) {
      skip_width = width - 4;
    }
    else if ((width >= 8 && height >= 8)) {
      skip_width = width - 8;
    }
  }

  return skip_width;
}

inline int get_vertical_skip(int width, int height, tr_type_t transform_type, bool is_lfnst)
{
  int skip_height = (transform_type != DCT2 && height == 32) ? 16 : (height > 32 ? height - 32 : 0);

  if (is_lfnst) {
    if ((width == 4 && height > 4) || (width > 4 && height == 4)) {
      skip_height = height - 4;
    }
    else if ((width >= 8 && height >= 8)) {
      skip_height = height - 8;
    }
  }

  return skip_height;
}

inline tr_type_t get_horizontal_transform_type(
  const cu_info_t* __restrict tu,
  int width,
  color_t color,
  bool is_isp,
  bool is_lfnst,
  bool explicit_mts,
  bool implicit_mts)
{
  if (color != COLOR_Y) {
    return DCT2;
  }

  if (is_isp && is_lfnst) {
    return DCT2;
  }

  if (implicit_mts || (is_isp && explicit_mts)) {
    if ((width >= 4) && (width <= 16)) {
      return DST7;
    }

    return DCT2;
  }

  // TODO: SBT handling.

  if (explicit_mts && tu->tr_idx > MTS_SKIP) {
    static const tr_type_t mts_subset_intra[4][2] = {{DST7, DST7}, {DCT8, DST7}, {DST7, DCT8}, {DCT8, DCT8}};
    return mts_subset_intra[tu->tr_idx - MTS_DST7_DST7][0];
  }

  return DCT2;
}

inline tr_type_t get_vertical_transform_type(
  const cu_info_t* __restrict tu,
  int height,
  color_t color,
  bool is_isp,
  bool is_lfnst,
  bool explicit_mts,
  bool implicit_mts)
{
  if (color != COLOR_Y) {
    return DCT2;
  }

  if (is_isp && is_lfnst) {
    return DCT2;
  }

  if (implicit_mts || (is_isp && explicit_mts)) {
    if ((height >= 4) && (height <= 16)) {
      return DST7;
    }
    return DCT2;
  }

  // TODO: SBT handling.

  if (explicit_mts && tu->tr_idx > MTS_SKIP) {
    static const tr_type_t mts_subset_intra[4][2] = {{DST7, DST7}, {DCT8, DST7}, {DST7, DCT8}, {DCT8, DCT8}};
    return mts_subset_intra[tu->tr_idx - MTS_DST7_DST7][1];
  }

  return DCT2;
}

#endif //STRATEGIES_TRANSFORM_COMMON_NEON_HPP_
