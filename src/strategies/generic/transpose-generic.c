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

#include "transpose-generic.h"
#include "strategyselector.h"

void uvg_transpose_8bit_generic(const void* __restrict src_, void* __restrict dst_, int width, int height)
{
  uint8_t* dst = dst_;
  const uint8_t* src = src_;

  if (width == 1 || height == 1) {
    memcpy(dst, src, (size_t)(width * height) * sizeof(uint8_t));
    return;
  }

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      dst[y + x * height] = src[x + y * width];
    }
  }
}

void uvg_transpose_16bit_generic(const void* __restrict src_, void* __restrict dst_, int width, int height)
{
  uint16_t* dst = dst_;
  const uint16_t* src = src_;

  if (width == 1 || height == 1) {
    memcpy(dst, src, (size_t)(width * height) * sizeof(uint16_t));
    return;
  }

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      dst[y + x * height] = src[x + y * width];
    }
  }
}

void uvg_transpose_px_generic(const uvg_pixel* __restrict src, uvg_pixel* __restrict dst, int width, int height)
{
  if (width == 1 || height == 1) {
    memcpy(dst, src, (size_t)(width * height) * sizeof(uvg_pixel));
    return;
  }

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      dst[y + x * height] = src[x + y * width];
    }
  }
}

int uvg_strategy_register_transpose_generic(void* opaque, uint8_t bitdepth)
{
  bool success = true;

  success &= uvg_strategyselector_register(opaque, "transpose_8bit", "generic", 0, &uvg_transpose_8bit_generic);
  success &= uvg_strategyselector_register(opaque, "transpose_16bit", "generic", 0, &uvg_transpose_16bit_generic);
  success &= uvg_strategyselector_register(opaque, "transpose_px", "generic", 0, &uvg_transpose_px_generic);

  return success;
}
