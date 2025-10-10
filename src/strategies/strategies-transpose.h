#ifndef STRATEGIES_TRANSPOSE_H_
#define STRATEGIES_TRANSPOSE_H_
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
 * Interface for matrix transposition functions.
 */

#include "global.h" // IWYU pragma: keep
#include "uvg266.h"

typedef void (transpose_8bit_func)(
  const void* __restrict src,
  void* __restrict dst,
  int width,
  int height
);

typedef void (transpose_16bit_func)(
  const void* __restrict src,
  void* __restrict dst,
  int width,
  int height
);

typedef void (transpose_px_func)(
  const uvg_pixel* __restrict src,
  uvg_pixel* __restrict dst,
  int width,
  int height
);

extern transpose_8bit_func* uvg_transpose_8bit;
extern transpose_16bit_func* uvg_transpose_16bit;
extern transpose_px_func* uvg_transpose_px;

int uvg_strategy_register_transpose(void* opaque, uint8_t bitdepth);

#define STRATEGIES_TRANSPOSE_EXPORTS \
  {"transpose_8bit", (void**)&uvg_transpose_8bit}, \
  {"transpose_16bit", (void**)&uvg_transpose_16bit}, \
  {"transpose_px", (void**)&uvg_transpose_px}, \

#endif //STRATEGIES_TRANSPOSE_H_
