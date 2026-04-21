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

#include "transform-neon.h"
#include "strategyselector.h"

int uvg_strategy_register_dct_neon(void* opaque, uint8_t bitdepth)
{
  bool success = true;
  success &= uvg_strategyselector_register(opaque, "mts_dct", "neon", 10, &uvg_mts_dct_neon);
  success &= uvg_strategyselector_register(opaque, "mts_idct", "neon", 10, &uvg_mts_idct_neon);
  success &= uvg_strategyselector_register(opaque, "dct_4x4", "neon", 10, &uvg_dct_4x4_neon);
  success &= uvg_strategyselector_register(opaque, "dct_8x8", "neon", 10, &uvg_dct_8x8_neon);
  success &= uvg_strategyselector_register(opaque, "dct_16x16", "neon", 10, &uvg_dct_16x16_neon);
  success &= uvg_strategyselector_register(opaque, "dct_32x32", "neon", 10, &uvg_dct_32x32_neon);
  success &= uvg_strategyselector_register(opaque, "idct_4x4", "neon", 10, &uvg_idct_4x4_neon);
  success &= uvg_strategyselector_register(opaque, "idct_8x8", "neon", 10, &uvg_idct_8x8_neon);
  success &= uvg_strategyselector_register(opaque, "idct_16x16", "neon", 10, &uvg_idct_16x16_neon);
  success &= uvg_strategyselector_register(opaque, "idct_32x32", "neon", 10, &uvg_idct_32x32_neon);
  return success;
}
