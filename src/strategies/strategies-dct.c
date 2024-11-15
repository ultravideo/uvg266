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
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

#include "strategies/strategies-dct.h"
#if defined(__AVX512F__)
#include "strategies/avx2/dct-avx2.h"
#endif
#include "strategies/generic/dct-generic.h"
#include "strategyselector.h"


// Define function pointers.
dct_func * uvg_fast_forward_dst_4x4 = 0;

dct_func * uvg_dct_4x4 = 0;
dct_func * uvg_dct_8x8 = 0;
dct_func * uvg_dct_16x16 = 0;
dct_func * uvg_dct_32x32 = 0;
dct_func * uvg_dct_non_square = 0;

dct_func * uvg_fast_inverse_dst_4x4 = 0;

dct_func * uvg_idct_4x4 = 0;
dct_func * uvg_idct_8x8= 0;
dct_func * uvg_idct_16x16 = 0;
dct_func * uvg_idct_32x32 = 0;

void(*uvg_mts_dct)(int8_t bitdepth,
  color_t color,
  const cu_info_t *tu,
  int8_t width,
  int8_t height,
  const int16_t *input,
  int16_t *output,
  const int8_t mts_type);

void(*uvg_mts_idct)(int8_t bitdepth,
  color_t color,
  const cu_info_t *tu,
  int8_t width,
  int8_t height,
  const int16_t *input,
  int16_t *output,
  const int8_t mts_type);


int uvg_strategy_register_dct(void* opaque, uint8_t bitdepth) {
  bool success = true;

  success &= uvg_strategy_register_dct_generic(opaque, bitdepth);
#if defined(__AVX512F__)
  if (uvg_g_hardware_flags.intel_flags.avx2) {
    success &= uvg_strategy_register_dct_avx2(opaque, bitdepth);
  }
#endif
  return success;
}


/**
 * \brief  Get a function that performs the transform for a block.
 *
 * \param width    Width of the region
 * \param color    Color plane
 * \param type     Prediction type
 *
 * \returns Pointer to the function.
 */
dct_func * uvg_get_dct_func(int8_t width, int8_t height, color_t color, cu_type_t type)
{
  if (width != height) {
    // Non-square block. Return generic dct for non-square blokcs.
    assert(false && "This should never be called at this point. Non-square stuff is done inside mts_dct function.");
    //return uvg_dct_non_square;
  }
  switch (width) {
  case 4:
    //if (color == COLOR_Y && type == CU_INTRA) {
    //  return uvg_fast_forward_dst_4x4;
    //} else {
      return uvg_dct_4x4;
    //}
  case 8:
    return uvg_dct_8x8;
  case 16:
    return uvg_dct_16x16;
  case 32:
    return uvg_dct_32x32;
  default:
    return NULL;
  }
}

/**
 * \brief  Get a function that performs the inverse transform for a block.
 *
 * \param width    Width of the region
 * \param color    Color plane
 * \param type     Prediction type
 *
 * \returns Pointer to the function.
 */
dct_func * uvg_get_idct_func(int8_t width, int8_t height, color_t color, cu_type_t type)
{
  if (width != height) {
    // Non-square block. Return generic dct for non-square blokcs.
    assert(false && "This should never be called at this point. Non-square stuff is done inside mts_idct function.");
    //return uvg_idct_non_square;
  }
  switch (width) {
  case 4:
    //if (color == COLOR_Y && type == CU_INTRA) {
    //  return uvg_fast_inverse_dst_4x4;
    //} else {
      return uvg_idct_4x4;
    //}
  case 8:
    return uvg_idct_8x8;
  case 16:
    return uvg_idct_16x16;
  case 32:
    return uvg_idct_32x32;
  default:
    return NULL;
  }
}
