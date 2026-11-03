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
 * AVX2 optimizations for DC intra prediction and the planar/DC PDPC filter.
 *
 * This file is C++ to allow compile-time specialisation of the kernels.
 */

#include "strategies/avx2/intra-dc-avx2.h"

extern "C" {
#include "cu.h"
#include "intra.h"
#include "strategyselector.h"
#include "tables.h"
#include "uvg266.h"
#include "uvg_math.h"
}

#if UVG_BIT_DEPTH == 8

#include <immintrin.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>


/**
* \brief Position Dependent Prediction Combination for Planar and DC modes.
* \param mode          Intra mode, 0 for planar, 1 for DC.
* \param cu_loc        Pointer to the CU location information.
* \param color         Color component.
* \param used_ref      Pointer to the used reference pixels.
* \param dst           Buffer of size MAX_PRED_WIDTH * MAX_PRED_WIDTH.
*/
// TODO: allegedly does not work with blocks with height 1 and 2. Test this.
// TODO: or just rework the whole thing. We might be able to optimize this further.
static void uvg_pdpc_planar_dc_avx2(
  const int mode,
  const cu_loc_t* const cu_loc,
  const color_t color,
  const uvg_intra_ref *const used_ref,
  uvg_pixel *const dst)
{
  assert(mode == 0 || mode == 1);  // planar or DC
  const int width = color == COLOR_Y ? cu_loc->width : cu_loc->chroma_width;
  const int height = color == COLOR_Y ? cu_loc->height : cu_loc->chroma_height;
  const int log2_width =  uvg_g_convert_to_log2[width];
  const int log2_height = uvg_g_convert_to_log2[height];

  __m256i shuf_mask_byte = _mm256_setr_epi8(
    0, -1, 0, -1, 0, -1, 0, -1,
    1, -1, 1, -1, 1, -1, 1, -1,
    2, -1, 2, -1, 2, -1, 2, -1,
    3, -1, 3, -1, 3, -1, 3, -1
  );

  __m256i shuf_mask_word = _mm256_setr_epi8(
    0, 1, 0, 1, 0, 1, 0, 1,
    2, 3, 2, 3, 2, 3, 2, 3,
    4, 5, 4, 5, 4, 5, 4, 5,
    6, 7, 6, 7, 6, 7, 6, 7
  );

  const int scale = ((log2_width - 2 + log2_height - 2 + 2) >> 2);

  // Same weights regardless of axis, compute once
  int16_t w[LCU_WIDTH];
  for (int i = 0; i < MAX(width, height); i += 4) {
    __m128i base = _mm_set1_epi32(i);
    __m128i offs = _mm_setr_epi32(0, 1, 2, 3);
    __m128i idxs = _mm_add_epi32(base, offs);
    __m128i unclipped = _mm_slli_epi32(idxs, 1);
    unclipped = _mm_srli_epi32(unclipped, scale);
    __m128i clipped = _mm_min_epi32( _mm_set1_epi32(31), unclipped);
    __m128i weights = _mm_srlv_epi32(_mm_set1_epi32(32), clipped);
    weights = _mm_packus_epi32(weights, weights);
    _mm_storel_epi64((__m128i*)&w[i], weights);
  }

  // Process in 4x4 blocks
  for (int y = 0; y < height; y += 4) {
    for (int x = 0; x < width; x += 4) {

      uint32_t dw_left;
      uint32_t dw_top;
      memcpy(&dw_left, &used_ref->left[y + 1], sizeof(dw_left));
      memcpy(&dw_top , &used_ref->top [x + 1], sizeof(dw_top));
      __m256i vleft = _mm256_set1_epi32(dw_left);
      __m256i vtop  = _mm256_set1_epi32(dw_top);
      vleft  = _mm256_shuffle_epi8(vleft, shuf_mask_byte);
      vtop = _mm256_cvtepu8_epi16(_mm256_castsi256_si128(vtop));

      __m128i vseq = _mm_setr_epi32(0, 1, 2, 3);
      __m128i vidx = _mm_slli_epi32(vseq, log2_width);
      __m128i vdst = _mm_i32gather_epi32((const int32_t*)(dst + y * width + x), vidx, 1);
      __m256i vdst16 = _mm256_cvtepu8_epi16(vdst);
      uint64_t quad_wL;
      uint64_t quad_wT;
      memcpy(&quad_wL, &w[x], sizeof(quad_wL));
      memcpy(&quad_wT, &w[y], sizeof(quad_wT));
      __m256i vwL = _mm256_set1_epi64x(quad_wL); 
      __m256i vwT = _mm256_set1_epi64x(quad_wT);
      vwT = _mm256_shuffle_epi8(vwT, shuf_mask_word);
      __m256i diff_left = _mm256_sub_epi16(vleft, vdst16);
      __m256i diff_top  = _mm256_sub_epi16(vtop , vdst16);
      __m256i prod_left = _mm256_mullo_epi16(vwL, diff_left);
      __m256i prod_top  = _mm256_mullo_epi16(vwT, diff_top);
      __m256i accu = _mm256_add_epi16(prod_left, prod_top);
      accu = _mm256_add_epi16(accu, _mm256_set1_epi16(32));
      accu = _mm256_srai_epi16(accu, 6);
      accu = _mm256_add_epi16(vdst16, accu);

      __m128i lo = _mm256_castsi256_si128(accu);
      __m128i hi = _mm256_extracti128_si256(accu, 1);
      vdst = _mm_packus_epi16(lo, hi);

      *(uint32_t*)(dst + (y + 0) * width + x) = _mm_extract_epi32(vdst, 0);
      *(uint32_t*)(dst + (y + 1) * width + x) = _mm_extract_epi32(vdst, 1);
      *(uint32_t*)(dst + (y + 2) * width + x) = _mm_extract_epi32(vdst, 2);
      *(uint32_t*)(dst + (y + 3) * width + x) = _mm_extract_epi32(vdst, 3);
    }
  }
}


#endif // UVG_BIT_DEPTH == 8

int uvg_strategy_register_intra_dc_avx2(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if UVG_BIT_DEPTH == 8
  if (bitdepth == 8) {
    success &= uvg_strategyselector_register(opaque, "pdpc_planar_dc", "avx2", 40, reinterpret_cast<void*>(&uvg_pdpc_planar_dc_avx2));
  }
#endif //UVG_BIT_DEPTH == 8
  return success;
}
