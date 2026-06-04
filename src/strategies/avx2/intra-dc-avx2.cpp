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


// Following functions calculate a sum of adjacent reference pixel values.
// Amount of pixels specified by function name, i.e. width 4 (w4) sums four adjacent pixels.
// These are used when only one reference side needs to be summed.
static INLINE int16_t dc_sum_w4_avx2(const uvg_pixel* const ref)
{
  __m128i vref = _mm_loadu_si128((const __m128i*)ref);
  __m128i vref16 = _mm_cvtepu8_epi16(vref);
  __m128i vhadd = _mm_hadd_epi16(vref16, vref16);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);

  int16_t sum = _mm_extract_epi16(vhadd, 0);
  return sum;
}

static INLINE int16_t dc_sum_w8_avx2(const uvg_pixel* const ref)
{
  __m128i vref = _mm_loadu_si128((const __m128i*)ref);
  __m128i vref16 = _mm_cvtepu8_epi16(vref);
  __m128i vhadd = _mm_hadd_epi16(vref16, vref16);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);

  
  int16_t sum = _mm_extract_epi16(vhadd, 0);
  return sum;
}

static INLINE int16_t dc_sum_w16_avx2(const uvg_pixel* const ref)
{
  const __m128i vshuf = _mm_setr_epi8(
    0x08, 0xff, 0x09, 0xff, 0x0a, 0xff, 0x0b, 0xff,
    0x0c, 0xff, 0x0d, 0xff, 0x0e, 0xff, 0x0f, 0xff
  );
  __m128i vrefa = _mm_loadu_si128((const __m128i*)&ref[0]);
  __m128i vref16b = _mm_shuffle_epi8(vrefa, vshuf);
  __m128i vref16a = _mm_cvtepu8_epi16(vrefa);

  __m128i vadd = _mm_add_epi16(vref16a, vref16b);
  __m128i vhadd = _mm_hadd_epi16(vadd, vadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);
  
  int16_t sum = _mm_extract_epi16(vhadd, 0);
  return sum;
}

static INLINE int16_t dc_sum_w32_avx2(const uvg_pixel* const ref)
{
  const __m128i vshuf = _mm_setr_epi8(
    0x08, 0xff, 0x09, 0xff, 0x0a, 0xff, 0x0b, 0xff,
    0x0c, 0xff, 0x0d, 0xff, 0x0e, 0xff, 0x0f, 0xff
  );
  __m128i vrefa = _mm_loadu_si128((const __m128i*) & ref[0]);
  __m128i vref16b = _mm_shuffle_epi8(vrefa, vshuf);
  __m128i vrefc = _mm_loadu_si128((const __m128i*) & ref[16]);
  __m128i vref16d = _mm_shuffle_epi8(vrefc, vshuf);

  __m128i vref16a = _mm_cvtepu8_epi16(vrefa);
  __m128i vref16c = _mm_cvtepu8_epi16(vrefc);

  __m128i vaddab = _mm_add_epi16(vref16a, vref16b);
  __m128i vaddcd = _mm_add_epi16(vref16c, vref16d);

  __m128i vadd = _mm_add_epi16(vaddab, vaddcd);
  __m128i vhadd = _mm_hadd_epi16(vadd, vadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);
  
  int16_t sum = _mm_extract_epi16(vhadd, 0);
  return sum;
}

static INLINE int16_t dc_sum_w64_avx2(const uvg_pixel* const ref)
{
  const __m128i vshuf = _mm_setr_epi8(
    0x08, 0xff, 0x09, 0xff, 0x0a, 0xff, 0x0b, 0xff,
    0x0c, 0xff, 0x0d, 0xff, 0x0e, 0xff, 0x0f, 0xff
  );
  __m128i vrefa = _mm_loadu_si128((const __m128i*) &ref[0]);
  __m128i vref16b = _mm_shuffle_epi8(vrefa, vshuf);
  __m128i vrefc = _mm_loadu_si128((const __m128i*) &ref[16]);
  __m128i vref16d = _mm_shuffle_epi8(vrefc, vshuf);
  __m128i vrefe = _mm_loadu_si128((const __m128i*) &ref[32]);
  __m128i vref16f = _mm_shuffle_epi8(vrefe, vshuf);
  __m128i vrefg = _mm_loadu_si128((const __m128i*) &ref[48]);
  __m128i vref16h = _mm_shuffle_epi8(vrefg, vshuf);

  __m128i vref16a = _mm_cvtepu8_epi16(vrefa);
  __m128i vref16c = _mm_cvtepu8_epi16(vrefc);
  __m128i vref16e = _mm_cvtepu8_epi16(vrefe);
  __m128i vref16g = _mm_cvtepu8_epi16(vrefg);

  __m128i vadd_ab = _mm_add_epi16(vref16a, vref16b);
  __m128i vadd_cd = _mm_add_epi16(vref16c, vref16d);
  __m128i vadd_ef = _mm_add_epi16(vref16e, vref16f);
  __m128i vadd_gh = _mm_add_epi16(vref16g, vref16h);

  __m128i vadd_abcd = _mm_add_epi16(vadd_ab, vadd_cd);
  __m128i vadd_efgh = _mm_add_epi16(vadd_ef, vadd_gh);

  __m128i vadd = _mm_add_epi16(vadd_abcd, vadd_efgh);

  __m128i vhadd = _mm_hadd_epi16(vadd, vadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);

  int16_t sum = _mm_extract_epi16(vhadd, 0);
  return sum;
}

// These functions do the same as above, but for square blocks where
// both references are summed
static INLINE int16_t dc_sum_square_w4_avx2(const uvg_pixel* const ref_top, const uvg_pixel* const ref_left)
{
  __m128i vrefa = _mm_loadu_si128((const __m128i*)ref_top);
  __m128i vrefb = _mm_loadu_si128((const __m128i*)ref_left);
  __m128i vref16a = _mm_cvtepu8_epi16(vrefa);
  __m128i vref16b = _mm_cvtepu8_epi16(vrefb);

  __m128i vadd = _mm_add_epi16(vref16a, vref16b);

  __m128i vhadd = _mm_hadd_epi16(vadd, vadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);

  int16_t sum = _mm_extract_epi16(vhadd, 0);
  return sum;
}

static INLINE int16_t dc_sum_square_w8_avx2(const uvg_pixel* const ref_top, const uvg_pixel* const ref_left)
{
  __m128i vrefa = _mm_loadu_si128((const __m128i*)ref_top);
  __m128i vrefb = _mm_loadu_si128((const __m128i*)ref_left);
  __m128i vref16a = _mm_cvtepu8_epi16(vrefa);
  __m128i vref16b = _mm_cvtepu8_epi16(vrefb);

  __m128i vadd = _mm_add_epi16(vref16a, vref16b);

  __m128i vhadd = _mm_hadd_epi16(vadd, vadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);

  int16_t sum = _mm_extract_epi16(vhadd, 0);
  return sum;
}

static INLINE int16_t dc_sum_square_w16_avx2(const uvg_pixel* const ref_top, const uvg_pixel* const ref_left)
{
  const __m128i vshuf = _mm_setr_epi8(
    0x08, 0xff, 0x09, 0xff, 0x0a, 0xff, 0x0b, 0xff,
    0x0c, 0xff, 0x0d, 0xff, 0x0e, 0xff, 0x0f, 0xff
  );
  __m128i vref_topa = _mm_loadu_si128((const __m128i*)ref_top);
  __m128i vref16_topb = _mm_shuffle_epi8(vref_topa, vshuf);
  __m128i vref_lefta = _mm_loadu_si128((const __m128i*)ref_left);
  __m128i vref16_leftb = _mm_shuffle_epi8(vref_lefta, vshuf);
  
  __m128i vref16_topa = _mm_cvtepu8_epi16(vref_topa);
  __m128i vref16_lefta = _mm_cvtepu8_epi16(vref_lefta);

  __m128i vaddtop = _mm_add_epi16(vref16_topa, vref16_topb);
  __m128i vaddleft = _mm_add_epi16(vref16_lefta, vref16_leftb);

  __m128i vadd = _mm_add_epi16(vaddtop, vaddleft);

  __m128i vhadd = _mm_hadd_epi16(vadd, vadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);

  int16_t sum = _mm_extract_epi16(vhadd, 0);
  return sum;
}

static INLINE int16_t dc_sum_square_w32_avx2(const uvg_pixel* const ref_top, const uvg_pixel* const ref_left)
{
  const __m128i vshuf = _mm_setr_epi8(
    0x08, 0xff, 0x09, 0xff, 0x0a, 0xff, 0x0b, 0xff,
    0x0c, 0xff, 0x0d, 0xff, 0x0e, 0xff, 0x0f, 0xff
  );
  __m128i vref_topa = _mm_loadu_si128((const __m128i*) &ref_top[0]);
  __m128i vref16_topb = _mm_shuffle_epi8(vref_topa, vshuf);
  __m128i vref_topc = _mm_loadu_si128((const __m128i*) &ref_top[16]);
  __m128i vref16_topd = _mm_shuffle_epi8(vref_topc, vshuf);
  
  __m128i vref_lefta = _mm_loadu_si128((const __m128i*) &ref_left[0]);
  __m128i vref16_leftb = _mm_shuffle_epi8(vref_lefta, vshuf);
  __m128i vref_leftc = _mm_loadu_si128((const __m128i*) &ref_left[16]);
  __m128i vref16_leftd = _mm_shuffle_epi8(vref_leftc, vshuf);

  __m128i vref16_topa = _mm_cvtepu8_epi16(vref_topa);
  __m128i vref16_topc = _mm_cvtepu8_epi16(vref_topc);
  __m128i vref16_lefta = _mm_cvtepu8_epi16(vref_lefta);
  __m128i vref16_leftc = _mm_cvtepu8_epi16(vref_leftc);

  __m128i vaddtop_ab = _mm_add_epi16(vref16_topa, vref16_topb);
  __m128i vaddtop_cd = _mm_add_epi16(vref16_topc, vref16_topd);
  __m128i vaddleft_ab = _mm_add_epi16(vref16_lefta, vref16_leftb);
  __m128i vaddleft_cd = _mm_add_epi16(vref16_leftc, vref16_leftd);

  __m128i vaddtop_abcd = _mm_add_epi16(vaddtop_ab, vaddtop_cd);
  __m128i vaddleft_abcd = _mm_add_epi16(vaddleft_ab, vaddleft_cd);

  __m128i vadd = _mm_add_epi16(vaddtop_abcd, vaddleft_abcd);

  __m128i vhadd = _mm_hadd_epi16(vadd, vadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);

  int16_t sum = _mm_extract_epi16(vhadd, 0);
  return sum;
}

static INLINE int16_t dc_sum_square_w64_avx2(const uvg_pixel* const ref_top, const uvg_pixel* const ref_left)
{
  const __m128i vshuf = _mm_setr_epi8(
    0x08, 0xff, 0x09, 0xff, 0x0a, 0xff, 0x0b, 0xff,
    0x0c, 0xff, 0x0d, 0xff, 0x0e, 0xff, 0x0f, 0xff
  );

  // Sum top ref
  __m128i vrefa = _mm_loadu_si128((const __m128i*) &ref_top[0]);
  __m128i vrefb = _mm_shuffle_epi8(vrefa, vshuf);
  __m128i vrefc = _mm_loadu_si128((const __m128i*) &ref_top[16]);
  __m128i vrefd = _mm_shuffle_epi8(vrefc, vshuf);
  __m128i vrefe = _mm_loadu_si128((const __m128i*) &ref_top[32]);
  __m128i vreff = _mm_shuffle_epi8(vrefe, vshuf);
  __m128i vrefg = _mm_loadu_si128((const __m128i*) &ref_top[48]);
  __m128i vrefh = _mm_shuffle_epi8(vrefg, vshuf);

  vrefa = _mm_cvtepu8_epi16(vrefa);
  vrefc = _mm_cvtepu8_epi16(vrefc);
  vrefe = _mm_cvtepu8_epi16(vrefe);
  vrefg = _mm_cvtepu8_epi16(vrefg);

  __m128i vadd_ab = _mm_add_epi16(vrefa, vrefb);
  __m128i vadd_cd = _mm_add_epi16(vrefc, vrefd);
  __m128i vadd_ef = _mm_add_epi16(vrefe, vreff);
  __m128i vadd_gh = _mm_add_epi16(vrefg, vrefh);

  __m128i vadd_abcd = _mm_add_epi16(vadd_ab, vadd_cd);
  __m128i vadd_efgh = _mm_add_epi16(vadd_ef, vadd_gh);

  __m128i vaddtop = _mm_add_epi16(vadd_abcd, vadd_efgh);

  // Sum left ref
  vrefa = _mm_loadu_si128((const __m128i*) &ref_left[0]);
  vrefb = _mm_shuffle_epi8(vrefa, vshuf);
  vrefc = _mm_loadu_si128((const __m128i*) &ref_left[16]);
  vrefd = _mm_shuffle_epi8(vrefc, vshuf);
  vrefe = _mm_loadu_si128((const __m128i*) &ref_left[32]);
  vreff = _mm_shuffle_epi8(vrefe, vshuf);
  vrefg = _mm_loadu_si128((const __m128i*) &ref_left[48]);
  vrefh = _mm_shuffle_epi8(vrefg, vshuf);

  vrefa = _mm_cvtepu8_epi16(vrefa);
  vrefc = _mm_cvtepu8_epi16(vrefc);
  vrefe = _mm_cvtepu8_epi16(vrefe);
  vrefg = _mm_cvtepu8_epi16(vrefg);

  vadd_ab = _mm_add_epi16(vrefa, vrefb);
  vadd_cd = _mm_add_epi16(vrefc, vrefd);
  vadd_ef = _mm_add_epi16(vrefe, vreff);
  vadd_gh = _mm_add_epi16(vrefg, vrefh);

  vadd_abcd = _mm_add_epi16(vadd_ab, vadd_cd);
  vadd_efgh = _mm_add_epi16(vadd_ef, vadd_gh);

  __m128i vaddleft = _mm_add_epi16(vadd_abcd, vadd_efgh);

  // Add top and left sums
  __m128i vadd = _mm_add_epi16(vaddtop, vaddleft);
  __m128i vhadd = _mm_hadd_epi16(vadd, vadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);
  vhadd = _mm_hadd_epi16(vhadd, vhadd);

  int16_t sum = _mm_extract_epi16(vhadd, 0);
  return sum;
}

static void uvg_intra_pred_dc_avx2(
  const cu_loc_t* const cu_loc,
  const color_t color,
  const uvg_pixel* const ref_top,
  const uvg_pixel* const ref_left,
  uvg_pixel* const out_block,
  const uint8_t multi_ref_idx)
{
  const int width = color == COLOR_Y ? cu_loc->width : cu_loc->chroma_width;
  const int height = color == COLOR_Y ? cu_loc->height : cu_loc->chroma_height;

  int16_t sum = 0;
  const uvg_pixel * const refp_left = &ref_left[1 + multi_ref_idx];
  const uvg_pixel * const refp_top = &ref_top[1 + multi_ref_idx];
  const uvg_pixel * const refp = width > height ? refp_top : refp_left;
  const int dim = width >= height ? width : height;
  // Only one loop is done for non-square blocks.
  // In case of non-square blocks, only the longer reference is summed.
  if (width != height) {
    switch (dim) {
      case 4: sum = dc_sum_w4_avx2(refp); break;
      case 8: sum = dc_sum_w8_avx2(refp); break;
      case 16: sum = dc_sum_w16_avx2(refp); break;
      case 32: sum = dc_sum_w32_avx2(refp); break;
      case 64: sum = dc_sum_w64_avx2(refp); break;
      default:
        assert(false && "Intra DC AVX2: Illegal width.\n");
    }
  } 
  else { // Square blocks
    switch (dim) {
    case 4: sum = dc_sum_square_w4_avx2(refp_top, refp_left); break;
    case 8: sum = dc_sum_square_w8_avx2(refp_top, refp_left); break;
    case 16: sum = dc_sum_square_w16_avx2(refp_top, refp_left); break;
    case 32: sum = dc_sum_square_w32_avx2(refp_top, refp_left); break;
    case 64: sum = dc_sum_square_w64_avx2(refp_top, refp_left); break;
    default:
      assert(false && "Intra DC AVX2: Illegal block size.\n");
    }
  }

  // JVET_K0122
  const int denom = width == height ? width << 1 : MAX(width, height);
  const int divShift = uvg_math_floor_log2(denom);
  const int divOffset = denom >> 1;

  const uvg_pixel dc_val = (sum + divOffset) >> divShift;
  //const uvg_pixel dc_val = (sum + width) >> (log2_width + 1);
  const uint16_t block_size = width * height;
  __m256i dc_val_vec = _mm256_set1_epi8(dc_val);

  switch (block_size)
  {
    case 16:
    case 32: _mm256_store_si256((__m256i*)(out_block), dc_val_vec); break;
    case 64: _mm256_store_si256((__m256i*)(out_block), dc_val_vec); _mm256_store_si256((__m256i*)(out_block + 32), dc_val_vec); break;
    case 128: _mm256_store_si256((__m256i*)(out_block), dc_val_vec); _mm256_store_si256((__m256i*)(out_block + 32), dc_val_vec); _mm256_store_si256((__m256i*)(out_block + 64), dc_val_vec); _mm256_store_si256((__m256i*)(out_block + 96), dc_val_vec); break;
    default: {
      for (int i = 0; i < block_size; i+=256) {
        _mm256_store_si256((__m256i*)(out_block + i), dc_val_vec);
        _mm256_store_si256((__m256i*)(out_block + 32 + i), dc_val_vec);
        _mm256_store_si256((__m256i*)(out_block + 64 + i), dc_val_vec);
        _mm256_store_si256((__m256i*)(out_block + 96 + i), dc_val_vec);
        _mm256_store_si256((__m256i*)(out_block + 128 +i), dc_val_vec);
        _mm256_store_si256((__m256i*)(out_block + 160 + i), dc_val_vec);
        _mm256_store_si256((__m256i*)(out_block + 192 + i), dc_val_vec);
        _mm256_store_si256((__m256i*)(out_block + 224 + i), dc_val_vec);
      }
    }
  }
}


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
template<int N>
static inline constexpr int log2_const()
{
  static_assert(N == 4 || N == 8 || N == 16 || N == 32 || N == 64,
    "PDPC size must be power of two from 4 to 64");

  if constexpr (N == 4)  return 2;
  if constexpr (N == 8)  return 3;
  if constexpr (N == 16) return 4;
  if constexpr (N == 32) return 5;
  if constexpr (N == 64) return 6;
}

static inline uint32_t load_u32_unaligned(const void* p)
{
  uint32_t v;
  memcpy(&v, p, sizeof(v));
  return v;
}

static inline void store_u32_unaligned(void* p, uint32_t v)
{
  memcpy(p, &v, sizeof(v));
}

template<int W, int H>
static inline void pdpc_make_weights(int16_t* w)
{
  constexpr int log2_width = log2_const<W>();
  constexpr int log2_height = log2_const<H>();
  constexpr int scale = ((log2_width - 2 + log2_height - 2 + 2) >> 2);
  constexpr int max_size = W > H ? W : H;

  for (int i = 0; i < max_size; ++i) {
    const int shift = (2 * i) >> scale;
    w[i] = (int16_t)(shift >= 6 ? 0 : (32 >> shift));
  }
}
template<int W, int H>
static inline void pdpc_16x2_avx2(
  const uvg_pixel* left,
  const uvg_pixel* top,
  const int16_t* w,
  uvg_pixel* dst,
  const int x,
  const int y)
{
  static_assert(W >= 16);
  static_assert(H >= 2);

  uvg_pixel* dst0 = dst + (y + 0) * W + x;
  uvg_pixel* dst1 = dst + (y + 1) * W + x;

  // Load 16 destination pixels from two rows.
  const __m128i dst0_8 = _mm_loadu_si128((const __m128i*)dst0);
  const __m128i dst1_8 = _mm_loadu_si128((const __m128i*)dst1);

  const __m256i vdst0 = _mm256_cvtepu8_epi16(dst0_8);
  const __m256i vdst1 = _mm256_cvtepu8_epi16(dst1_8);

  // Same top samples for both rows.
  const __m128i top_8 = _mm_loadu_si128((const __m128i*)(top + x + 1));
  const __m256i vtop = _mm256_cvtepu8_epi16(top_8);

  // Same horizontal weights for both rows.
  const __m256i vwL = _mm256_loadu_si256((const __m256i*)(w + x));

  // Row-specific left references.
  const __m256i vleft0 = _mm256_set1_epi16((int16_t)left[y + 1]);
  const __m256i vleft1 = _mm256_set1_epi16((int16_t)left[y + 2]);

  // Row-specific vertical weights.
  const __m256i vwT0 = _mm256_set1_epi16(w[y + 0]);
  const __m256i vwT1 = _mm256_set1_epi16(w[y + 1]);

  const __m256i diff_left0 = _mm256_sub_epi16(vleft0, vdst0);
  const __m256i diff_left1 = _mm256_sub_epi16(vleft1, vdst1);

  const __m256i diff_top0 = _mm256_sub_epi16(vtop, vdst0);
  const __m256i diff_top1 = _mm256_sub_epi16(vtop, vdst1);

  __m256i accu0 = _mm256_add_epi16(
    _mm256_mullo_epi16(vwL, diff_left0),
    _mm256_mullo_epi16(vwT0, diff_top0)
  );

  __m256i accu1 = _mm256_add_epi16(
    _mm256_mullo_epi16(vwL, diff_left1),
    _mm256_mullo_epi16(vwT1, diff_top1)
  );

  const __m256i round = _mm256_set1_epi16(32);

  accu0 = _mm256_add_epi16(accu0, round);
  accu1 = _mm256_add_epi16(accu1, round);

  accu0 = _mm256_srai_epi16(accu0, 6);
  accu1 = _mm256_srai_epi16(accu1, 6);

  accu0 = _mm256_add_epi16(vdst0, accu0);
  accu1 = _mm256_add_epi16(vdst1, accu1);

  // Packs inside 128-bit lanes:
  //   before permute:
  //     qword0 = row0[0..7]
  //     qword1 = row1[0..7]
  //     qword2 = row0[8..15]
  //     qword3 = row1[8..15]
  //
  //   after 0xD8 = [0, 2, 1, 3]:
  //     lower 128 = row0[0..15]
  //     upper 128 = row1[0..15]
  __m256i packed = _mm256_packus_epi16(accu0, accu1);
  packed = _mm256_permute4x64_epi64(packed, 0xD8);

  if constexpr (W == 16) {
    _mm256_storeu_si256((__m256i*)dst0, packed);
  }
  else {
    _mm_storeu_si128((__m128i*)dst0, _mm256_castsi256_si128(packed));
    _mm_storeu_si128((__m128i*)dst1, _mm256_extracti128_si256(packed, 1));
  }
}

template<int W, int H>
static inline void pdpc_8x4_avx2(
  const uvg_pixel* left,
  const uvg_pixel* top,
  const int16_t* w,
  uvg_pixel* dst,
  const int y)
{
  static_assert(W == 8);
  static_assert(H >= 4);

  uvg_pixel* dst0 = dst + (y + 0) * W;
  uvg_pixel* dst1 = dst + (y + 1) * W;
  uvg_pixel* dst2 = dst + (y + 2) * W;
  uvg_pixel* dst3 = dst + (y + 3) * W;

  const __m128i dst0_8 = _mm_loadl_epi64((const __m128i*)dst0);
  const __m128i dst1_8 = _mm_loadl_epi64((const __m128i*)dst1);
  const __m128i dst2_8 = _mm_loadl_epi64((const __m128i*)dst2);
  const __m128i dst3_8 = _mm_loadl_epi64((const __m128i*)dst3);

  // Vector A:
  //   low  lane: row0[0..7]
  //   high lane: row1[0..7]
  //
  // Vector B:
  //   low  lane: row2[0..7]
  //   high lane: row3[0..7]
  const __m256i vdst01 = _mm256_inserti128_si256(
    _mm256_cvtepu8_epi16(dst0_8),
    _mm256_castsi256_si128(_mm256_cvtepu8_epi16(dst1_8)),
    1
  );

  const __m256i vdst23 = _mm256_inserti128_si256(
    _mm256_cvtepu8_epi16(dst2_8),
    _mm256_castsi256_si128(_mm256_cvtepu8_epi16(dst3_8)),
    1
  );

  const __m128i top_8 = _mm_loadl_epi64((const __m128i*)(top + 1));
  __m256i vtop = _mm256_cvtepu8_epi16(top_8);
  vtop = _mm256_permute2x128_si256(vtop, vtop, 0x00);

  const __m128i wL_128 = _mm_loadu_si128((const __m128i*)w);
  const __m256i vwL = _mm256_set_m128i(wL_128, wL_128);

  const __m256i vleft01 = _mm256_set_m128i(
    _mm_set1_epi16((int16_t)left[y + 2]),
    _mm_set1_epi16((int16_t)left[y + 1])
  );

  const __m256i vleft23 = _mm256_set_m128i(
    _mm_set1_epi16((int16_t)left[y + 4]),
    _mm_set1_epi16((int16_t)left[y + 3])
  );

  const __m256i vwT01 = _mm256_set_m128i(
    _mm_set1_epi16(w[y + 1]),
    _mm_set1_epi16(w[y + 0])
  );

  const __m256i vwT23 = _mm256_set_m128i(
    _mm_set1_epi16(w[y + 3]),
    _mm_set1_epi16(w[y + 2])
  );

  const __m256i round = _mm256_set1_epi16(32);

  __m256i accu01 = _mm256_add_epi16(
    _mm256_mullo_epi16(vwL, _mm256_sub_epi16(vleft01, vdst01)),
    _mm256_mullo_epi16(vwT01, _mm256_sub_epi16(vtop, vdst01))
  );

  __m256i accu23 = _mm256_add_epi16(
    _mm256_mullo_epi16(vwL, _mm256_sub_epi16(vleft23, vdst23)),
    _mm256_mullo_epi16(vwT23, _mm256_sub_epi16(vtop, vdst23))
  );

  accu01 = _mm256_add_epi16(_mm256_srai_epi16(_mm256_add_epi16(accu01, round), 6), vdst01);
  accu23 = _mm256_add_epi16(_mm256_srai_epi16(_mm256_add_epi16(accu23, round), 6), vdst23);

  // packed contains:
  //   qword0 = row0[0..7]
  //   qword1 = row2[0..7]
  //   qword2 = row1[0..7]
  //   qword3 = row3[0..7]
  //
  // because packing is lane-local:
  //   low lane  packs accu01.low  + accu23.low  => row0,row2
  //   high lane packs accu01.high + accu23.high => row1,row3
  __m256i packed = _mm256_packus_epi16(accu01, accu23);

  // Reorder qwords from [row0,row2,row1,row3] to [row0,row1,row2,row3].
  packed = _mm256_permute4x64_epi64(packed, 0xD8);

  // Since W == 8, four rows are contiguous in memory.
  _mm256_storeu_si256((__m256i*)dst0, packed);
}

template<int W, int H>
static inline void pdpc_4x4_avx2(
  const uvg_pixel* left,
  const uvg_pixel* top,
  const int16_t* w,
  uvg_pixel* dst,
  const int y,
  const __m256i shuf_mask_byte,
  const __m256i shuf_mask_word)
{
  static_assert(W == 4);

  const uint32_t r0 = load_u32_unaligned(dst + (y + 0) * W);
  const uint32_t r1 = load_u32_unaligned(dst + (y + 1) * W);
  const uint32_t r2 = load_u32_unaligned(dst + (y + 2) * W);
  const uint32_t r3 = load_u32_unaligned(dst + (y + 3) * W);

  const __m128i rows = _mm_setr_epi32((int)r0, (int)r1, (int)r2, (int)r3);
  const __m256i vdst = _mm256_cvtepu8_epi16(rows);

  uint32_t dw_left;
  uint32_t dw_top;
  memcpy(&dw_left, left + y + 1, sizeof(dw_left));
  memcpy(&dw_top, top + 1, sizeof(dw_top));

  __m256i vleft = _mm256_set1_epi32((int)dw_left);
  vleft = _mm256_shuffle_epi8(vleft, shuf_mask_byte);

  __m256i vtop = _mm256_set1_epi32((int)dw_top);
  vtop = _mm256_cvtepu8_epi16(_mm256_castsi256_si128(vtop));

  uint64_t quad_wL;
  uint64_t quad_wT;
  memcpy(&quad_wL, w, sizeof(quad_wL));
  memcpy(&quad_wT, w + y, sizeof(quad_wT));

  __m256i vwL = _mm256_set1_epi64x((int64_t)quad_wL);

  __m256i vwT = _mm256_set1_epi64x((int64_t)quad_wT);
  vwT = _mm256_shuffle_epi8(vwT, shuf_mask_word);

  const __m256i diff_left = _mm256_sub_epi16(vleft, vdst);
  const __m256i diff_top = _mm256_sub_epi16(vtop, vdst);

  __m256i accu = _mm256_add_epi16(
    _mm256_mullo_epi16(vwL, diff_left),
    _mm256_mullo_epi16(vwT, diff_top)
  );

  accu = _mm256_add_epi16(accu, _mm256_set1_epi16(32));
  accu = _mm256_srai_epi16(accu, 6);
  accu = _mm256_add_epi16(vdst, accu);

  const __m128i lo = _mm256_castsi256_si128(accu);
  const __m128i hi = _mm256_extracti128_si256(accu, 1);
  const __m128i packed = _mm_packus_epi16(lo, hi);

  store_u32_unaligned(dst + (y + 0) * W, (uint32_t)_mm_extract_epi32(packed, 0));
  store_u32_unaligned(dst + (y + 1) * W, (uint32_t)_mm_extract_epi32(packed, 1));
  store_u32_unaligned(dst + (y + 2) * W, (uint32_t)_mm_extract_epi32(packed, 2));
  store_u32_unaligned(dst + (y + 3) * W, (uint32_t)_mm_extract_epi32(packed, 3));
}
template<int W, int H>
static inline void uvg_pdpc_planar_dc_avx2_size(
  const uvg_intra_ref* const used_ref,
  uvg_pixel* const dst)
{
  alignas(32) int16_t w[64];
  pdpc_make_weights<W, H>(w);

  if constexpr (W >= 16) {
    for (int y = 0; y < H; y += 2) {
      for (int x = 0; x < W; x += 16) {
        pdpc_16x2_avx2<W, H>(
          used_ref->left,
          used_ref->top,
          w,
          dst,
          x,
          y
        );
        if (y > 12) break;
      }
    }
  }
  else if constexpr (W == 8) {
    for (int y = 0; y < H; y += 4) {
      pdpc_8x4_avx2<W, H>(
        used_ref->left,
        used_ref->top,
        w,
        dst,
        y
      );
    }
  }
  else {
    static_assert(W == 4);

    const __m256i shuf_mask_byte = _mm256_setr_epi8(
      0, -1, 0, -1, 0, -1, 0, -1,
      1, -1, 1, -1, 1, -1, 1, -1,
      2, -1, 2, -1, 2, -1, 2, -1,
      3, -1, 3, -1, 3, -1, 3, -1
    );

    const __m256i shuf_mask_word = _mm256_setr_epi8(
      0, 1, 0, 1, 0, 1, 0, 1,
      2, 3, 2, 3, 2, 3, 2, 3,
      4, 5, 4, 5, 4, 5, 4, 5,
      6, 7, 6, 7, 6, 7, 6, 7
    );

    for (int y = 0; y < H; y += 4) {
      pdpc_4x4_avx2<W, H>(
        used_ref->left,
        used_ref->top,
        w,
        dst,
        y,
        shuf_mask_byte,
        shuf_mask_word
      );
    }
  }
}

#define DISPATCH_PDPC_HEIGHT(WIDTH)                         \
  switch (height) {                                         \
    case 4:  uvg_pdpc_planar_dc_avx2_size<WIDTH, 4 >(used_ref, dst); break; \
    case 8:  uvg_pdpc_planar_dc_avx2_size<WIDTH, 8 >(used_ref, dst); break; \
    case 16: uvg_pdpc_planar_dc_avx2_size<WIDTH, 16>(used_ref, dst); break; \
    case 32: uvg_pdpc_planar_dc_avx2_size<WIDTH, 32>(used_ref, dst); break; \
    case 64: uvg_pdpc_planar_dc_avx2_size<WIDTH, 64>(used_ref, dst); break; \
    default: assert(0); break;                              \
  }

static void uvg_pdpc_planar_dc_avx2(
  const int mode,
  const cu_loc_t* const cu_loc,
  const color_t color,
  const uvg_intra_ref* const used_ref,
  uvg_pixel* const dst)
{
  assert(mode == 0 || mode == 1);

  const int width = color == COLOR_Y ? cu_loc->width : cu_loc->chroma_width;
  const int height = color == COLOR_Y ? cu_loc->height : cu_loc->chroma_height;

  assert(width == 4 || width == 8 || width == 16 || width == 32 || width == 64);
  assert(height == 4 || height == 8 || height == 16 || height == 32 || height == 64);

  switch (width) {
  case 4:  DISPATCH_PDPC_HEIGHT(4);  break;
  case 8:  DISPATCH_PDPC_HEIGHT(8);  break;
  case 16: DISPATCH_PDPC_HEIGHT(16); break;
  case 32: DISPATCH_PDPC_HEIGHT(32); break;
  case 64: DISPATCH_PDPC_HEIGHT(64); break;
  default: assert(0); break;
  }
}


#endif // UVG_BIT_DEPTH == 8

int uvg_strategy_register_intra_dc_avx2(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if UVG_BIT_DEPTH == 8
  if (bitdepth == 8) {
    success &= uvg_strategyselector_register(opaque, "pdpc_planar_dc", "avx2", 40, reinterpret_cast<void*>(&uvg_pdpc_planar_dc_avx2));
    success &= uvg_strategyselector_register(opaque, "intra_dc", "avx2", 40, reinterpret_cast<void*>(&uvg_intra_pred_dc_avx2));
  }
#endif //UVG_BIT_DEPTH == 8
  return success;
}
