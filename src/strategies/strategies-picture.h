#ifndef STRATEGIES_PICTURE_H_
#define STRATEGIES_PICTURE_H_
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
 * Interface for distortion metric functions.
 */

#include "global.h" // IWYU pragma: keep
#include "inter.h"
#include "uvg266.h"
#include "encoderstate.h"
#include "strategies/optimized_sad_func_ptr_t.h"


typedef uvg_pixel (*pred_buffer)[32 * 32];

// Function macro for defining hadamard calculating functions
// for fixed size blocks. They calculate hadamard for integer
// multiples of 8x8 with the 8x8 hadamard function.
#define SATD_NxN(suffix, n) \
/* Declare the function in advance, hopefully reducing the probability that the
 * macro expands to something unexpected and silently breaks things. */ \
static cost_pixel_nxn_func satd_ ## n ## x ## n ## _ ## suffix;\
static unsigned satd_ ## n ## x ## n ## _ ## suffix ( \
    const uvg_pixel * const block1, \
    const uvg_pixel * const block2) \
{ \
  unsigned sum = 0; \
  for (unsigned y = 0; y < (n); y += 8) { \
    unsigned row = y * (n); \
    for (unsigned x = 0; x < (n); x += 8) { \
      sum += satd_8x8_subblock_ ## suffix(&block1[row + x], (n), &block2[row + x], (n)); \
    } \
  } \
  return sum >> (UVG_BIT_DEPTH - 8); \
}


// Function macro for defining hadamard calculating functions for dynamic size
// blocks. They calculate hadamard for integer multiples of 8x8 with the 8x8
// hadamard function.
#define SATD_ANY_SIZE(suffix) \
  static cost_pixel_any_size_func satd_any_size_ ## suffix; \
  static unsigned satd_any_size_ ## suffix ( \
      int width, int height, \
      const uvg_pixel *block1, int stride1, \
      const uvg_pixel *block2, int stride2) \
  { \
    unsigned sum = 0; \
    if (width % 8 != 0) { \
      /* Process the first column using 4x4 blocks. */ \
      for (int y = 0; y < height; y += 4) { \
        sum += uvg_satd_4x4_subblock_ ## suffix(&block1[y * stride1], stride1, \
                                                &block2[y * stride2], stride2); \
      } \
      block1 += 4; \
      block2 += 4; \
      width -= 4; \
    } \
    if (height % 8 != 0) { \
      /* Process the first row using 4x4 blocks. */ \
      for (int x = 0; x < width; x += 4) { \
        sum += uvg_satd_4x4_subblock_ ## suffix(&block1[x], stride1, \
                                                &block2[x], stride2); \
      } \
      block1 += 4 * stride1; \
      block2 += 4 * stride2; \
      height -= 4; \
    } \
    /* The rest can now be processed with 8x8 blocks. */ \
    for (int y = 0; y < height; y += 8) { \
      const uvg_pixel *row1 = &block1[y * stride1]; \
      const uvg_pixel *row2 = &block2[y * stride2]; \
      for (int x = 0; x < width; x += 8) { \
        sum += satd_8x8_subblock_ ## suffix(&row1[x], stride1, \
                                            &row2[x], stride2); \
      } \
    } \
    return sum >> (UVG_BIT_DEPTH - 8); \
  }

typedef unsigned(reg_sad_func)(const uvg_pixel *const data1, const uvg_pixel *const data2,
  const int width, const int height,
  const unsigned stride1, const unsigned stride2);
typedef unsigned (cost_pixel_nxn_func)(const uvg_pixel *block1, const uvg_pixel *block2);
typedef unsigned (cost_pixel_any_size_func)(
    int width, int height,
    const uvg_pixel *block1, int stride1,
    const uvg_pixel *block2, int stride2
);
typedef void (cost_pixel_nxn_multi_func)(const pred_buffer preds, const uvg_pixel *orig, unsigned num_modes, unsigned *costs_out);
typedef void (cost_pixel_any_size_multi_func)(int width, int height, const uvg_pixel **preds, const int stride, const uvg_pixel *orig, const int orig_stride, unsigned num_modes, unsigned *costs_out, int8_t *valid);

typedef unsigned (pixels_calc_ssd_func)(const uvg_pixel *const ref, const uvg_pixel *const rec, const int ref_stride, const int rec_stride, const int width, const int height);
typedef optimized_sad_func_ptr_t (get_optimized_sad_func)(int32_t);
typedef uint32_t (ver_sad_func)(const uvg_pixel *pic_data, const uvg_pixel *ref_data,
                                int32_t block_width, int32_t block_height,
                                uint32_t pic_stride);
typedef uint32_t (hor_sad_func)(const uvg_pixel *pic_data, const uvg_pixel *ref_data,
                                int32_t width, int32_t height, uint32_t pic_stride,
                                uint32_t ref_stride, uint32_t left, uint32_t right);

typedef void (inter_recon_bipred_func)(lcu_t * const lcu,
  const yuv_t *const px_L0,
  const yuv_t *const px_L1,
  const yuv_im_t *const im_L0,
  const yuv_im_t *const im_L1,
  const unsigned pu_x,
  const unsigned pu_y,
  const unsigned pu_w,
  const unsigned pu_h,
  const unsigned im_flags_L0,
  const unsigned im_flags_L1,
  const bool predict_luma,
  const bool predict_chroma);

typedef double (pixel_var_func)(const uvg_pixel *buf, const uint32_t len);

typedef void (generate_residual_func)(const uvg_pixel* ref_in, const uvg_pixel* pred_in, int16_t* residual, int width, int height, int ref_stride, int pred_stride);


extern const uint32_t uvg_crc_table[256];

typedef uint32_t(crc32c_4x4_func)(const uvg_pixel *buf, uint32_t pic_stride);
typedef uint32_t(crc32c_8x8_func)(const uvg_pixel *buf, uint32_t pic_stride);

// Declare function pointers.
extern crc32c_4x4_func * uvg_crc32c_4x4;
extern crc32c_8x8_func * uvg_crc32c_8x8;

extern reg_sad_func * uvg_reg_sad;

extern cost_pixel_nxn_func * uvg_sad_4x4;
extern cost_pixel_nxn_func * uvg_sad_8x8;
extern cost_pixel_nxn_func * uvg_sad_16x16;
extern cost_pixel_nxn_func * uvg_sad_32x32;
extern cost_pixel_nxn_func * uvg_sad_64x64;

extern cost_pixel_nxn_func * uvg_satd_4x4;
extern cost_pixel_nxn_func * uvg_satd_8x8;
extern cost_pixel_nxn_func * uvg_satd_16x16;
extern cost_pixel_nxn_func * uvg_satd_32x32;
extern cost_pixel_nxn_func * uvg_satd_64x64;
extern cost_pixel_any_size_func *uvg_satd_any_size;
extern cost_pixel_any_size_func *uvg_satd_any_size_vtm;

extern cost_pixel_nxn_multi_func * uvg_sad_4x4_dual;
extern cost_pixel_nxn_multi_func * uvg_sad_8x8_dual;
extern cost_pixel_nxn_multi_func * uvg_sad_16x16_dual;
extern cost_pixel_nxn_multi_func * uvg_sad_32x32_dual;
extern cost_pixel_nxn_multi_func * uvg_sad_64x64_dual;

extern cost_pixel_nxn_multi_func * uvg_satd_4x4_dual;
extern cost_pixel_nxn_multi_func * uvg_satd_8x8_dual;
extern cost_pixel_nxn_multi_func * uvg_satd_16x16_dual;
extern cost_pixel_nxn_multi_func * uvg_satd_32x32_dual;
extern cost_pixel_nxn_multi_func * uvg_satd_64x64_dual;

extern cost_pixel_any_size_multi_func *uvg_satd_any_size_quad;

extern pixels_calc_ssd_func *uvg_pixels_calc_ssd;

extern inter_recon_bipred_func * uvg_bipred_average;

extern get_optimized_sad_func *uvg_get_optimized_sad;
extern ver_sad_func *uvg_ver_sad;
extern hor_sad_func *uvg_hor_sad;

extern pixel_var_func *uvg_pixel_var;

extern generate_residual_func* uvg_generate_residual;

int uvg_strategy_register_picture(void* opaque, uint8_t bitdepth);
cost_pixel_nxn_multi_func * uvg_pixels_get_satd_dual_func(unsigned width, unsigned height);
cost_pixel_nxn_multi_func * uvg_pixels_get_sad_dual_func(unsigned width, unsigned height);

#if UVG_BIT_DEPTH == 8 
#define STRATEGIES_PICTURE_EXPORTS \
  {"crc32c_4x4", (void**) &uvg_crc32c_4x4}, \
  {"crc32c_8x8", (void **)&uvg_crc32c_8x8}, \
  {"reg_sad", (void**) &uvg_reg_sad}, \
  {"sad_4x4", (void**) &uvg_sad_4x4}, \
  {"sad_8x8", (void**) &uvg_sad_8x8}, \
  {"sad_16x16", (void**) &uvg_sad_16x16}, \
  {"sad_32x32", (void**) &uvg_sad_32x32}, \
  {"sad_64x64", (void**) &uvg_sad_64x64}, \
  {"satd_4x4", (void**) &uvg_satd_4x4}, \
  {"satd_8x8", (void**) &uvg_satd_8x8}, \
  {"satd_16x16", (void**) &uvg_satd_16x16}, \
  {"satd_32x32", (void**) &uvg_satd_32x32}, \
  {"satd_64x64", (void**) &uvg_satd_64x64}, \
  {"satd_any_size", (void**) &uvg_satd_any_size}, \
  {"satd_any_size_vtm", (void**) &uvg_satd_any_size_vtm}, \
  {"sad_4x4_dual", (void**) &uvg_sad_4x4_dual}, \
  {"sad_8x8_dual", (void**) &uvg_sad_8x8_dual}, \
  {"sad_16x16_dual", (void**) &uvg_sad_16x16_dual}, \
  {"sad_32x32_dual", (void**) &uvg_sad_32x32_dual}, \
  {"sad_64x64_dual", (void**) &uvg_sad_64x64_dual}, \
  {"satd_4x4_dual", (void**) &uvg_satd_4x4_dual}, \
  {"satd_8x8_dual", (void**) &uvg_satd_8x8_dual}, \
  {"satd_16x16_dual", (void**) &uvg_satd_16x16_dual}, \
  {"satd_32x32_dual", (void**) &uvg_satd_32x32_dual}, \
  {"satd_64x64_dual", (void**) &uvg_satd_64x64_dual}, \
  {"satd_any_size_quad", (void**) &uvg_satd_any_size_quad}, \
  {"pixels_calc_ssd", (void**) &uvg_pixels_calc_ssd}, \
  {"bipred_average", (void**) &uvg_bipred_average}, \
  {"get_optimized_sad", (void**) &uvg_get_optimized_sad}, \
  {"ver_sad", (void**) &uvg_ver_sad}, \
  {"hor_sad", (void**) &uvg_hor_sad}, \
  {"pixel_var", (void**) &uvg_pixel_var}, \
  {"generate_residual", (void**) &uvg_generate_residual}, \

#else
#define STRATEGIES_PICTURE_EXPORTS \
  {"reg_sad", (void**) &uvg_reg_sad}, \
  {"sad_4x4", (void**) &uvg_sad_4x4}, \
  {"sad_8x8", (void**) &uvg_sad_8x8}, \
  {"sad_16x16", (void**) &uvg_sad_16x16}, \
  {"sad_32x32", (void**) &uvg_sad_32x32}, \
  {"sad_64x64", (void**) &uvg_sad_64x64}, \
  {"satd_4x4", (void**) &uvg_satd_4x4}, \
  {"satd_8x8", (void**) &uvg_satd_8x8}, \
  {"satd_16x16", (void**) &uvg_satd_16x16}, \
  {"satd_32x32", (void**) &uvg_satd_32x32}, \
  {"satd_64x64", (void**) &uvg_satd_64x64}, \
  {"satd_any_size", (void**) &uvg_satd_any_size}, \
  {"satd_any_size_vtm", (void**) &uvg_satd_any_size_vtm}, \
  {"sad_4x4_dual", (void**) &uvg_sad_4x4_dual}, \
  {"sad_8x8_dual", (void**) &uvg_sad_8x8_dual}, \
  {"sad_16x16_dual", (void**) &uvg_sad_16x16_dual}, \
  {"sad_32x32_dual", (void**) &uvg_sad_32x32_dual}, \
  {"sad_64x64_dual", (void**) &uvg_sad_64x64_dual}, \
  {"satd_4x4_dual", (void**) &uvg_satd_4x4_dual}, \
  {"satd_8x8_dual", (void**) &uvg_satd_8x8_dual}, \
  {"satd_16x16_dual", (void**) &uvg_satd_16x16_dual}, \
  {"satd_32x32_dual", (void**) &uvg_satd_32x32_dual}, \
  {"satd_64x64_dual", (void**) &uvg_satd_64x64_dual}, \
  {"satd_any_size_quad", (void**) &uvg_satd_any_size_quad}, \
  {"pixels_calc_ssd", (void**) &uvg_pixels_calc_ssd}, \
  {"bipred_average", (void**) &uvg_bipred_average}, \
  {"get_optimized_sad", (void**) &uvg_get_optimized_sad}, \
  {"ver_sad", (void**) &uvg_ver_sad}, \
  {"hor_sad", (void**) &uvg_hor_sad}, \
  {"pixel_var", (void**) &uvg_pixel_var}, \
  {"generate_residual", (void**) &uvg_generate_residual}, \

#endif 



#endif //STRATEGIES_PICTURE_H_
