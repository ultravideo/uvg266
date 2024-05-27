#ifndef TRANSFORM_H_
#define TRANSFORM_H_
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
 * \ingroup Reconstruction
 * \file
 * Quantization and transform functions.
 */

#include "cu.h"
#include "encoder.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep

extern const uint8_t uvg_g_chroma_scale[58];
extern const int16_t uvg_g_inv_quant_scales[2][6];
extern const int16_t uvg_g_quant_scales[2][6];

#define COEFF_ORDER_LINEAR 0
#define COEFF_ORDER_CU 1

void uvg_transformskip(const encoder_control_t *const encoder, int16_t *block,int16_t *coeff, int8_t width, int8_t height);
void uvg_itransformskip(const encoder_control_t *const encoder, int16_t *block,int16_t *coeff, int8_t width, int8_t height);

void uvg_transform2d(const encoder_control_t * const encoder,
                     int16_t *block,
                     int16_t *coeff,
                     int8_t block_width,
                     int8_t block_height,
                     color_t color,
                     const cu_info_t *tu);

void uvg_itransform2d(const encoder_control_t * const encoder,
                      int16_t *block,
                      int16_t *coeff,
                      int8_t block_width,
                      int8_t block_height,
                      color_t color,
                      const cu_info_t *tu);


int32_t uvg_get_scaled_qp(color_t color, int8_t qp, int8_t qp_offset, int8_t const* const chroma_scale);

void uvg_derive_lfnst_constraints(
  cu_info_t* const pred_cu,
  bool* constraints,
  const coeff_t* coeff,
  const int width,
  const int height,
  const vector2d_t * const ,
  color_t color);

typedef struct {
  double best_u_cost;
  double best_v_cost;
  double best_combined_cost;
  int best_u_index;
  int best_v_index;
  int best_combined_index;
  uint64_t u_distortion;
  uint64_t v_distortion;
  double   u_bits;
  double   v_bits;
} uvg_chorma_ts_out_t;

void uvg_quantize_lcu_residual(
  encoder_state_t *const state,
  const bool luma,
  const bool chroma,
  const bool jccr,
  const cu_loc_t* cu_loc,
  cu_info_t *cur_cu,
  lcu_t* lcu,
  bool early_skip,
  enum uvg_tree_type tree_type);

void uvg_chroma_transform_search(
  encoder_state_t* const state,
  lcu_t* const lcu,
  cabac_data_t* temp_cabac,
  const cu_loc_t* const cu_loc,
  const int offset,
  cu_info_t* pred_cu,
  uvg_pixel u_pred[1024],
  uvg_pixel v_pred[1024],
  int16_t u_resi[1024],
  int16_t v_resi[1024],
  uvg_chorma_ts_out_t* chorma_ts_out,
  enum uvg_tree_type tree_type);

enum uvg_chroma_transforms {
  DCT7_CHROMA = 0,
  CHROMA_TS = 4,
  NO_RESIDUAL = 8,
  JCCR_1 = 1,
  JCCR_2 = 2,
  JCCR_3 = 3,
};

void uvg_fwd_lfnst(
  const cu_info_t* const cur_cu,
  const int width,
  const int height,
  const color_t color,
  const uint16_t lfnst_idx,
  coeff_t *coeffs,
  enum uvg_tree_type tree_type,
  int8_t luma_mode);

void uvg_inv_lfnst(
  const cu_info_t* cur_cu,
  const int width,
  const int height,
  const color_t color,
  const uint16_t lfnst_idx,
  coeff_t* coeffs,
  enum uvg_tree_type tree_type,
  int8_t luma_mode);

#endif
