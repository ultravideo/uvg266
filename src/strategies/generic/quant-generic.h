#ifndef STRATEGIES_QUANT_GENERIC_H_
#define STRATEGIES_QUANT_GENERIC_H_
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
 * Generic C implementations of optimized functions.
 */

#include "cu.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "uvg266.h"
#include "tables.h"

int uvg_strategy_register_quant_generic(void* opaque, uint8_t bitdepth);
void uvg_quant_generic(
  const encoder_state_t * const state,
  coeff_t *coef,
  coeff_t *q_coef,
  int32_t width,
  int32_t height,
  color_t color,
  int8_t scan_idx,
  int8_t block_type,
  int8_t transform_skip,
  uint8_t lfnst_idx);

int uvg_quantize_residual_generic(encoder_state_t *const state,
  const cu_info_t *const cur_cu, const int width, const int height, const color_t color,
  const coeff_scan_order_t scan_order, const int use_trskip,
  const int in_stride, const int out_stride,
  const uvg_pixel *const ref_in, const uvg_pixel *const pred_in,
  uvg_pixel *rec_out, coeff_t *coeff_out,
  bool early_skip, int lmcs_chroma_adj, enum uvg_tree_type tree_type);

int uvg_quant_cbcr_residual_generic(
  encoder_state_t* const state,
  const cu_info_t* const cur_cu,
  const int width,
  const int height,
  const coeff_scan_order_t scan_order,
  const int in_stride, const int out_stride,
  const uvg_pixel* const u_ref_in,
  const uvg_pixel* const v_ref_in,
  const uvg_pixel* const u_pred_in,
  const uvg_pixel* const v_pred_in,
  uvg_pixel* u_rec_out,
  uvg_pixel* v_rec_out,
  coeff_t* coeff_out,
  bool early_skip,
  int lmcs_chroma_adj, 
  enum uvg_tree_type tree_type
);

#endif //STRATEGIES_QUANT_GENERIC_H_
