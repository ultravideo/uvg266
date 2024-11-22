#ifndef STRATEGIES_INTRA_H_
#define STRATEGIES_INTRA_H_
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
 * Interface for intra prediction functions.
 */

#include "cu.h"
#include "global.h" // IWYU pragma: keep
#include "intra.h"
#include "uvg266.h"


typedef void (angular_pred_func)(
  const cu_loc_t* const cu_loc,
  const int_fast8_t intra_mode,
  const int_fast8_t channel_type,
  const uvg_pixel *const in_ref_above,
  const uvg_pixel *const in_ref_left,
  uvg_pixel *const dst,
  const uint8_t multi_ref_idx,
  const uint8_t isp_mode,
  const int cu_dim);

typedef void (intra_pred_planar_func)(
  const cu_loc_t* const cu_loc,
  color_t color,
  const uvg_pixel *const ref_top,
  const uvg_pixel *const ref_left,
  uvg_pixel *const dst);

typedef void (intra_pred_filtered_dc_func)(
  const int_fast8_t log2_width,
  const uvg_pixel *const ref_top,
  const uvg_pixel *const ref_left,
  uvg_pixel *const out_block,
  const uint8_t multi_ref_idx);

typedef void (pdpc_planar_dc_func)(
  const int mode,
  const cu_loc_t* const cu_loc,
  const color_t color,
  const uvg_intra_ref *const used_ref,
  uvg_pixel *const dst);

typedef void(mip_pred_func)(
  const uvg_intra_references * const refs,
  const uint16_t                     pred_block_width,
  const uint16_t                     pred_block_height,
  uvg_pixel                         *dst,
  const int                          mip_mode,
  const bool                         mip_transp);

// Declare function pointers.
extern angular_pred_func * uvg_angular_pred;
extern intra_pred_planar_func * uvg_intra_pred_planar;
extern intra_pred_filtered_dc_func * uvg_intra_pred_filtered_dc;
extern pdpc_planar_dc_func * uvg_pdpc_planar_dc;
extern mip_pred_func *uvg_mip_predict;

int uvg_strategy_register_intra(void* opaque, uint8_t bitdepth);


#define STRATEGIES_INTRA_EXPORTS \
  {"angular_pred", (void**) &uvg_angular_pred}, \
  {"intra_pred_planar", (void**) &uvg_intra_pred_planar}, \
  {"intra_pred_filtered_dc", (void**) &uvg_intra_pred_filtered_dc}, \
  {"pdpc_planar_dc", (void**) &uvg_pdpc_planar_dc}, \
  {"mip_predict", (void**) &uvg_mip_predict},



#endif //STRATEGIES_INTRA_H_
