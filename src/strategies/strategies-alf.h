#pragma once
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
 * Interface for alf functions.
 */

#include "encoder.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "uvg266.h"
#include "alf.h"


// Declare function pointers.
typedef void (alf_derive_classification_blk_func)(encoder_state_t * const state,
  const int shift,
  const int n_height,
  const int n_width,
  const int blk_pos_x,
  const int blk_pos_y,
  const int blk_dst_x,
  const int blk_dst_y,
  const int vb_ctu_height,
  int vb_pos);

typedef void (alf_filter_5x5_blk_func)(encoder_state_t* const state,
  const uvg_pixel* src_pixels,
  uvg_pixel* dst_pixels,
  const int src_stride,
  const int dst_stride,
  const short* filter_set,
  const int16_t* fClipSet,
  clp_rng clp_rng,
  const int width,
  const int height,
  int x_pos,
  int y_pos,
  int blk_dst_x,
  int blk_dst_y,
  int vb_pos,
  const int vb_ctu_height);

typedef void (alf_filter_7x7_blk_func)(encoder_state_t* const state,
  const uvg_pixel* src_pixels,
  uvg_pixel* dst_pixels,
  const int src_stride,
  const int dst_stride,
  const short* filter_set,
  const int16_t* fClipSet,
  clp_rng clp_rng,
  const int width,
  const int height,
  int x_pos,
  int y_pos,
  int blk_dst_x,
  int blk_dst_y,
  int vb_pos,
  const int vb_ctu_height);

typedef void (alf_get_blk_stats_func)(encoder_state_t* const state,
  channel_type channel,
  alf_covariance* alf_covariance,
  alf_classifier** g_classifier,
  uvg_pixel* org,
  int32_t org_stride,
  uvg_pixel* rec,
  int32_t rec_stride,
  const int x_pos,
  const int y_pos,
  const int x_dst,
  const int y_dst,
  const int width,
  const int height,
  int vb_ctu_height,
  int vb_pos,
  short alf_clipping_values[MAX_NUM_CHANNEL_TYPE][MAX_ALF_NUM_CLIPPING_VALUES]);

// Declare function pointers.
extern alf_derive_classification_blk_func * uvg_alf_derive_classification_blk;
extern alf_filter_5x5_blk_func* uvg_alf_filter_5x5_blk;
extern alf_filter_7x7_blk_func* uvg_alf_filter_7x7_blk;
extern alf_get_blk_stats_func* uvg_alf_get_blk_stats;

int uvg_strategy_register_alf(void* opaque, uint8_t bitdepth);


#define STRATEGIES_ALF_EXPORTS \
  {"alf_derive_classification_blk", (void**) &uvg_alf_derive_classification_blk}, \
  {"alf_filter_5x5_blk", (void**) &uvg_alf_filter_5x5_blk}, \
  {"alf_filter_7x7_blk", (void**) &uvg_alf_filter_7x7_blk}, \
  {"alf_get_blk_stats", (void**) &uvg_alf_get_blk_stats}, \
 

