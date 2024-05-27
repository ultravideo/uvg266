#ifndef INTER_H_
#define INTER_H_
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
 * Inter prediction.
 */

#include "cu.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "image.h"
#include "uvg266.h"

extern const int8_t uvg_g_imv_to_prec[4];

typedef struct {
  mv_t mv[2][2];
  uint16_t mer[2];
  uint8_t dir;
  uint8_t ref[2]; // index to L0/L1
} inter_merge_cand_t;

void uvg_change_precision(int src, int dst, mv_t* hor, mv_t* ver);
void uvg_change_precision_vector2d(int src, int dst, vector2d_t* mv);
void uvg_round_precision(int src, int dst, mv_t* hor, mv_t* ver);
void uvg_round_precision_vector2d(int src, int dst, vector2d_t* mv);

void uvg_inter_recon_cu(
  const encoder_state_t * const state,
  lcu_t *lcu,
  bool predict_luma,
  bool predict_chroma,
  const cu_loc_t* const cu_loc);

void uvg_inter_pred_pu(
  const encoder_state_t * const state,
  lcu_t *lcu,
  bool predict_luma,
  bool predict_chroma,
  const cu_loc_t* const cu_loc);

void uvg_hmvp_add_mv(const encoder_state_t* const state, uint32_t pic_x, uint32_t pic_y, uint32_t block_width, uint32_t block_height, const cu_info_t* cu);

void uvg_inter_recon_bipred(
  const encoder_state_t * const state,
  const uvg_picture * ref1,
  const uvg_picture * ref2,
  mv_t mv_param[2][2],
  lcu_t* lcu,
  bool predict_luma,
  bool predict_chroma,
  const cu_loc_t* const cu_loc);


void uvg_inter_get_mv_cand(
  const encoder_state_t * const state,
  mv_t mv_cand[2][2],
  const cu_info_t *const cur_cu,
  lcu_t *lcu,
  int8_t reflist,
  const cu_loc_t* const cu_loc);

void uvg_inter_get_mv_cand_cua(
  const encoder_state_t * const state,
  mv_t mv_cand[2][2],
  const cu_info_t* cur_cu,
  int8_t reflist,
  const cu_loc_t* const cu_loc);

uint8_t uvg_inter_get_merge_cand(
  const encoder_state_t * const state,
  const cu_loc_t* const cu_loc,
  inter_merge_cand_t mv_cand[MRG_MAX_NUM_CANDS],
  lcu_t *lcu);
#endif
