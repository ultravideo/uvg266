#ifndef SEARCH_INTRA_H_
#define SEARCH_INTRA_H_
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
 * \ingroup Compression
 * \file
 * Intra prediction parameter search.
 */

#include "cu.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "intra.h"

double uvg_luma_mode_bits(const encoder_state_t *state, const cu_info_t* const cur_cu, const cu_loc_t*
                          const cu_loc,
                          const lcu_t* lcu);
                       
double uvg_chroma_mode_bits(const encoder_state_t *state,
                        int8_t chroma_mode, int8_t luma_mode);

int8_t uvg_search_cu_intra_chroma(
  encoder_state_t * const state,
  const cu_loc_t* const cu_loc,
  lcu_t *lcu,
  intra_search_data_t* best_cclm,
  int8_t luma_mode,
  enum uvg_tree_type tree_type,
  bool is_separate);

void uvg_search_cu_intra(
  encoder_state_t * const state,
  intra_search_data_t* search_data,
  lcu_t *lcu,
  enum uvg_tree_type tree_type,
  const cu_loc_t* const cu_loc);

#endif // SEARCH_INTRA_H_
