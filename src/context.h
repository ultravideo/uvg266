#ifndef CONTEXT_H_
#define CONTEXT_H_
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
 * \ingroup CABAC
 * \file
 * Context derivation for CABAC.
 */

#include "cabac.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep


// Functions
void uvg_ctx_init(cabac_ctx_t* ctx, int32_t qp, int32_t init_value, uint8_t rate);
void uvg_init_contexts(encoder_state_t *state, int8_t QP, int8_t slice);

void uvg_context_copy(encoder_state_t *const target_state, const encoder_state_t *const source_state);

uint32_t uvg_context_get_sig_coeff_group( uint32_t *sig_coeff_group_flag,uint32_t pos_x, uint32_t pos_y,int32_t width, int32_t height);
uint32_t uvg_context_get_sig_coeff_group_ts(uint32_t* sig_coeff_group_flag, uint32_t pos_x, uint32_t pos_y, int32_t width);
uint32_t uvg_context_get_sig_ctx_idx_abs(const coeff_t* coeff, uint32_t pos_x, uint32_t pos_y,
                                         uint32_t width, uint32_t height, int8_t type, 
                                         int32_t* temp_diag, int32_t* temp_sum);

uint32_t uvg_context_get_sig_ctx_idx_abs_ts(const coeff_t* coeff, int32_t pos_x, int32_t pos_y,
                                             uint32_t width);
uint32_t uvg_sign_ctx_id_abs_ts(const coeff_t* coeff, int32_t pos_x, int32_t pos_y, int32_t width, int bdpcm);
int32_t uvg_derive_mod_coeff(int rightPixel, int belowPixel, coeff_t absCoeff, int bdpcm);
unsigned uvg_lrg1_ctx_id_abs_ts(const coeff_t* coeff, int32_t pos_x, int32_t pos_y, int32_t width, int bdpcm);


uint32_t uvg_abs_sum(const coeff_t* coeff, uint32_t pos_x, uint32_t pos_y,
                     uint32_t height, uint32_t width, uint32_t baselevel);

uint32_t uvg_go_rice_par_abs(const coeff_t* coeff, uint32_t pos_x, uint32_t pos_y,
                             uint32_t width, uint32_t height, uint32_t baselevel);

#define CNU 35
#define DWS 8

#endif
