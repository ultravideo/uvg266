#ifndef CONTEXT_H_
#define CONTEXT_H_
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2015 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * Kvazaar is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
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
void kvz_ctx_init(cabac_ctx_t* ctx, int32_t qp, int32_t init_value, uint8_t rate);
void kvz_init_contexts(encoder_state_t *state, int8_t QP, int8_t slice);

void kvz_context_copy(encoder_state_t * target_state, const encoder_state_t * source_state);

uint32_t kvz_context_get_sig_coeff_group( uint32_t *sig_coeff_group_flag,uint32_t pos_x, uint32_t pos_y,int32_t width);

uint32_t kvz_context_get_sig_ctx_idx_abs(const coeff_t* coeff, int32_t pos_x, int32_t pos_y,
                                         uint32_t height, uint32_t width, int8_t type, 
                                         int32_t* temp_diag, int32_t* temp_sum);

uint32_t kvz_context_get_sig_ctx_idx_abs_ts(const coeff_t* coeff, int32_t pos_x, int32_t pos_y,
                                             uint32_t width);
uint32_t kvz_sign_ctx_id_abs_ts(const coeff_t* coeff, int32_t pos_x, int32_t pos_y, int32_t width, int bdpcm);
int32_t kvz_derive_mod_coeff(int rightPixel, int belowPixel, coeff_t absCoeff, int bdpcm);
unsigned kvz_lrg1_ctx_id_abs_ts(const coeff_t* coeff, int32_t pos_x, int32_t pos_y, int32_t width, int bdpcm);


uint32_t kvz_abs_sum(const coeff_t* coeff, int32_t pos_x, int32_t pos_y,
                     uint32_t height, uint32_t width, uint32_t baselevel);

uint32_t kvz_go_rice_par_abs(const coeff_t* coeff, int32_t pos_x, int32_t pos_y,
                             uint32_t height, uint32_t width, uint32_t baselevel);

#define CNU 35
#define DWS 8

#endif
