#ifndef ENCODING_RESUME_H_
#define ENCODING_RESUME_H_
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
 * \file
 * Allow skipping search for specific frames to speed up debugging by reading previous data from file
 */

#include "global.h" // IWYU pragma: keep

#ifdef UVG_ENCODING_RESUME

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "cu.h"
#include "encoderstate.h"
#include "uvg266.h"

#define RESUME_DIRNAME "./resume_data"

#define RESUME_SLICE_COND UVG_SLICE_I
//#define RESUME_POC_LTE_COND 0
//#define RESUME_X_LTE_COND 0
//#define RESUME_Y_LTE_COND 0
//#define RESUME_NO_CHROMA_COND 1

#define RESUME_SUB_FRAME_POC_COND 16
#define RESUME_SUB_FRAME_LCU_IND_LT_COND 791


bool uvg_can_resume_encoding(const encoder_state_t * const state, const int x, const int y, const bool chroma_only);
void uvg_process_resume_encoding(const encoder_state_t * const state, const int x, const int y, const bool chroma_only, double * const cost, lcu_t* const lcu, bool read_mode);

#endif

#if !defined(UVG_ENCODING_RESUME)

#endif



#endif //ENCODER_RESUME_H_
