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
#pragma once

#ifdef UVG_DEBUG_PRINT_YUVIEW_CSV

enum {
  DBG_YUVIEW_CU_TYPE          =  0,
  DBG_YUVIEW_INTRADIR_LUMA    =  1,
  DBG_YUVIEW_INTRADIR_CHROMA  =  2,
  DBG_YUVIEW_REFIDX_SKIP_L0   =  3,
  DBG_YUVIEW_REFIDX_SKIP_L1   =  4,
  DBG_YUVIEW_REFIDX_MERGE_L0  =  5,
  DBG_YUVIEW_REFIDX_MERGE_L1  =  6,
  DBG_YUVIEW_REFIDX_INTER_L0  =  7,
  DBG_YUVIEW_REFIDX_INTER_L1  =  8,
  DBG_YUVIEW_MVSKIP_L0        =  9,
  DBG_YUVIEW_MVSKIP_L1        = 10,
  DBG_YUVIEW_MVMERGE_L0       = 11,
  DBG_YUVIEW_MVMERGE_L1       = 12,
  DBG_YUVIEW_MVINTER_L0       = 13,
  DBG_YUVIEW_MVINTER_L1       = 14,
  DBG_YUVIEW_NUM_SIG_COEFF_Y  = 15,
  DBG_YUVIEW_NUM_SIG_COEFF_U  = 16,
  DBG_YUVIEW_NUM_SIG_COEFF_V  = 17,
  DBG_YUVIEW_TR_SKIP          = 18,
  DBG_YUVIEW_MRL              = 19,
  DBG_YUVIEW_NUM_ITEMS        = 20,
  
};

typedef struct encoder_state_t encoder_state_t;
typedef struct encoder_control_t encoder_control_t;

void uvg_dbg_yuview_init(const encoder_control_t* const encoder, char* filename, char* sequence);
void uvg_dbg_yuview_add_vector(int poc, int x, int y, int width, int height, int type, int x_vec, int y_vec);
void uvg_dbg_yuview_add(int poc, int x, int y, int width, int height, int type, int val);
void uvg_dbg_yuview_finish_frame(int poc);
void uvg_dbg_yuview_cleanup();

#define DBG_YUVIEW_INIT(_encoder, _filename, _sequence) uvg_dbg_yuview_init(_encoder, _filename, _sequence);
#define DBG_YUVIEW_MV(_poc, _type, _x, _y, _width, _height, _x_vec, _y_vec) uvg_dbg_yuview_add_vector(_poc, _x, _y, _width, _height, _type, _x_vec, _y_vec);
#define DBG_YUVIEW_VALUE(_poc, _type, _x, _y, _width, _height, _val) uvg_dbg_yuview_add(_poc, _x, _y, _width, _height, _type, _val);
#define DBG_YUVIEW_FINISH_FRAME(_poc) uvg_dbg_yuview_finish_frame(_poc);
#define DBG_YUVIEW_CLEANUP() uvg_dbg_yuview_cleanup();
#else
#define DBG_YUVIEW_INIT(_encoder, _filename, _sequence)
#define DBG_YUVIEW_MV(_poc, _type, _x, _y, _width, _height, _x_vec, _y_vec)
#define DBG_YUVIEW_VALUE(_poc, _type, _x, _y, _width, _height, _val)
#define DBG_YUVIEW_FINISH_FRAME(_poc)
#define DBG_YUVIEW_CLEANUP()
#endif //UVG_DEBUG_PRINT_YUVIEW_CSV

#ifdef UVG_DEBUG_PRINT_THREADING_INFO
void uvg_dbg_encoder_state_dump_graphviz(const encoder_state_t* const state)
#endif //UVG_DEBUG_PRINT_THREADING_INFO

#ifdef UVG_DEBUG_PRINT_MV_INFO
typedef struct encoder_state_t encoder_state_t;

#include <stdint.h>
#include "cu.h"

void uvg_print_merge_vectors(const encoder_state_t* const state, uint32_t pic_x, uint32_t pic_y, uint32_t block_width, uint32_t block_height, cu_info_t* cu);
#define DBG_PRINT_MV(_state, _pic_x, _pic_y, _block_width, _block_height, _cu) uvg_print_merge_vectors(_state, _pic_x, _pic_y, _block_width, _block_height, _cu);
#else
#define DBG_PRINT_MV(_state, _pic_x, _pic_y, _block_width, _block_height, _cu)
#endif // UVG_DEBUG_PRINT_MV_INFO