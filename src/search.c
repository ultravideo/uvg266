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

#include "search.h"

#include <limits.h>
#include <string.h>

#include "cabac.h"
#include "cu.h"
#include "encoder.h"
#include "encode_coding_tree.h"
#include "filter.h"
#include "imagelist.h"
#include "inter.h"
#include "intra.h"
#include "rate_control.h"
#include "uvg266.h"
#include "rdo.h"
#include "search_inter.h"
#include "search_intra.h"
#include "search_ibc.h"
#include "threadqueue.h"
#include "transform.h"
#include "videoframe.h"
#include "strategies/strategies-picture.h"
#include "strategies/strategies-quant.h"
#include "reshape.h"

#ifdef UVG_ENCODING_RESUME
#include "encoding_resume.h"
#endif // UVG_ENCODING_RESUME


#define IN_FRAME(x, y, width, height, block_width, block_height) \
  ((x) >= 0 && (y) >= 0 \
  && (x) + (block_width) <= (width) \
  && (y) + (block_height) <= (height))

// Cost threshold for doing intra search in inter frames with --rd=0.
static const int INTRA_THRESHOLD = 8;


static INLINE void copy_cu_info(lcu_t *from, lcu_t *to, const cu_loc_t* const cu_loc, enum uvg_tree_type
                                tree_type)
{
  const int y_limit = (cu_loc->local_y + cu_loc->height);
  const int x_limit = (cu_loc->local_x + cu_loc->width);
  for   (int y = cu_loc->local_y ; y < y_limit; y += SCU_WIDTH) {
    for (int x = cu_loc->local_x ; x < x_limit; x += SCU_WIDTH) {
      *LCU_GET_CU_AT_PX(to, x, y) = *LCU_GET_CU_AT_PX(from, x, y);
    }
  }
}


static INLINE void initialize_partial_work_tree(
  const encoder_state_t* const state,
  lcu_t* from,
  lcu_t *to,
  const cu_loc_t * const cu_loc,
  const cu_loc_t* const
  chroma_loc,
  const enum uvg_tree_type tree_type) {

  const int y_limit = MIN(LCU_WIDTH,  state->tile->frame->height - cu_loc->y / 64 * 64);
  const int x_limit = MIN(LCU_WIDTH, state->tile->frame->width - cu_loc->x / 64 * 64);

  if (cu_loc->local_x == 0) {
    to->left_ref = from->left_ref;
    *LCU_GET_TOP_RIGHT_CU(to) = *LCU_GET_TOP_RIGHT_CU(from);
  }
  else {
    if(tree_type != UVG_CHROMA_T) {
      uvg_pixels_blit(from->rec.y, to->rec.y, cu_loc->local_x, LCU_WIDTH, LCU_WIDTH, LCU_WIDTH);
    }
    if(tree_type != UVG_LUMA_T && from->ref.chroma_format != UVG_CSP_400) {
      uvg_pixels_blit(from->rec.u, to->rec.u, chroma_loc->local_x / 2, LCU_WIDTH_C, LCU_WIDTH_C, LCU_WIDTH_C);
      uvg_pixels_blit(from->rec.v, to->rec.v, chroma_loc->local_x / 2, LCU_WIDTH_C, LCU_WIDTH_C, LCU_WIDTH_C);
    }
  }

  if (cu_loc->local_y == 0) {
    to->top_ref = from->top_ref;
    *LCU_GET_TOP_RIGHT_CU(to) = *LCU_GET_TOP_RIGHT_CU(from);
  }
  else {
    if (tree_type != UVG_CHROMA_T) {
      uvg_pixels_blit(&from->rec.y[cu_loc->local_x], &to->rec.y[cu_loc->local_x], 
        LCU_WIDTH - cu_loc->local_x, cu_loc->local_y,
        LCU_WIDTH, LCU_WIDTH);
    }
    if (tree_type != UVG_LUMA_T && from->ref.chroma_format != UVG_CSP_400) {
      uvg_pixels_blit(&from->rec.u[chroma_loc->local_x / 2], &to->rec.u[chroma_loc->local_x / 2],
        LCU_WIDTH_C - chroma_loc->local_x / 2, chroma_loc->local_y / 2,
        LCU_WIDTH_C, LCU_WIDTH_C);
      uvg_pixels_blit(&from->rec.v[chroma_loc->local_x / 2], &to->rec.v[chroma_loc->local_x / 2],
        LCU_WIDTH_C - chroma_loc->local_x / 2, chroma_loc->local_y / 2,
        LCU_WIDTH_C, LCU_WIDTH_C);
    }
  }

  if (tree_type == UVG_CHROMA_T) {
    // These are needed for CCLM
    uvg_pixels_blit(from->rec.y, to->rec.y, MIN(cu_loc->local_x + cu_loc->width * 2, LCU_WIDTH), MIN(cu_loc->local_y + cu_loc->height * 2, LCU_WIDTH), LCU_WIDTH, LCU_WIDTH);
  }

  to->ref.chroma_format = from->ref.chroma_format;
  to->rec.chroma_format = from->rec.chroma_format;

  if (tree_type != UVG_CHROMA_T) {
    const int offset = cu_loc->local_x + cu_loc->local_y * LCU_WIDTH;
    uvg_pixels_blit(&from->ref.y[offset], &to->ref.y[offset], cu_loc->width, cu_loc->height, LCU_WIDTH, LCU_WIDTH);
  }

  if(tree_type != UVG_LUMA_T && from->ref.chroma_format != UVG_CSP_400) {
    const int offset = chroma_loc->local_x / 2 + chroma_loc->local_y / 2 * LCU_WIDTH_C;
    uvg_pixels_blit(&from->ref.u[offset], &to->ref.u[offset], chroma_loc->chroma_width, chroma_loc->chroma_height, LCU_WIDTH_C, LCU_WIDTH_C);
    uvg_pixels_blit(&from->ref.v[offset], &to->ref.v[offset], chroma_loc->chroma_width, chroma_loc->chroma_height, LCU_WIDTH_C, LCU_WIDTH_C);
  }
  if(chroma_loc->local_y != cu_loc->local_y || chroma_loc->local_x != cu_loc->local_x && tree_type == UVG_BOTH_T) {
    for (int y = chroma_loc->local_y; y < chroma_loc->local_y + chroma_loc->height; y += SCU_WIDTH) {
      for (int x = chroma_loc->local_x; x < chroma_loc->local_x + chroma_loc->width; x += SCU_WIDTH) {
        memset(LCU_GET_CU_AT_PX(to, x, y), 0, sizeof(cu_info_t));
      }
    }
    
  }

  const int y_start = (cu_loc->local_y) - 4;
  const int x_start = (cu_loc->local_x) - 4;
  for (int y = y_start; y < y_limit; y += SCU_WIDTH) {
    *LCU_GET_CU_AT_PX(to, x_start, y) = *LCU_GET_CU_AT_PX(from, x_start, y);
  }
  for (int x = x_start; x < x_limit; x += SCU_WIDTH) {
    *LCU_GET_CU_AT_PX(to, x, y_start) = *LCU_GET_CU_AT_PX(from, x, y_start);
  }

  for (int y = cu_loc->local_y; y < y_limit; y += SCU_WIDTH) {
    for (int x = cu_loc->local_x ; x < x_limit; x += SCU_WIDTH) {
      memset(LCU_GET_CU_AT_PX(to, x, y), 0, sizeof(cu_info_t));
    }
  }

  if(chroma_loc->local_y != cu_loc->local_y || chroma_loc->local_x != cu_loc->local_x && tree_type == UVG_BOTH_T) {
    const int y_start = (chroma_loc->local_y) - 4;
    const int x_start = (chroma_loc->local_x) - 4;
    for (int y = y_start; y < y_limit; y += SCU_WIDTH) {
      *LCU_GET_CU_AT_PX(to, x_start, y) = *LCU_GET_CU_AT_PX(from, x_start, y);
    }
    for (int x = x_start; x < y_limit; x += SCU_WIDTH) {
      *LCU_GET_CU_AT_PX(to, x, y_start) = *LCU_GET_CU_AT_PX(from, x, y_start);
    }

    for(int y = chroma_loc->local_y; y < chroma_loc->local_y + chroma_loc->height; y += SCU_WIDTH) {
      for (int x = chroma_loc->local_x; x < chroma_loc->local_x + chroma_loc->width; x += SCU_WIDTH) {
        if(x >= cu_loc->local_x && y>= cu_loc->local_y) continue;
        *LCU_GET_CU_AT_PX(to, x, y) = *LCU_GET_CU_AT_PX(from, x, y);
      }      
    }

    if (chroma_loc->local_x == 0) {
      to->left_ref = from->left_ref;
      *LCU_GET_TOP_RIGHT_CU(to) = *LCU_GET_TOP_RIGHT_CU(from);      
    }
    if (chroma_loc->local_y == 0) {
      to->top_ref = from->top_ref;
      *LCU_GET_TOP_RIGHT_CU(to) = *LCU_GET_TOP_RIGHT_CU(from);      
    }
    if (x_limit != LCU_WIDTH) {
      for (int y = y_start; y < y_limit; y += SCU_WIDTH) {
        memset(LCU_GET_CU_AT_PX(to, x_limit, y), 0, sizeof(cu_info_t));
      }
    }
    if (y_limit != LCU_WIDTH) {
      for (int x = x_start; x < x_limit; x += SCU_WIDTH) {
        memset(LCU_GET_CU_AT_PX(to, x, y_limit), 0, sizeof(cu_info_t));
      }
    }
  }
  else {
    if (x_limit != LCU_WIDTH) {
      for (int y = y_start; y < y_limit; y += SCU_WIDTH) {
        memset(LCU_GET_CU_AT_PX(to, x_limit, y), 0, sizeof(cu_info_t));
      }
    }
    if (y_limit != LCU_WIDTH) {
      for (int x = x_start; x < x_limit; x += SCU_WIDTH) {
        memset(LCU_GET_CU_AT_PX(to, x, y_limit), 0, sizeof(cu_info_t));
      }
    }
  }
}

static INLINE void copy_cu_pixels(
  lcu_t *from,
  lcu_t *to,
  const cu_loc_t* const cu_loc,
  enum uvg_tree_type
  tree_type)
{
  const int x_local = cu_loc->local_x;
  const int y_local = cu_loc->local_y;
  const int luma_index = x_local + y_local * LCU_WIDTH;
  const int chroma_index =  (x_local / 2) + (y_local / 2) * LCU_WIDTH_C;

  if(tree_type != UVG_CHROMA_T) {
    uvg_pixels_blit(&from->rec.y[luma_index], &to->rec.y[luma_index],
                    cu_loc->width, cu_loc->height, LCU_WIDTH, LCU_WIDTH);
  }
  if (from->rec.chroma_format != UVG_CSP_400 && tree_type != UVG_LUMA_T) {
    uvg_pixels_blit(&from->rec.u[chroma_index], &to->rec.u[chroma_index],
                    cu_loc->chroma_width, cu_loc->chroma_height, LCU_WIDTH_C, LCU_WIDTH_C);
    uvg_pixels_blit(&from->rec.v[chroma_index], &to->rec.v[chroma_index],
                    cu_loc->chroma_width, cu_loc->chroma_height, LCU_WIDTH_C, LCU_WIDTH_C);
  }
}

// ISP_TODO: this needs to work with the new coeff cu orderr
static INLINE void copy_cu_coeffs(const cu_loc_t *cu_loc, lcu_t *from, lcu_t *to, bool joint, enum
                                  uvg_tree_type tree_type)
{
  if (tree_type != UVG_CHROMA_T) {
    //const int luma_z = xy_to_zorder(LCU_WIDTH, cu_loc->x, cu_loc->y);
    const int idx = (cu_loc->x % LCU_WIDTH) + ((cu_loc->y % LCU_WIDTH) * LCU_WIDTH);
    copy_coeffs(&from->coeff.y[idx], &to->coeff.y[idx], cu_loc->width, cu_loc->height, LCU_WIDTH);
    
  }

  if (from->rec.chroma_format != UVG_CSP_400 && tree_type != UVG_LUMA_T) {
    //const int chroma_z = xy_to_zorder(LCU_WIDTH_C, cu_loc->x >> (tree_type != UVG_CHROMA_T), cu_loc->y >> (tree_type != UVG_CHROMA_T));
    const int chroma_x = (cu_loc->x >> 1);
    const int chroma_y = (cu_loc->y >> 1);

    const int idx = (chroma_x % LCU_WIDTH_C) + ((chroma_y % LCU_WIDTH_C) * LCU_WIDTH_C);
    copy_coeffs(&from->coeff.u[idx], &to->coeff.u[idx], cu_loc->chroma_width, cu_loc->chroma_height, LCU_WIDTH_C);
    copy_coeffs(&from->coeff.v[idx], &to->coeff.v[idx], cu_loc->chroma_width, cu_loc->chroma_height, LCU_WIDTH_C);
    if (joint) {
      copy_coeffs(&from->coeff.joint_uv[idx], &to->coeff.joint_uv[idx], cu_loc->chroma_width, cu_loc->chroma_height, LCU_WIDTH_C);
    }
  }
}


static void lcu_fill_chroma_cu_info(lcu_t* lcu, const cu_loc_t* const cu_loc);
/**
 * Copy all non-reference CU data from next level to current level.
 */
static void work_tree_copy_up(
  lcu_t *from,
  lcu_t* to,
  bool joint,
  enum
  uvg_tree_type tree_type,
  const cu_loc_t* const cu_loc,
  const cu_loc_t* const chroma_loc)
{
  copy_cu_info  (from, to, cu_loc, tree_type);
  copy_cu_pixels(from, to, cu_loc, cu_loc != chroma_loc && tree_type == UVG_LUMA_T ? UVG_LUMA_T : tree_type);
  copy_cu_coeffs(cu_loc, from, to, joint, cu_loc != chroma_loc && tree_type == UVG_LUMA_T ? UVG_LUMA_T : tree_type);
  if (chroma_loc && tree_type != UVG_LUMA_T) {
    copy_cu_pixels(from, to, chroma_loc, UVG_CHROMA_T);
    copy_cu_coeffs(chroma_loc, from, to, joint, UVG_CHROMA_T);

    for (int y = chroma_loc->local_y; y < chroma_loc->local_y + chroma_loc->height; y += 4) {
      for (int x = chroma_loc->local_x; x < chroma_loc->local_x + chroma_loc->width; x += 4) {
        cu_info_t* to_cu = LCU_GET_CU_AT_PX(to, x, y);
        cu_info_t* from_cu = LCU_GET_CU_AT_PX(from, x, y);
        to_cu->intra.mode_chroma = from_cu->intra.mode_chroma;
        to_cu->joint_cb_cr = from_cu->joint_cb_cr;
        to_cu->cr_lfnst_idx = from_cu->cr_lfnst_idx;
        to_cu->chroma_deblocking = from_cu->chroma_deblocking;
        to_cu->log2_chroma_width = from_cu->log2_chroma_width;
        to_cu->log2_chroma_height = from_cu->log2_chroma_height;

        cbf_copy(&to_cu->cbf, from_cu->cbf, COLOR_U);
        cbf_copy(&to_cu->cbf, from_cu->cbf, COLOR_V);
      }
    }
  }
  
}


static void lcu_fill_cu_info(lcu_t *lcu, int x_local, int y_local, int width, int height, const cu_info_t *cu)
{
  // Set mode in every CU covered by part_mode in this depth.
  for (int y = y_local; y < y_local + height; y += SCU_WIDTH) {
    for (int x = x_local; x < x_local + width; x += SCU_WIDTH) {
      cu_info_t *to = LCU_GET_CU_AT_PX(lcu, x, y);
      to->type      = cu->type;
      to->qp        = cu->qp;
      to->split_tree = cu->split_tree;
      to->mode_type_tree = cu->mode_type_tree;
      //to->tr_idx    = cu->tr_idx;
      to->lfnst_idx = cu->lfnst_idx;
      to->cr_lfnst_idx = cu->cr_lfnst_idx;
      to->joint_cb_cr = cu->joint_cb_cr;
      to->lfnst_last_scan_pos = cu->lfnst_last_scan_pos;
      to->violates_lfnst_constrained_luma = cu->violates_lfnst_constrained_luma;
      to->violates_lfnst_constrained_chroma = cu->violates_lfnst_constrained_chroma;

      to->log2_height = cu->log2_height;
      to->log2_width = cu->log2_width;

      to->log2_chroma_height = cu->log2_chroma_height;
      to->log2_chroma_width = cu->log2_chroma_width;

      if (cu->type == CU_INTRA) {
        to->intra.mode        = cu->intra.mode;
        to->intra.mode_chroma = cu->intra.mode_chroma;
        to->intra.multi_ref_idx = cu->intra.multi_ref_idx;
        to->intra.mip_flag = cu->intra.mip_flag;
        to->intra.mip_is_transposed = cu->intra.mip_is_transposed;
        to->intra.isp_mode = cu->intra.isp_mode;
      } else {
        to->skipped   = cu->skipped;
        to->merged    = cu->merged;
        to->merge_idx = cu->merge_idx;
        to->inter     = cu->inter;
      }
    }
  }
}

static void lcu_fill_chroma_cu_info(lcu_t *lcu, const cu_loc_t * const cu_loc)
{
  // The bottom right cu will always have the chroma info
  cu_info_t *bottom_right = LCU_GET_CU_AT_PX(
    lcu,
    cu_loc->local_x + cu_loc->width - 1,
    cu_loc->local_y + cu_loc->height - 1);
  if(bottom_right->type != CU_INTRA) return;


  for(int y = cu_loc->local_y; y < cu_loc->local_y + cu_loc->height; y += 4 ) {
    for (int x = cu_loc->local_x; x < cu_loc->local_x + cu_loc->width; x += 4) {
      cu_info_t *cu         = LCU_GET_CU_AT_PX(lcu, x, y);
      cu->intra.mode_chroma = bottom_right->intra.mode_chroma;
      cu->joint_cb_cr       = bottom_right->joint_cb_cr;
      cu->cr_lfnst_idx      = bottom_right->cr_lfnst_idx;
      cu->log2_chroma_height = bottom_right->log2_chroma_height;
      cu->log2_chroma_width = bottom_right->log2_chroma_width;
      cu->type = bottom_right->type;
      cu->tr_skip |= bottom_right->tr_skip & 6;
    }
  }
}


static void lcu_fill_chroma_cbfs(lcu_t *lcu, const cu_loc_t * const chroma_loc, enum uvg_tree_type tree_type)
{
  uint8_t height = chroma_loc->height;
  uint8_t width =  chroma_loc->width;
  uint32_t x_local = chroma_loc->local_x;
  uint32_t y_local = chroma_loc->local_y;
  const int offset = ~((TR_MAX_WIDTH) - 1);
  // Set coeff flags in every CU covered by part_mode in this depth.
  for (uint32_t y = 0; y < height; y += SCU_WIDTH) {
    for (uint32_t x = 0; x < width; x += SCU_WIDTH) {
      // Use TU top-left CU to propagate coeff flags
      cu_info_t* cu_from = LCU_GET_CU_AT_PX(lcu, x_local + (x & offset), y_local + (y & offset));
      cu_info_t* cu_to = LCU_GET_CU_AT_PX(lcu, x_local + x, y_local + y);
      if (cu_from != cu_to) {
        cbf_copy(&cu_to->cbf, cu_from->cbf, COLOR_U);
        cbf_copy(&cu_to->cbf, cu_from->cbf, COLOR_V);
      }
    }
  }
  
}

static void lcu_fill_cbf(lcu_t *lcu, int x_local, unsigned y_local, unsigned width, unsigned height, const cu_info_t *cur_cu, enum
                         uvg_tree_type tree_type)
{
  // Set coeff flags in every CU covered by part_mode in this depth.
  for (uint32_t y = 0; y < height; y += SCU_WIDTH) {
    for (uint32_t x = 0; x < width; x += SCU_WIDTH) {
      // Use TU top-left CU to propagate coeff flags
      cu_info_t *cu_from = LCU_GET_CU_AT_PX(lcu, x_local + (x & ~(TR_MAX_WIDTH - 1)), y_local + (y & ~(TR_MAX_WIDTH - 1)));
      cu_info_t *cu_to   = LCU_GET_CU_AT_PX(lcu, x_local + x, y_local + y);
      if (cu_from != cu_to) {
        // Chroma and luma coeff data is needed for deblocking
        if(tree_type != UVG_CHROMA_T) cbf_copy(&cu_to->cbf, cu_from->cbf, COLOR_Y);
        if(tree_type != UVG_LUMA_T) cbf_copy(&cu_to->cbf, cu_from->cbf, COLOR_U);
        if (tree_type != UVG_LUMA_T)cbf_copy(&cu_to->cbf, cu_from->cbf, COLOR_V);
      }
    }
  }
}


//Calculates cost for all zero coeffs
static double cu_zero_coeff_cost(
  const encoder_state_t *state,
  lcu_t *work_tree,
  const cu_loc_t* const cu_loc,
  const int depth)
{
  lcu_t *const lcu = &work_tree[depth];

  const int y_local = cu_loc->local_y;
  const int x_local = cu_loc->local_x;

  const int luma_index = y_local * LCU_WIDTH + x_local;
  const int chroma_index = (y_local / 2) * LCU_WIDTH_C + (x_local / 2);

  double ssd = 0.0;
  ssd += UVG_LUMA_MULT * uvg_pixels_calc_ssd(
    &lcu->ref.y[luma_index], &lcu->rec.y[luma_index],
    LCU_WIDTH, LCU_WIDTH, cu_loc->width, cu_loc->height
    );
  if (y_local % 8 == 0 && x_local % 8 == 0 && state->encoder_control->chroma_format != UVG_CSP_400) {
    ssd += UVG_CHROMA_MULT * uvg_pixels_calc_ssd(
      &lcu->ref.u[chroma_index], &lcu->rec.u[chroma_index],
      LCU_WIDTH_C, LCU_WIDTH_C, cu_loc->chroma_width, cu_loc->chroma_height
      );
    ssd += UVG_CHROMA_MULT * uvg_pixels_calc_ssd(
      &lcu->ref.v[chroma_index], &lcu->rec.v[chroma_index],
      LCU_WIDTH_C, LCU_WIDTH_C, cu_loc->chroma_width, cu_loc->chroma_height
      );
  }
  // Save the pixels at a lower level of the working tree.
  copy_cu_pixels(lcu, &work_tree[depth + 1], cu_loc, UVG_BOTH_T);

  return ssd;
}


static void downsample_cclm_rec(encoder_state_t *state, int x, int y, int width, int height, uvg_pixel *y_rec, uvg_pixel extra_pixel) {
  if (!state->encoder_control->cfg.cclm) return;
  int x_scu = SUB_SCU(x);
  int y_scu = SUB_SCU(y);
  y_rec += x_scu + y_scu * LCU_WIDTH;
  const int stride = state->tile->frame->rec->stride;
  const int stride2 = (((state->tile->frame->width + 7) & ~7) + FRAME_PADDING_LUMA);

  for (int y_ = 0; y_ < height && y_ * 2 + y < state->tile->frame->height; y_++) {
    for (int x_ = 0; x_ < width; x_++) {
      int s = 4;
      s += y_rec[2 * x_] * 2;
      s += y_rec[2 * x_ + 1];
      // If we are at the edge of the CTU read the pixel from the frame reconstruct buffer,
      // *except* when we are also at the edge of the frame, in which case we want to duplicate
      // the edge pixel
      s += !x_scu && !x_ && x ? state->tile->frame->rec->y[x - 1 + (y + y_ * 2) * stride] : y_rec[2 * x_ - ((x_ + x) > 0)];
      s += y_rec[2 * x_ + LCU_WIDTH] * 2;
      s += y_rec[2 * x_ + 1 + LCU_WIDTH];
      s += !x_scu && !x_ && x ? state->tile->frame->rec->y[x - 1 + (y + y_ * 2 + 1) * stride] : y_rec[2 * x_ - ((x_ + x) > 0) + LCU_WIDTH];
      int index = x / 2 + x_ + (y / 2 + y_ )* stride2 / 2;
      state->tile->frame->cclm_luma_rec[index] = s >> 3;
    }
    y_rec += LCU_WIDTH * 2;
  }
  if((y + height * 2) % 64 == 0) {
    int line = y / 64 * stride2 / 2;
    y_rec -= LCU_WIDTH;
    for (int i = 0; i < width && i + x / 2 < stride2 / 2; ++i) {
      int s = 2;
      s += y_rec[i * 2] * 2;
      s += y_rec[i * 2 + 1];
      s += !x_scu && !i && x ? extra_pixel : y_rec[i * 2 - ((i + x) > 0)] ;
      state->tile->frame->cclm_luma_rec_top_line[i + x / 2 + line] = s >> 2;
    }
  }
}


/**
* Calculate RD cost for a Coding Unit.
* \return Cost of block
* \param ref_cu  CU used for prediction parameters.
*
* Calculates the RDO cost of a single CU that will not be split further.
* Takes into account SSD of reconstruction and the cost of encoding whatever
* prediction unit data needs to be coded.
*/
double uvg_cu_rd_cost_luma(
  const encoder_state_t *const state,
  const cu_loc_t* const cu_loc,
  const cu_info_t *const pred_cu,
  lcu_t *const lcu,
  uint8_t isp_cbf)
{
  const int skip_residual_coding = pred_cu->skipped || (pred_cu->type != CU_INTRA && pred_cu->cbf == 0);
  cabac_data_t* cabac = (cabac_data_t *)&state->search_cabac;
  
  // cur_cu is used for TU parameters.
  cu_info_t *const tr_cu = LCU_GET_CU_AT_PX(lcu, cu_loc->local_x, cu_loc->local_y);

  double coeff_bits = 0;
  double tr_tree_bits = 0;

  // Check that lcu is not in   

  if (cu_loc->width > TR_MAX_WIDTH || cu_loc->height > TR_MAX_WIDTH) {
    double sum = 0;
    // Recursively process sub-CUs.
    enum split_type split;
    if (cu_loc->width > TR_MAX_WIDTH && cu_loc->height > TR_MAX_WIDTH) {
      split = QT_SPLIT;
    }
    else if (cu_loc->width > TR_MAX_WIDTH) {
      split = BT_VER_SPLIT;
    }
    else {
      split = BT_HOR_SPLIT;
    }

    cu_loc_t split_cu_loc[4];
    const int split_count = uvg_get_split_locs(cu_loc, split, split_cu_loc,NULL);
    for (int i = 0; i < split_count; ++i) {
      sum += uvg_cu_rd_cost_luma(state, &split_cu_loc[i], pred_cu, lcu, isp_cbf);
    }

    return sum + tr_tree_bits * state->lambda;
  }

  const bool is_not_isp = pred_cu->type == CU_INTER || pred_cu->intra.isp_mode == ISP_MODE_NO_ISP;
  // Add transform_tree cbf_luma bit cost.
  if (is_not_isp) {
    const int depth = 6 - uvg_g_convert_to_log2[cu_loc->width];
    int is_set = cbf_is_set(pred_cu->cbf, COLOR_Y);
    if (pred_cu->type == CU_INTRA ||
      !PU_IS_TU(pred_cu) ||
      cbf_is_set(tr_cu->cbf, COLOR_U) ||
      cbf_is_set(tr_cu->cbf, COLOR_V))
    {
      cabac_ctx_t* ctx = &(cabac->ctx.qt_cbf_model_luma[0]);

      CABAC_FBITS_UPDATE(cabac, ctx, is_set, tr_tree_bits, "cbf_y_search");
    }

    if (is_set && state->encoder_control->cfg.trskip_enable 
      && cu_loc->width <= (1 << state->encoder_control->cfg.trskip_max_size)
      && cu_loc->height <= (1 << state->encoder_control->cfg.trskip_max_size)) {
      CABAC_FBITS_UPDATE(cabac, &cabac->ctx.transform_skip_model_luma, pred_cu->tr_idx == MTS_SKIP, tr_tree_bits, "transform_skip_flag");
    }
  }
  else {
    // TODO: 8x4 CUs
    const int split_limit = uvg_get_isp_split_num(cu_loc->width, cu_loc->height, pred_cu->intra.isp_mode, true);
    int luma_ctx = 2;
    const int split_limit_minus_one = split_limit - 1;
    for (int i = 0; i < split_limit; i++) {
      if (i != split_limit_minus_one || isp_cbf != 1 << split_limit_minus_one) {
        const int flag = (isp_cbf >> i) & 1;
        CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.qt_cbf_model_luma[luma_ctx]), flag, tr_tree_bits, "cbf_y_search");
        luma_ctx = 2 + flag;
      }
    }
  }

  // SSD between reconstruction and original
  int ssd = 0;
  if (!state->encoder_control->cfg.lossless) {
    int index = cu_loc->local_y * LCU_WIDTH + cu_loc->local_x;
    ssd = uvg_pixels_calc_ssd(&lcu->ref.y[index], &lcu->rec.y[index],
                                        LCU_WIDTH,          LCU_WIDTH,
                                        cu_loc->width, cu_loc->height);
  }


  if (!skip_residual_coding) {
    int8_t luma_scan_mode = SCAN_DIAG;
    if (is_not_isp) {
      //const coeff_t* coeffs = &lcu->coeff.y[xy_to_zorder(LCU_WIDTH, x_px, y_px)];
      const coeff_t* coeffs = lcu->coeff.y;

      coeff_bits += uvg_get_coeff_cost(state, coeffs, NULL, cu_loc, 0, luma_scan_mode, pred_cu->tr_idx == MTS_SKIP, COEFF_ORDER_CU);
    }
    else {
      int split_type = pred_cu->intra.isp_mode;
      int split_limit = uvg_get_isp_split_num(cu_loc->width, cu_loc->height, split_type, true);

      for (int i = 0; i < split_limit; ++i) {
        cu_loc_t split_loc;
        uvg_get_isp_split_loc(&split_loc, cu_loc->x, cu_loc->y,  cu_loc->width, cu_loc->height, i, split_type, true);
        const int part_x = split_loc.x;
        const int part_y = split_loc.y;

        // TODO: maybe just pass the cu_loc_t to these functions
        //const coeff_t* coeffs = &lcu->coeff.y[xy_to_zorder(LCU_WIDTH, part_x, part_y)];
        const coeff_t* coeffs = lcu->coeff.y;

        coeff_bits += uvg_get_coeff_cost(state, coeffs, NULL, &split_loc, 0, luma_scan_mode, pred_cu->tr_idx == MTS_SKIP, COEFF_ORDER_CU);
      }
    }
  }

  double bits = tr_tree_bits + coeff_bits;
  return (double)ssd * UVG_LUMA_MULT + bits * state->lambda;
}


double uvg_cu_rd_cost_chroma(
  const encoder_state_t *const state,
  cu_info_t *const pred_cu,
  lcu_t *const lcu,
  const cu_loc_t * const cu_loc)
{
  const vector2d_t lcu_px = { (cu_loc->local_x) / 2, (cu_loc->local_y) / 2 };
  cu_info_t *const tr_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x, lcu_px.y);
  const int skip_residual_coding = pred_cu->skipped || (pred_cu->type != CU_INTRA && pred_cu->cbf == 0);
  
  double tr_tree_bits = 0;
  double coeff_bits = 0;
  
  const int depth = 6 - uvg_g_convert_to_log2[cu_loc->width];
  int u_is_set = pred_cu->joint_cb_cr ? (pred_cu->joint_cb_cr & 2) >> 1 : cbf_is_set(pred_cu->cbf, COLOR_U);
  int v_is_set = pred_cu->joint_cb_cr ? (pred_cu->joint_cb_cr & 1) : cbf_is_set(pred_cu->cbf, COLOR_V);

  if (cu_loc->width > TR_MAX_WIDTH || cu_loc->height > TR_MAX_WIDTH) {
    double sum = 0;
    // Recursively process sub-CUs.
    enum split_type split;
    if (cu_loc->width > TR_MAX_WIDTH && cu_loc->height > TR_MAX_WIDTH) {
      split = QT_SPLIT;
    }
    else if (cu_loc->width > TR_MAX_WIDTH) {
      split = BT_VER_SPLIT;
    }
    else {
      split = BT_HOR_SPLIT;
    }

    cu_loc_t split_cu_loc[4];
    const int split_count = uvg_get_split_locs(cu_loc, split, split_cu_loc,NULL);
    for (int i = 0; i < split_count; ++i) {
      sum += uvg_cu_rd_cost_chroma(state, pred_cu, lcu, &split_cu_loc[i]);
    }

    return sum + tr_tree_bits * state->lambda;
  }
  
  if (!skip_residual_coding) {
    cabac_data_t* cabac = (cabac_data_t*)&state->search_cabac;
    cabac_ctx_t* ctx = &(cabac->ctx.qt_cbf_model_cb[0]);
    cabac->cur_ctx = ctx;
    CABAC_FBITS_UPDATE(cabac, ctx, u_is_set, tr_tree_bits, "cbf_cb_search");
    
    ctx = &(cabac->ctx.qt_cbf_model_cr[u_is_set]);
    CABAC_FBITS_UPDATE(cabac, ctx, v_is_set, tr_tree_bits, "cbf_cb_search");
    
  }



  if (state->encoder_control->cfg.jccr) {
    int cbf_mask = u_is_set * 2 + v_is_set - 1;
    cabac_ctx_t* ctx = NULL;
    if (cbf_mask != -1) {
      cabac_data_t* cabac = (cabac_data_t*)&state->search_cabac;
      ctx = &(cabac->ctx.joint_cb_cr[cbf_mask]);
      CABAC_FBITS_UPDATE(cabac, ctx, 0, tr_tree_bits, "cbf_cb_search");
    }
  }

  // Chroma SSD
  int ssd = 0;
  if (!state->encoder_control->cfg.lossless) {
    int index = lcu_px.y * LCU_WIDTH_C + lcu_px.x;
    int ssd_u = uvg_pixels_calc_ssd(&lcu->ref.u[index], &lcu->rec.u[index],
                                    LCU_WIDTH_C,         LCU_WIDTH_C,
                                    cu_loc->chroma_width, cu_loc->chroma_height);
    int ssd_v = uvg_pixels_calc_ssd(&lcu->ref.v[index], &lcu->rec.v[index],
                                    LCU_WIDTH_C,        LCU_WIDTH_C,
                                    cu_loc->chroma_width, cu_loc->chroma_height);
    ssd = ssd_u + ssd_v;
  }

  if (!skip_residual_coding) {
    int8_t scan_order = uvg_get_scan_order(pred_cu->type, pred_cu->intra.mode_chroma, depth);

    // We need the rounded & shifted coordinates for the chroma coeff calculation
    cu_loc_t chroma_loc;
    uvg_cu_loc_ctor(&chroma_loc, lcu_px.x, lcu_px.y, cu_loc->width, cu_loc->height);

    if((pred_cu->joint_cb_cr & 3) == 0){
      coeff_bits += uvg_get_coeff_cost(state, lcu->coeff.u, NULL, &chroma_loc, 2, scan_order, 0, COEFF_ORDER_CU);
      coeff_bits += uvg_get_coeff_cost(state, lcu->coeff.v, NULL, &chroma_loc, 2, scan_order, 0, COEFF_ORDER_CU);
    }
    else {
      coeff_bits += uvg_get_coeff_cost(state, lcu->coeff.joint_uv, NULL, &chroma_loc, 2, scan_order, 0, COEFF_ORDER_CU);
      
    }
  }


  double bits = tr_tree_bits + coeff_bits;

  return (double)ssd * UVG_CHROMA_MULT + bits * state->c_lambda;
}

static double cu_rd_cost_tr_split_accurate(
  const encoder_state_t* const state,
  const cu_info_t* const pred_cu,
  lcu_t* const lcu,
  enum uvg_tree_type tree_type,
  uint8_t isp_cbf,
  const cu_loc_t* const cu_loc,
  const cu_loc_t* const chroma_loc,
  bool has_chroma) {
  const int width = cu_loc->width;
  const int height = cu_loc->height; // TODO: height for non-square blocks
  
  const int skip_residual_coding = pred_cu->skipped || (pred_cu->type != CU_INTRA && pred_cu->cbf == 0);
  // cur_cu is used for TU parameters.
  cu_info_t* const tr_cu = LCU_GET_CU_AT_PX(lcu, cu_loc->local_x, cu_loc->local_y);

  double coeff_bits = 0;
  double luma_bits = 0;
  double chroma_bits = 0;
  
  const int cb_flag_u = tr_cu->joint_cb_cr ? tr_cu->joint_cb_cr >> 1 : cbf_is_set(tr_cu->cbf, COLOR_U);
  const int cb_flag_v = tr_cu->joint_cb_cr ? tr_cu->joint_cb_cr & 1 : cbf_is_set(tr_cu->cbf, COLOR_V);

  cabac_data_t* cabac = (cabac_data_t*)&state->search_cabac;

  {
    int cbf = cbf_is_set_any(tr_cu->cbf);
    // Only need to signal coded block flag if not skipped or merged
    // skip = no coded residual, merge = coded residual
    if (pred_cu->type != CU_INTRA && (!pred_cu->merged)) {
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_qt_root_cbf_model), cbf, luma_bits, "rqt_root_cbf");
    }

  }
  
  if (cu_loc->width > TR_MAX_WIDTH || cu_loc->height > TR_MAX_WIDTH) {
    double sum = 0;
    enum split_type split;
    if(cu_loc->width > TR_MAX_WIDTH && cu_loc->height > TR_MAX_WIDTH) {
      split = QT_SPLIT;
    } else if(cu_loc->width > TR_MAX_WIDTH) {
      split = BT_VER_SPLIT;
    } else {
      split = BT_HOR_SPLIT;
    }
    
    cu_loc_t split_cu_loc[4];
    const int split_count= uvg_get_split_locs(cu_loc, split, split_cu_loc,NULL);
    cu_loc_t split_chroma_cu_loc[4];
    if (chroma_loc) {
      uvg_get_split_locs(chroma_loc, split, split_chroma_cu_loc, NULL);
    }
    for (int i = 0; i < split_count; ++i) {
      sum += cu_rd_cost_tr_split_accurate(state, pred_cu, lcu, tree_type, isp_cbf, &split_cu_loc[i], chroma_loc ? &split_chroma_cu_loc[i] : NULL, has_chroma);
    }
    return sum + luma_bits * state->lambda;
  }

  has_chroma = state->encoder_control->chroma_format != UVG_CSP_400 && has_chroma && tree_type != UVG_LUMA_T;
  if (!skip_residual_coding && has_chroma) {
    CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.qt_cbf_model_cb[0]), cb_flag_u, chroma_bits, "cbf_cb");  
    CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.qt_cbf_model_cr[cb_flag_u]), cb_flag_v, chroma_bits, "cbf_cr");    
  }

  const int cb_flag_y = cbf_is_set(tr_cu->cbf, COLOR_Y) && tree_type != UVG_CHROMA_T;

  const bool is_isp = !(pred_cu->type != CU_INTRA || pred_cu->intra.isp_mode == ISP_MODE_NO_ISP);
  // Add transform_tree cbf_luma bit cost.
  if (!is_isp) {
    const int is_tr_split = cu_loc->width > TR_MAX_WIDTH || cu_loc->height > TR_MAX_WIDTH;
    if ((pred_cu->type == CU_INTRA ||
      is_tr_split ||
      cb_flag_u ||
      cb_flag_v)
      && !skip_residual_coding && tree_type != UVG_CHROMA_T)
    {
      cabac_ctx_t* ctx = &(cabac->ctx.qt_cbf_model_luma[0]);

      CABAC_FBITS_UPDATE(cabac, ctx, cb_flag_y, luma_bits, "cbf_y_search");
    }
  }
  else {
    // TODO: 8x4 CUs
    const int split_limit = uvg_get_isp_split_num(width, height, pred_cu->intra.isp_mode, true);
    int luma_ctx = 2;
    const int split_limit_minus_one = split_limit - 1;
    for (int i = 0; i < split_limit; i++) {
      if (i != split_limit_minus_one || isp_cbf != 1 << split_limit_minus_one) {
        const int flag = (isp_cbf >> i) & 1;
        CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.qt_cbf_model_luma[luma_ctx]), flag, luma_bits, "cbf_y_search");
        luma_ctx = 2 + flag;
      }
    }
  }

  if (cb_flag_y || cb_flag_u || cb_flag_v) {
    // TODO qp_delta_sign_flag

    if ((cb_flag_u || cb_flag_v) && has_chroma && state->encoder_control->cfg.jccr) {
      CABAC_FBITS_UPDATE(cabac, &cabac->ctx.joint_cb_cr[cb_flag_u * 2 + cb_flag_v - 1], tr_cu->joint_cb_cr != 0, chroma_bits, "tu_joint_cbcr_residual_flag");
    }
  }


  // SSD between reconstruction and original
  unsigned luma_ssd = 0;
  if (!state->encoder_control->cfg.lossless && tree_type != UVG_CHROMA_T) {
    int index = cu_loc->local_x + LCU_WIDTH * cu_loc->local_y;
    luma_ssd = uvg_pixels_calc_ssd(&lcu->ref.y[index], &lcu->rec.y[index],
      LCU_WIDTH, LCU_WIDTH,
      width, height);
  }
  // Chroma transform skip enable/disable is non-normative, so we need to count the chroma
  // tr-skip bits even when we are never using it.
  const bool can_use_tr_skip = state->encoder_control->cfg.trskip_enable
                               && width <= (1 << state->encoder_control->cfg.trskip_max_size)
                               && height <= (1 << state->encoder_control->cfg.trskip_max_size)
                               && !is_isp;

  if(cb_flag_y || is_isp){
    if (can_use_tr_skip) {
      CABAC_FBITS_UPDATE(cabac, &cabac->ctx.transform_skip_model_luma, tr_cu->tr_idx == MTS_SKIP, luma_bits, "transform_skip_flag");
    }
    int8_t luma_scan_mode = SCAN_DIAG;
    if (pred_cu->type != CU_INTRA || pred_cu->intra.isp_mode == ISP_MODE_NO_ISP) {
      //const coeff_t* coeffs = &lcu->coeff.y[xy_to_zorder(LCU_WIDTH, x_px, y_px)];
      const coeff_t* coeffs = lcu->coeff.y;

      coeff_bits += uvg_get_coeff_cost(state, coeffs, tr_cu, cu_loc, 0, luma_scan_mode, pred_cu->tr_idx == MTS_SKIP, COEFF_ORDER_CU);
    }
    else {
      int split_type = pred_cu->intra.isp_mode;
      int split_limit = uvg_get_isp_split_num(width, height, split_type, true);

      for (int i = 0; i < split_limit; ++i) {
        cu_loc_t split_loc;
        uvg_get_isp_split_loc(&split_loc, cu_loc->x, cu_loc->y, width, height, i, split_type, true);
        const int part_x = split_loc.x;
        const int part_y = split_loc.y;

        // TODO: maybe just pass the cu_loc_t to these functions
        //const coeff_t* coeffs = &lcu->coeff.y[xy_to_zorder(LCU_WIDTH, part_x, part_y)];
        const coeff_t* coeffs = lcu->coeff.y;

        coeff_bits += uvg_get_coeff_cost(state, coeffs, tr_cu, &split_loc, 0, luma_scan_mode, pred_cu->tr_idx == MTS_SKIP, COEFF_ORDER_CU);
      }
    }
  }

  const bool is_local_sep_tree = (cu_loc->width != chroma_loc->width || cu_loc->height != chroma_loc->height) && state->encoder_control->chroma_format != UVG_CSP_400;

  if(is_local_sep_tree || tree_type == UVG_LUMA_T) {

    if (uvg_is_lfnst_allowed(state, tr_cu, is_local_sep_tree ? UVG_LUMA_T : tree_type, COLOR_Y, cu_loc, lcu)) {
      const int lfnst_idx = tr_cu->lfnst_idx;
      CABAC_FBITS_UPDATE(
        cabac,
        &cabac->ctx.lfnst_idx_model[1],
        lfnst_idx != 0,
        luma_bits,
        "lfnst_idx");
      if (lfnst_idx > 0) {
        CABAC_FBITS_UPDATE(
          cabac,
          &cabac->ctx.lfnst_idx_model[2],
          lfnst_idx == 2,
          luma_bits,
          "lfnst_idx");
      }
    }
    tr_cu->lfnst_last_scan_pos = false;
  }

  unsigned chroma_ssd = 0;
  if(has_chroma) {
    cu_loc_t temp_chroma_loc;
    const vector2d_t lcu_px = { chroma_loc->local_x >> 1, chroma_loc->local_y >> 1};
    uvg_cu_loc_ctor(&temp_chroma_loc, lcu_px.x, lcu_px.y, chroma_loc->width, chroma_loc->height);
    const int chroma_width  = chroma_loc->chroma_width;
    const int chroma_height = chroma_loc->chroma_height; 
    int8_t scan_order = SCAN_DIAG;
    //const unsigned index = xy_to_zorder(LCU_WIDTH_C, lcu_px.x, lcu_px.y);

    const bool chroma_can_use_tr_skip = state->encoder_control->cfg.trskip_enable
      && chroma_width <= (1 << state->encoder_control->cfg.trskip_max_size)
      && chroma_height <= (1 << state->encoder_control->cfg.trskip_max_size);
    if(pred_cu->joint_cb_cr == 0) {
      if (!state->encoder_control->cfg.lossless) {
        int index = lcu_px.y * LCU_WIDTH_C + lcu_px.x;
        unsigned ssd_u = uvg_pixels_calc_ssd(&lcu->ref.u[index], &lcu->rec.u[index],
          LCU_WIDTH_C, LCU_WIDTH_C, chroma_width, chroma_height) * state->chroma_weights[1];
        unsigned ssd_v = uvg_pixels_calc_ssd(&lcu->ref.v[index], &lcu->rec.v[index],
          LCU_WIDTH_C, LCU_WIDTH_C, chroma_width, chroma_height) * state->chroma_weights[2];
        chroma_ssd = ssd_u + ssd_v;
      }
      if(chroma_can_use_tr_skip && cb_flag_u) {
        CABAC_FBITS_UPDATE(cabac, &cabac->ctx.transform_skip_model_chroma, tr_cu->tr_skip & 2, chroma_bits, "transform_skip_flag");        
      }
      if(chroma_can_use_tr_skip && cb_flag_v) {
        CABAC_FBITS_UPDATE(cabac, &cabac->ctx.transform_skip_model_chroma, tr_cu->tr_skip & 4, chroma_bits, "transform_skip_flag");        
      }
      chroma_bits += uvg_get_coeff_cost(state, lcu->coeff.u, tr_cu, &temp_chroma_loc, COLOR_U, scan_order, tr_cu->tr_skip & 2, COEFF_ORDER_CU);
      chroma_bits += uvg_get_coeff_cost(state, lcu->coeff.v, tr_cu, &temp_chroma_loc, COLOR_V, scan_order, tr_cu->tr_skip & 4, COEFF_ORDER_CU);
      
    }
    else {
      {
        int index = lcu_px.y * LCU_WIDTH_C + lcu_px.x;
        int ssd_u_joint = uvg_pixels_calc_ssd(&lcu->ref.u[index], &lcu->rec.u[index],
          LCU_WIDTH_C, LCU_WIDTH_C, chroma_width, chroma_height) * state->chroma_weights[3];
        int ssd_v_joint = uvg_pixels_calc_ssd(&lcu->ref.v[index], &lcu->rec.v[index],
          LCU_WIDTH_C, LCU_WIDTH_C, chroma_width, chroma_height) * state->chroma_weights[3];
        chroma_ssd = ssd_u_joint + ssd_v_joint;
      }
      if (chroma_can_use_tr_skip) {
        CABAC_FBITS_UPDATE(cabac, &cabac->ctx.transform_skip_model_chroma, tr_cu->tr_skip & 2, chroma_bits, "transform_skip_flag");
      }
      chroma_bits += uvg_get_coeff_cost(state, lcu->coeff.joint_uv, tr_cu, &temp_chroma_loc, COLOR_U, scan_order, 0, COEFF_ORDER_CU);
    }
  }

  const bool is_chroma_tree = is_local_sep_tree || tree_type == UVG_CHROMA_T;
  if (uvg_is_lfnst_allowed(state, tr_cu, is_local_sep_tree ? UVG_CHROMA_T : tree_type, is_chroma_tree ? COLOR_UV : COLOR_Y, is_chroma_tree ? chroma_loc : cu_loc, lcu) && tree_type != UVG_LUMA_T) {
    const int lfnst_idx = is_chroma_tree ? tr_cu->cr_lfnst_idx : tr_cu->lfnst_idx;
    CABAC_FBITS_UPDATE(
      cabac,
      &cabac->ctx.lfnst_idx_model[is_chroma_tree],
      lfnst_idx != 0,
      luma_bits,
      "lfnst_idx");
    if (lfnst_idx > 0) {
      CABAC_FBITS_UPDATE(
        cabac,
        &cabac->ctx.lfnst_idx_model[2],
        lfnst_idx == 2,
        luma_bits,
        "lfnst_idx");
    }
  }
  tr_cu->lfnst_last_scan_pos = false;
  tr_cu->violates_lfnst_constrained_luma = false;
  tr_cu->violates_lfnst_constrained_chroma = false;
  if (uvg_is_mts_allowed(state, tr_cu, cu_loc) && tree_type != UVG_CHROMA_T) {

    bool symbol = tr_cu->tr_idx != 0;
    int ctx_idx = 0;
    CABAC_FBITS_UPDATE(cabac, &cabac->ctx.mts_idx_model[ctx_idx], symbol, luma_bits, "mts_idx");

    ctx_idx++;
    for (int i = 0; i < 3 && symbol; i++, ctx_idx++)
    {
      symbol = tr_cu->tr_idx > i + MTS_DST7_DST7 ? 1 : 0;
      CABAC_FBITS_UPDATE(cabac, &cabac->ctx.mts_idx_model[ctx_idx], symbol, luma_bits, "mts_idx");
    }
    tr_cu->mts_last_scan_pos = false;
    tr_cu->violates_mts_coeff_constraint = false;
  }

  double bits = luma_bits + coeff_bits;
  return luma_ssd * UVG_LUMA_MULT + chroma_ssd * UVG_CHROMA_MULT  + (bits + chroma_bits) * state->lambda;
}


// Return estimate of bits used to code prediction mode of cur_cu.
static double calc_mode_bits(
  const encoder_state_t *state,
  const lcu_t *lcu,
  const cu_info_t * cur_cu,
  const cu_loc_t* const cu_loc)
{
  assert(cur_cu->type == CU_INTRA);

  double mode_bits = uvg_luma_mode_bits(state, cur_cu, cu_loc, lcu);

  if (((cu_loc->width == 4 && cu_loc->x % 8 && cu_loc->y % 8) || (cu_loc->width != 4)) && state->encoder_control->chroma_format != UVG_CSP_400) {
    mode_bits += uvg_chroma_mode_bits(state, cur_cu->intra.mode_chroma, cur_cu->intra.mode);
  }

  return mode_bits;
}


// TODO: replace usages of this by the uvg_sort_indices_by_cost function.
/**
 * \brief Sort modes and costs to ascending order according to costs.
 */
void uvg_sort_modes(int8_t *__restrict modes, double *__restrict costs, uint8_t length)
{
  // Length for intra is always between 5 and 23, and is either 21, 17, 9 or 8 about
  // 60% of the time, so there should be no need for anything more complex
  // than insertion sort.
  // Length for merge is 5 or less.
  for (uint8_t i = 1; i < length; ++i) {
    const double cur_cost = costs[i];
    const int8_t cur_mode = modes[i];
    uint8_t j = i;
    while (j > 0 && cur_cost < costs[j - 1]) {
      costs[j] = costs[j - 1];
      modes[j] = modes[j - 1];
      --j;
    }
    costs[j] = cur_cost;
    modes[j] = cur_mode;
  }
}

/**
 * \brief Sort modes and costs to ascending order according to costs.
 */
void uvg_sort_modes_intra_luma(int8_t *__restrict modes, int8_t *__restrict trafo, double *__restrict costs, uint8_t length)
{
  // Length for intra is always between 5 and 23, and is either 21, 17, 9 or 8 about
  // 60% of the time, so there should be no need for anything more complex
  // than insertion sort.
  // Length for merge is 5 or less.
  for (uint8_t i = 1; i < length; ++i) {
    const double cur_cost = costs[i];
    const int8_t cur_mode = modes[i];
    const int8_t cur_tr = trafo[i];
    uint8_t j = i;
    while (j > 0 && cur_cost < costs[j - 1]) {
      costs[j] = costs[j - 1];
      modes[j] = modes[j - 1];
      trafo[j] = trafo[j - 1];
      --j;
    }
    costs[j] = cur_cost;
    modes[j] = cur_mode;
    trafo[j] = cur_tr;
  }
}

/**
 * \brief Sort keys (indices) to ascending order according to costs.
 */
void uvg_sort_keys_by_cost(unit_stats_map_t *__restrict map)
{
  // Size of sorted arrays is expected to be "small". No need for faster algorithm.
  for (uint8_t i = 1; i < map->size; ++i) {
    const int8_t cur_indx = map->keys[i];
    const double cur_cost = map->cost[cur_indx];
    uint8_t j = i;
    while (j > 0 && cur_cost < map->cost[map->keys[j - 1]]) {
      map->keys[j] = map->keys[j - 1];
      --j;
    }
    map->keys[j] = cur_indx;
  }
}


static void mark_deblocking(const cu_loc_t* const cu_loc, const cu_loc_t* const chroma_loc, lcu_t* lcu, enum uvg_tree_type tree_type, bool has_chroma, const bool is_separate_tree, int x_local, int y_local)
{
  if(tree_type != UVG_CHROMA_T) {
    if(cu_loc->x) {
      for (int x = cu_loc->local_x; x < cu_loc->local_x + cu_loc->width; x += TR_MAX_WIDTH) {
        for (int y = cu_loc->local_y; y < cu_loc->local_y + cu_loc->height; y += SCU_WIDTH) {
          LCU_GET_CU_AT_PX(lcu, x, y)->luma_deblocking |= EDGE_VER;
          if(!is_separate_tree && tree_type == UVG_BOTH_T) LCU_GET_CU_AT_PX(lcu, x, y)->chroma_deblocking |= EDGE_VER;
        }
      }
    }
    else if(cu_loc->width == 64) {
      for (int y = cu_loc->local_y; y < cu_loc->local_y + cu_loc->height; y += SCU_WIDTH) {
        LCU_GET_CU_AT_PX(lcu, TR_MAX_WIDTH, y)->luma_deblocking |= EDGE_VER;
        if (!is_separate_tree && tree_type == UVG_BOTH_T) LCU_GET_CU_AT_PX(lcu, TR_MAX_WIDTH, y)->chroma_deblocking |= EDGE_VER;
      }        
    }

    if(cu_loc->y) {
      for (int y = cu_loc->local_y; y < cu_loc->local_y + cu_loc->height; y += TR_MAX_WIDTH) {
        for (int x = cu_loc->local_x; x < cu_loc->local_x + cu_loc->width; x += SCU_WIDTH) {
          LCU_GET_CU_AT_PX(lcu, x, y)->luma_deblocking |= EDGE_HOR;
          if (!is_separate_tree && tree_type == UVG_BOTH_T) LCU_GET_CU_AT_PX(lcu, x, y)->chroma_deblocking |= EDGE_HOR;
        }
      }
    }
    else if (cu_loc->height == 64) {
      for (int x = cu_loc->local_x; x < cu_loc->local_x + cu_loc->width; x += SCU_WIDTH) {
        LCU_GET_CU_AT_PX(lcu, x, TR_MAX_WIDTH)->luma_deblocking |= EDGE_HOR;
        if (!is_separate_tree && tree_type == UVG_BOTH_T) LCU_GET_CU_AT_PX(lcu, x, TR_MAX_WIDTH)->chroma_deblocking |= EDGE_HOR;
      }
    }

    if(is_separate_tree && has_chroma) {
      if (chroma_loc->x) {
        for (int x = chroma_loc->local_x; x < chroma_loc->local_x + chroma_loc->width; x += TR_MAX_WIDTH) {
          for (int y = chroma_loc->local_y; y < chroma_loc->local_y + chroma_loc->height; y += SCU_WIDTH) {
            LCU_GET_CU_AT_PX(lcu, x, y)->chroma_deblocking |= EDGE_VER;
          }
        }
      }
      else if(cu_loc->width == 64) {
        for (int y = chroma_loc->local_y; y < chroma_loc->local_y + chroma_loc->height; y += SCU_WIDTH) {
          LCU_GET_CU_AT_PX(lcu, TR_MAX_WIDTH, y)->chroma_deblocking |= EDGE_VER;
        }          
      }

      if (chroma_loc->y) {
        for (int y = chroma_loc->local_y; y < chroma_loc->local_y + chroma_loc->height; y += TR_MAX_WIDTH) {
          for (int x = chroma_loc->local_x; x < chroma_loc->local_x + chroma_loc->width; x += SCU_WIDTH) {
            LCU_GET_CU_AT_PX(lcu, x, y)->chroma_deblocking |= EDGE_HOR;
          }
        }
      }
      else if (cu_loc->height == 64) {
        for (int x = chroma_loc->local_x; x < chroma_loc->local_x + chroma_loc->width; x += SCU_WIDTH) {
          LCU_GET_CU_AT_PX(lcu, x, TR_MAX_WIDTH)->chroma_deblocking |= EDGE_HOR;
        }
      }
    }
  }
  else {

    if (chroma_loc->x) {
      for (int x = x_local; x < x_local + chroma_loc->width; x += TR_MAX_WIDTH) {
        for (int y = y_local; y < y_local + chroma_loc->height; y += SCU_WIDTH) {
          LCU_GET_CU_AT_PX(lcu, x, y)->chroma_deblocking |= EDGE_VER;
        }
      }
    }
    else if(chroma_loc->width == 64) {
      for (int y = y_local; y < y_local + chroma_loc->height; y += SCU_WIDTH) {
        LCU_GET_CU_AT_PX(lcu, TR_MAX_WIDTH, y)->chroma_deblocking |= EDGE_VER;
      }        
    }

    if(chroma_loc->y) {
      for (int y = y_local; y < y_local + chroma_loc->height; y += TR_MAX_WIDTH) {
        for (int x = x_local; x < x_local + chroma_loc->width; x += SCU_WIDTH) {
          LCU_GET_CU_AT_PX(lcu, x, y)->chroma_deblocking |= EDGE_HOR;
        }
      }        
    }
    else if (chroma_loc->height == 64) {
      for (int x = x_local; x < x_local + chroma_loc->width; x += SCU_WIDTH) {
        LCU_GET_CU_AT_PX(lcu, x, TR_MAX_WIDTH)->chroma_deblocking |= EDGE_HOR;
      }
    }
  }
}

static bool check_for_early_termission(const int cu_width, const int cu_height, const cu_info_t* const cur_cu, int x_local, int y_local, const
                                       bool* improved,
                                       int cbf,
                                       lcu_t* split_lcu,
                                       int split_type,
                                       const bool* can_split)
{
  // Best no split has no residual and same direction bt didn't improve so don't try tt
  // 3.11
  if (
    !cbf && ((!improved[BT_VER_SPLIT] && split_type == TT_VER_SPLIT) ||
             (!improved[BT_HOR_SPLIT] && split_type == TT_HOR_SPLIT)))
    return true;


  // 3.8
  if (split_type == TT_HOR_SPLIT && can_split[BT_HOR_SPLIT]) {
    bool can_skip = true;
    for (int x_scu = x_local; x_scu < x_local + cu_width; x_scu += 4) {
      can_skip &=
        LCU_GET_CU_AT_PX(&split_lcu[BT_HOR_SPLIT - 1], x_scu, y_local)->log2_height == cur_cu->log2_height - 1 &&
        LCU_GET_CU_AT_PX(&split_lcu[BT_HOR_SPLIT - 1], x_scu, y_local + cu_height / 2)->log2_height == cur_cu->log2_height - 1;
    }
    if (can_skip) return true;
  }
  if (split_type == TT_VER_SPLIT && can_split[BT_VER_SPLIT]) {
    bool can_skip = true;
    for (int y_scu = y_local; y_scu < y_local + cu_height; y_scu += 4) {
      can_skip &=
        LCU_GET_CU_AT_PX(&split_lcu[BT_VER_SPLIT - 1], x_local, y_scu)->log2_width == cur_cu->log2_width - 1 &&
        LCU_GET_CU_AT_PX(&split_lcu[BT_VER_SPLIT - 1], x_local + cu_width / 2, y_scu)->log2_width == cur_cu->log2_width - 1;
    }
    if (can_skip) return true;
  }
  return false;
}

/**
 * Search every mode from 0 to MAX_PU_DEPTH and return cost of best mode.
 * - The recursion is started at depth 0 and goes in Z-order to MAX_PU_DEPTH.
 * - Data structure work_tree is maintained such that the neighbouring SCUs
 *   and pixels to the left and up of current CU are the final CUs decided
 *   via the search. This is done by copying the relevant data to all
 *   relevant levels whenever a decision is made whether to split or not.
 * - All the final data for the LCU gets eventually copied to depth 0, which
 *   will be the final output of the recursion.
 */
static double search_cu(
  encoder_state_t* const state,
  const cu_loc_t* const cu_loc,
  const cu_loc_t* const chroma_loc,
  lcu_t* lcu,
  enum uvg_tree_type tree_type,
  const split_tree_t split_tree,
  bool has_chroma,
  enum mode_type mode_type)
{
  const int depth = split_tree.current_depth;
  const encoder_control_t* ctrl = state->encoder_control;
  const videoframe_t * const frame = state->tile->frame;
  const int cu_width = cu_loc->width;
  const int cu_height =  cu_loc->height;
  const int x = cu_loc->x;
  const int y = cu_loc->y;
  const int luma_width = cu_loc->width;
  const int luma_height = cu_loc->height;
  
  const bool is_separate_tree = chroma_loc == NULL || cu_loc->height != chroma_loc->height || cu_loc->width != chroma_loc->width;
  assert(cu_width >= 4);
  double cost = MAX_DOUBLE;
  double inter_zero_coeff_cost = MAX_DOUBLE;
  double inter_bitcost = MAX_INT;
  cu_info_t *cur_cu;
  cabac_data_t pre_search_cabac;
  memcpy(&pre_search_cabac, &state->search_cabac, sizeof(pre_search_cabac));

  const uint32_t ctu_row = (cu_loc->y >> LOG2_LCU_WIDTH);
  const uint32_t ctu_row_mul_five = ctu_row * MAX_NUM_HMVP_CANDS;

  cu_info_t hmvp_lut[MAX_NUM_HMVP_CANDS];
  uint8_t hmvp_lut_size = state->tile->frame->hmvp_size[ctu_row];
  cu_info_t hmvp_lut_ibc[MAX_NUM_HMVP_CANDS];
  uint8_t hmvp_lut_size_ibc = state->tile->frame->hmvp_size_ibc[ctu_row];

  // Store original HMVP lut before search and restore after, since it's modified
  if (state->frame->slicetype != UVG_SLICE_I) memcpy(hmvp_lut, &state->tile->frame->hmvp_lut[ctu_row_mul_five], sizeof(cu_info_t) * MAX_NUM_HMVP_CANDS);
  if(state->encoder_control->cfg.ibc) memcpy(hmvp_lut_ibc, &state->tile->frame->hmvp_lut_ibc[ctu_row_mul_five], sizeof(cu_info_t) * MAX_NUM_HMVP_CANDS);

  struct {
    int32_t min;
    int32_t max;
  } pu_depth_inter, pu_depth_intra;
  
  int x_local = SUB_SCU(x);
  int y_local = SUB_SCU(y);

  int32_t frame_width = frame->width;
  int32_t frame_height = frame->height;
  // Stop recursion if the CU is completely outside the frame.
  if (x >= frame_width || y >= frame_height) {
    // Return zero cost because this CU does not have to be coded.
    return 0;
  }

  int gop_layer = ctrl->cfg.gop_len != 0 ? ctrl->cfg.gop[state->frame->gop_offset].layer - 1 : 0;

  // Assign correct depth limit
  constraint_t* constr = state->constraint;
  if(constr->ml_intra_depth_ctu) {
    pu_depth_intra.min = constr->ml_intra_depth_ctu->_mat_upper_depth[(x_local >> 3) + (y_local >> 3) * 8];
    pu_depth_intra.max = constr->ml_intra_depth_ctu->_mat_lower_depth[(x_local >> 3) + (y_local >> 3) * 8];
  }
  else {
    pu_depth_intra.min = ctrl->cfg.pu_depth_intra.min[gop_layer] >= 0 ? ctrl->cfg.pu_depth_intra.min[gop_layer] : ctrl->cfg.pu_depth_intra.min[0];
    pu_depth_intra.max = ctrl->cfg.pu_depth_intra.max[gop_layer] >= 0 ? ctrl->cfg.pu_depth_intra.max[gop_layer] : ctrl->cfg.pu_depth_intra.max[0];
  }

  pu_depth_inter.min = ctrl->cfg.pu_depth_inter.min[gop_layer] >= 0 ? ctrl->cfg.pu_depth_inter.min[gop_layer] : ctrl->cfg.pu_depth_inter.min[0];
  pu_depth_inter.max = ctrl->cfg.pu_depth_inter.max[gop_layer] >= 0 ? ctrl->cfg.pu_depth_inter.max[gop_layer] : ctrl->cfg.pu_depth_inter.max[0];

  cur_cu = LCU_GET_CU_AT_PX(lcu, x_local, y_local);
  memset(cur_cu, 0, sizeof(cu_info_t));
  // Assign correct depth
  cur_cu->type = CU_NOTSET;
  cur_cu->qp = state->qp;
  cur_cu->split_tree = split_tree.split_tree;
  cur_cu->mode_type_tree = split_tree.mode_type_tree | mode_type << (split_tree.current_depth * 2);
  cur_cu->log2_width = uvg_g_convert_to_log2[cu_width];
  cur_cu->log2_height = uvg_g_convert_to_log2[cu_height];

  if(chroma_loc) {
    cur_cu->log2_chroma_height = uvg_g_convert_to_log2[chroma_loc->chroma_height];
    cur_cu->log2_chroma_width = uvg_g_convert_to_log2[chroma_loc->chroma_width];
  }

  //Account for SCIPU. All sub-blocks of a SCIPU need to be either all intra or all inter
  bool is_scipu_inter = mode_type == MODE_TYPE_INTER;
  bool is_constrain_mode = mode_type != MODE_TYPE_ALL;


  intra_search_data_t intra_search = {0};

  const bool completely_inside = x + luma_width <= frame_width && y + luma_height <= frame_height;
  // If the CU is completely inside the frame at this depth, search for
  // prediction modes at this depth.
  if ( completely_inside)
  {
    int cu_width_inter_min = LCU_WIDTH >> pu_depth_inter.max;
    bool can_use_inter =
      state->frame->slicetype != UVG_SLICE_I &&
      split_tree.current_depth <= MAX_DEPTH &&
      (
        WITHIN(split_tree.current_depth, pu_depth_inter.min, pu_depth_inter.max) ||
        // When the split was forced because the CTU is partially outside the
        // frame, we permit inter coding even if pu_depth_inter would
        // otherwise forbid it.
        (x & ~(cu_width_inter_min - 1)) + cu_width_inter_min > frame_width ||
        (y & ~(cu_width_inter_min - 1)) + cu_width_inter_min > frame_height
        ) && (cu_loc->width != 4 || cu_loc->height != 4) && // 4x4 inter not allowed
      ((state->encoder_control->chroma_format != UVG_CSP_400) ? (cu_loc->chroma_height * cu_loc->chroma_width) >= 16 : true) && // Don't allow blocks with <16 chroma samples for now (TODO: investigate small chroma blocks)
      (is_constrain_mode ? is_scipu_inter : true); //Account for SCIPU

    if (can_use_inter) {
      double mode_cost;
      double mode_bitcost;
      uvg_search_cu_inter(state,
                          cu_loc, lcu,
                          &mode_cost,
                          &mode_bitcost);
      if (mode_cost < cost) {
        cost = mode_cost;
        inter_bitcost = mode_bitcost;
        cur_cu->type = CU_INTER;
      }
    }

    // Try to skip intra search in rd==0 mode.
    // This can be quite severe on bdrate. It might be better to do this
    // decision after reconstructing the inter frame.
    bool skip_intra = (state->encoder_control->cfg.rdo == 0
                      && cur_cu->type != CU_NOTSET
                      && cost / (cu_width * cu_width) < INTRA_THRESHOLD)
                      || (ctrl->cfg.early_skip && cur_cu->skipped);

    int32_t cu_width_intra_min = LCU_WIDTH >> pu_depth_intra.max;
    bool can_use_intra =
      (WITHIN(split_tree.current_depth, pu_depth_intra.min, pu_depth_intra.max) ||
        // When the split was forced because the CTU is partially outside
        // the frame, we permit intra coding even if pu_depth_intra would
        // otherwise forbid it.
        (x & ~(cu_width_intra_min - 1)) + cu_width_intra_min > frame_width ||
        (y & ~(cu_width_intra_min - 1)) + cu_width_intra_min > frame_height) &&
      !(state->encoder_control->cfg.force_inter && state->frame->slicetype != UVG_SLICE_I) &&
      (is_constrain_mode ? !is_scipu_inter : true); //Account for SCIPU;

    intra_search.cost = 0;
    if (can_use_intra && !skip_intra) {
      intra_search.pred_cu = *cur_cu;
      if(tree_type != UVG_CHROMA_T) {
        uvg_search_cu_intra(state, &intra_search, lcu, is_separate_tree ? UVG_LUMA_T : tree_type, cu_loc);
      }
#ifdef COMPLETE_PRED_MODE_BITS
      // Technically counting these bits would be correct, however counting
      // them universally degrades quality so this block is disabled by default
      if(state->frame->slicetype != UVG_SLICE_I) {
        double pred_mode_type_bits = 0;
        CABAC_FBITS_UPDATE(&state->search_cabac, &state->search_cabac.ctx.cu_pred_mode_model, 1, pred_mode_type_bits, "pred_mode_flag");
        CABAC_FBITS_UPDATE(&state->search_cabac, &state->search_cabac.ctx.cu_skip_flag_model[uvg_get_skip_context(x, y, lcu, NULL)], 0, pred_mode_type_bits, "skip_flag");
        intra_cost += pred_mode_type_bits * state->lambda;
      }
#endif
      if (state->encoder_control->cfg.cclm && tree_type != UVG_CHROMA_T && state->encoder_control->chroma_format != UVG_CSP_400) {
        if(intra_search.pred_cu.intra.isp_mode == ISP_MODE_NO_ISP) {
          uvg_intra_recon_cu(state,
                             &intra_search, cu_loc,
                             &intra_search.pred_cu, lcu,
                             tree_type,
                             true,
                             false);
        }
        else {
          cabac_data_t temp_cabac;
          memcpy(&temp_cabac, &state->search_cabac, sizeof(cabac_data_t));
          state->search_cabac.update = 1;
          uvg_recon_and_estimate_cost_isp(
            state,
            cu_loc,
            0,
            &intra_search,
            lcu,
            NULL
          );
          memcpy(&state->search_cabac, &temp_cabac, sizeof(cabac_data_t));
        }

        downsample_cclm_rec(
          state, x, y, cu_width / 2, cu_height / 2, lcu->rec.y, lcu->left_ref.y[64]
        );
      }
      double intra_cost = intra_search.cost;
      if (intra_cost < cost && tree_type != UVG_LUMA_T) {
        int8_t intra_mode = intra_search.pred_cu.intra.mode;
        
        if ((has_chroma || tree_type == UVG_CHROMA_T)
          && state->encoder_control->chroma_format != UVG_CSP_400) {

          intra_search.pred_cu.joint_cb_cr = 0;
          if(tree_type == UVG_CHROMA_T || is_separate_tree) {
            intra_mode = uvg_get_co_located_luma_mode(
                    chroma_loc, cu_loc, &intra_search.pred_cu, is_separate_tree ? lcu : NULL,
                    tree_type == UVG_CHROMA_T ? state->tile->frame->cu_array : NULL,
                    UVG_CHROMA_T);
            state->collocated_luma_mode = intra_mode;
            intra_search.pred_cu.type = CU_INTRA;
          } else  if (intra_search.pred_cu.intra.mip_flag) {
            intra_mode = 0;
          }
          intra_search.pred_cu.intra.mode_chroma = intra_mode;
          if (ctrl->cfg.rdo >= 2 || ctrl->cfg.jccr || ctrl->cfg.lfnst) {
            uvg_search_cu_intra_chroma(state, chroma_loc, lcu, &intra_search, intra_mode, tree_type, is_separate_tree);
          }
          else if (!intra_search.pred_cu.intra.mip_flag) {
            intra_search.pred_cu.intra.mode_chroma = intra_mode;
          }
          else {
            intra_search.pred_cu.intra.mode_chroma = 0;
          }
          state->quant_blocks[2].needs_init = true;
          uvg_intra_recon_cu(state,
                             &intra_search, chroma_loc,
                             &intra_search.pred_cu, lcu,
                             is_separate_tree ? UVG_CHROMA_T : tree_type,
                             false,
                             true);
          if(tree_type != UVG_CHROMA_T) {
            intra_cost += uvg_cu_rd_cost_chroma(state, &intra_search.pred_cu, lcu, chroma_loc);
          }
          else {
            intra_cost = intra_search.cost;
          }
          intra_search.pred_cu.violates_lfnst_constrained_chroma = false;
          intra_search.pred_cu.lfnst_last_scan_pos = false;
        }
        else {
          intra_search.pred_cu.intra.mode_chroma = intra_mode;
        }
      }
      if (intra_cost < cost) {
        cost = intra_cost;
        *cur_cu = intra_search.pred_cu;
        cur_cu->type = CU_INTRA;
      }
    }

    // Simple IBC search
    if (can_use_intra //&& state->frame->slicetype == UVG_SLICE_I
         && state->encoder_control->cfg.ibc 
         && cost > 1000
         && cu_width > 4
         && (x >= cu_width || y >= cu_width)
         && !cur_cu->skipped) {

      cu_info_t backup_cu = *cur_cu;

      double mode_cost;
      double mode_bitcost;
      uvg_search_cu_ibc(state,
                        cu_loc,
                        lcu,
                        &mode_cost, &mode_bitcost);
      if (mode_cost < cost) {
        cost = mode_cost;
        inter_bitcost = mode_bitcost;
        cur_cu->type = CU_IBC;
        cur_cu->inter.mv_dir = 1;
        cur_cu->joint_cb_cr = 0;
      } else {
        *cur_cu = backup_cu;
      }
    }

    // Reconstruct best mode because we need the reconstructed pixels for
    // mode search of adjacent CUs.
    if (cur_cu->type == CU_INTRA) {

      bool recon_chroma = true;
      bool recon_luma = tree_type != UVG_CHROMA_T && cur_cu->intra.isp_mode == ISP_MODE_NO_ISP;
      if (is_separate_tree || !has_chroma || state->encoder_control->chroma_format == UVG_CSP_400 || tree_type == UVG_LUMA_T || cu_loc->chroma_height % 4 == 2) {
        recon_chroma = false; 
      }
      lcu_fill_cu_info(lcu, x_local, y_local, cu_width, cu_height, cur_cu);
      if (!state->encoder_control->cfg.cclm && cur_cu->intra.isp_mode != ISP_MODE_NO_ISP) {
        uvg_recon_and_estimate_cost_isp(
          state,
          cu_loc,
          0,
          &intra_search,
          lcu,
          NULL
        );
      }
      else {
        uvg_intra_recon_cu(state,
          &intra_search, cu_loc,
          NULL, lcu,
          tree_type,
          recon_luma, recon_chroma);        
      }


      if((!recon_chroma && state->encoder_control->chroma_format != UVG_CSP_400 && tree_type != UVG_LUMA_T) 
        || tree_type == UVG_CHROMA_T) {
        intra_search.pred_cu.intra.mode_chroma = cur_cu->intra.mode_chroma;
        if(tree_type != UVG_CHROMA_T) {
          lcu_fill_chroma_cu_info(
            lcu,
            chroma_loc);
        }
        uvg_intra_recon_cu(state,
                           &intra_search, chroma_loc,
                           NULL, lcu,
                           UVG_CHROMA_T,
                           false,
                           true);
        lcu_fill_chroma_cbfs(
          lcu,
          chroma_loc,
          tree_type);
      } else {
        assert(cur_cu->cr_lfnst_idx == 0 && "If we don't have separate tree chroma lfnst index must be 0");
      }

      // Set isp split cbfs here
      const int split_type = intra_search.pred_cu.intra.isp_mode;
      const int split_num = split_type == ISP_MODE_NO_ISP || tree_type == UVG_CHROMA_T ? 0 : uvg_get_isp_split_num(cu_width, cu_height, split_type, true);

      const int cbf_cb = cbf_is_set(cur_cu->cbf, COLOR_U);
      const int cbf_cr = cbf_is_set(cur_cu->cbf, COLOR_V);
      const int jccr = cur_cu->joint_cb_cr;
      for (int i = 0; i < split_num; ++i) {
        cu_loc_t isp_loc;
        uvg_get_isp_split_loc(&isp_loc, x, y, cu_width, cu_height, i, split_type, true);
        // Fetching from CU array does not work for dimensions less than 4
        // Fetch proper x, y coords for isp blocks
        int tmp_x = isp_loc.x;
        int tmp_y = isp_loc.y;
        uvg_get_isp_cu_arr_coords(&tmp_x, &tmp_y, MAX(cu_width, cu_height));
        cu_info_t* split_cu = LCU_GET_CU_AT_PX(lcu, tmp_x % LCU_WIDTH, tmp_y % LCU_WIDTH);
        bool cur_cbf = (intra_search.best_isp_cbfs >> i) & 1;
        cbf_clear(&split_cu->cbf, COLOR_Y);
        cbf_clear(&split_cu->cbf, COLOR_U);
        cbf_clear(&split_cu->cbf, COLOR_V);
        if (cur_cbf) {
          cbf_set(&split_cu->cbf, COLOR_Y);
        }
        if(cbf_cb) cbf_set(&split_cu->cbf, COLOR_U);
        if(cbf_cr) cbf_set(&split_cu->cbf, COLOR_V);
        split_cu->joint_cb_cr = jccr;
      }
      lcu_fill_cu_info(lcu, x_local, y_local, cu_width, cu_height, cur_cu);


    } else if (cur_cu->type == CU_INTER || cur_cu->type == CU_IBC) {

      if (!cur_cu->skipped) {

        if (!cur_cu->merged) {
            if (cur_cu->inter.mv_dir & 1) uvg_round_precision(INTERNAL_MV_PREC, 2, &cur_cu->inter.mv[0][0], &cur_cu->inter.mv[0][1]);
            if (cur_cu->inter.mv_dir & 2) uvg_round_precision(INTERNAL_MV_PREC, 2, &cur_cu->inter.mv[1][0], &cur_cu->inter.mv[1][1]);
        }

        const bool has_chroma = state->encoder_control->chroma_format != UVG_CSP_400;
        uvg_inter_recon_cu(state, lcu, true, has_chroma, cu_loc);

        if (ctrl->cfg.zero_coeff_rdo && !ctrl->cfg.lossless && !ctrl->cfg.rdoq_enable && false) {
          //Calculate cost for zero coeffs
          // inter_zero_coeff_cost = cu_zero_coeff_cost(state, work_tree, cu_loc, split_tree.current_depth) + inter_bitcost * state->lambda;

        }
        cu_loc_t loc;
        uvg_cu_loc_ctor(&loc, x, y, cu_width, cu_height);
        uvg_quantize_lcu_residual(state,
                                  true, has_chroma && !cur_cu->joint_cb_cr,
                                  cur_cu->joint_cb_cr, &loc,
                                  NULL,
                                  lcu,
                                  false,
                                  tree_type);

        int cbf = cbf_is_set_any(cur_cu->cbf);

        if (cur_cu->merged && !cbf) {
          cur_cu->merged = 0;
          cur_cu->skipped = 1;
          // Selecting skip reduces bits needed to code the CU
          int skip_ctx = uvg_get_skip_context(x, y, lcu, NULL, NULL);
          inter_bitcost = CTX_ENTROPY_FBITS(&state->search_cabac.ctx.cu_skip_flag_model[skip_ctx], 1);
          inter_bitcost += CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.cu_merge_idx_ext_model), cur_cu->merge_idx != 0);
          inter_bitcost += cur_cu->merge_idx;        
        }
      }
      lcu_fill_cu_info(lcu, x_local, y_local, cu_width, cu_height, cur_cu);
      lcu_fill_cbf(lcu, x_local, y_local, cu_width, cu_height, cur_cu, UVG_BOTH_T);
    }
  }
  
  // The cabac functions assume chroma locations whereas the search uses luma locations
  // for the chroma tree, therefore we need to shift the chroma coordinates here for
  // passing to the bit cost calculating functions.
  cu_loc_t separate_tree_chroma_loc = *cu_loc;
  separate_tree_chroma_loc.y >>= 1;
  separate_tree_chroma_loc.x >>= 1;
  separate_tree_chroma_loc.width >>= 1;
  separate_tree_chroma_loc.height >>= 1;

  if (cur_cu->type == CU_INTRA || cur_cu->type == CU_INTER || cur_cu->type == CU_IBC) {
    double bits = 0;
    cabac_data_t* cabac  = &state->search_cabac;
    cabac->update = 1;
    
    bits += uvg_mock_encode_coding_unit(
      state,
      cabac,
      cu_loc,
      is_separate_tree && !has_chroma ? NULL : chroma_loc,
      lcu,
      cur_cu,
      tree_type,
      split_tree);

    
    cost = bits * state->lambda;

    cost += cu_rd_cost_tr_split_accurate(state, cur_cu, lcu, tree_type, intra_search.best_isp_cbfs, cu_loc, chroma_loc, has_chroma);
    //fprintf(stderr, "%4d %4d %2d %2d %d %d %f\n", x, y, cu_width, cu_height, has_chroma, cur_cu->split_tree, cost);
    
    //if (ctrl->cfg.zero_coeff_rdo && inter_zero_coeff_cost <= cost) {
    //  cost = inter_zero_coeff_cost;

    //  // Restore saved pixels from lower level of the working tree.
    //  copy_cu_pixels(&work_tree[split_tree.current_depth + 1], lcu, cu_loc, tree_type);

    //  if (cur_cu->merged) {
    //    cur_cu->merged = 0;
    //    cur_cu->skipped = 1;
    //    lcu_fill_cu_info(lcu, x_local, y_local, cu_width, cu_height, cur_cu);
    //  }

    //  cur_cu->cbf = 0;
    //  lcu_fill_cbf(lcu, x_local, y_local, cu_width, cur_cu);
    //}
    cabac->update = 0;

    mark_deblocking(
      cu_loc,
      chroma_loc,
      lcu,
      tree_type,
      has_chroma,
      is_separate_tree,
      x_local,
      y_local);
    if (cur_cu->type == CU_INTRA && cur_cu->intra.isp_mode != ISP_MODE_NO_ISP && tree_type != UVG_CHROMA_T) {
      const int split_num = uvg_get_isp_split_num( cu_width, cu_height, cur_cu->intra.isp_mode,true);
      for (int i = 1; i < split_num; i++) {
        cu_loc_t isp_loc;
        uvg_get_isp_split_loc(
          &isp_loc,
          x,
          y,
          cu_width,
          cu_height,
          i,
          cur_cu->intra.isp_mode,
          true);
        if (isp_loc.x % 4 || isp_loc.y % 4) continue;
        mark_deblocking(
          &isp_loc,
          chroma_loc,
          lcu,
          UVG_LUMA_T,
          false,
          false,
          isp_loc.local_x,
          isp_loc.local_y);
      }
    }
  } 

  bool can_split_cu =
    // If the CU is partially outside the frame, we need to split it even
    // if pu_depth_intra and pu_depth_inter would not permit it.
    cur_cu->type == CU_NOTSET ||
    (split_tree.current_depth < pu_depth_intra.max && !(state->encoder_control->cfg.force_inter&& state->frame->slicetype != UVG_SLICE_I)) ||
    (state->frame->slicetype != UVG_SLICE_I &&
      split_tree.current_depth < pu_depth_inter.max);

  if(state->encoder_control->cabac_debug_file) {
    fprintf(state->encoder_control->cabac_debug_file, "S %4d %4d %9d %d", x, y, split_tree.split_tree, tree_type);
    fwrite(&state->search_cabac.ctx, 1,  sizeof(state->search_cabac.ctx), state->encoder_control->cabac_debug_file);
  }

  bool can_split[6];
  bool is_implicit = uvg_get_possible_splits(state, cu_loc, split_tree, tree_type, mode_type, can_split);

  const int slice_type = state->frame->is_irap ? (tree_type == UVG_CHROMA_T ? 2 : 0) : 1;
  const int max_btd = state->encoder_control->cfg.max_btt_depth[slice_type];
  int minimum_split_amount;
  switch (slice_type) {
  case 0: minimum_split_amount = pu_depth_intra.min - split_tree.current_depth; break;
  case 1: minimum_split_amount = MIN(pu_depth_intra.min, pu_depth_inter.min) - split_tree.current_depth; break;
  case 2: minimum_split_amount = pu_depth_intra.min - split_tree.current_depth; break;
    default:
      assert(0 && "Incorrect_slice_type");
  }
  if(minimum_split_amount > max_btd && !is_implicit && can_split[1]) {
    // If search should not be performed at depths that cannot be reached after a maximum mtt split amount
    // we are in trouble, therefore prevent mtt splits in such situation
    can_split[2] = can_split[3] = can_split[4] = can_split[5] = false;
  }

  can_split_cu &= can_split[1] || can_split[2] || can_split[3] || can_split[4] || can_split[5];

  bool improved[6] = {false};

  // If skip mode was selected for the block, skip further search.
  // Skip mode means there's no coefficients in the block, so splitting
  // might not give any better results but takes more time to do.
  // It is ok to interrupt the search as soon as it is known that
  // the split costs at least as much as not splitting.
  int cbf = cbf_is_set_any(cur_cu->cbf);

  // 3.13
  if ((cu_height < 32 || cu_width < 32) && cur_cu->type != CU_NOTSET  && !cbf && split_tree.mtt_depth > 1 && tree_type != UVG_CHROMA_T) {
    can_split_cu = false;
  }

  if (can_split_cu && (cur_cu->type == CU_NOTSET || cbf || state->encoder_control->cfg.cu_split_termination == UVG_CU_SPLIT_TERMINATION_OFF || true)) {
    lcu_t * split_lcu = MALLOC(lcu_t, 5);
    enum split_type best_split = 0;
    double best_split_cost = MAX_DOUBLE;
    cabac_data_t post_seach_cabac;
    cabac_data_t best_split_cabac;
    memcpy(&post_seach_cabac, &state->search_cabac, sizeof(post_seach_cabac));
    
    cu_info_t best_split_hmvp_lut[MAX_NUM_HMVP_CANDS];
    uint8_t best_split_hmvp_lut_size = state->tile->frame->hmvp_size[ctu_row];
    cu_info_t best_split_hmvp_lut_ibc[MAX_NUM_HMVP_CANDS];
    uint8_t best_split_hmvp_lut_size_ibc = state->tile->frame->hmvp_size_ibc[ctu_row];

    // Recursively split all the way to max search depth.
    for (int split_type = QT_SPLIT; split_type <= TT_VER_SPLIT; ++split_type) {
      if (!can_split[split_type])
        continue;
      split_tree_t new_split = {
        split_tree.split_tree | split_type << (split_tree.current_depth * 3),
        split_tree.mode_type_tree, // | mode_type << (split_tree.current_depth * 2),
        split_tree.current_depth + 1,
        split_tree.mtt_depth + (split_type != QT_SPLIT),
        split_tree.implicit_mtt_depth + (split_type != QT_SPLIT && is_implicit),
        0,
        split_tree.scipu_cb_depth
      };

      if (completely_inside && check_for_early_termission(
            cu_width,
            cu_height,
            cur_cu,
            x_local,
            y_local,
            improved,
            cbf,
            split_lcu,
            split_type,
            can_split)) {
        can_split[split_type] = false;
        continue;
      }

      double split_cost = 0.0;
      memcpy(&state->search_cabac, &pre_search_cabac, sizeof(post_seach_cabac));


      double split_bits = 0;

      enum mode_type split_mode_type = mode_type;
      //TODO: use for searching both modes in scipu
      const enum mode_type_condition mode_type_cond = uvg_derive_mode_type_cond(cu_loc, state->frame->slicetype, tree_type, state->encoder_control->chroma_format, split_type, mode_type);
      if (mode_type_cond == MODE_TYPE_INFER) split_mode_type = MODE_TYPE_INTRA;
      //new_split.mode_type_tree = split_tree.mode_type_tree | split_mode_type << (split_tree.current_depth * 2);

      if (cur_cu->log2_height + cur_cu->log2_width > 4) {

        state->search_cabac.update = 1;
        // Add cost of cu_split_flag.
        const cu_info_t* left_cu = NULL, * above_cu = NULL;
        if (x) {
          if (x_local || tree_type != UVG_CHROMA_T) {
            left_cu = LCU_GET_CU_AT_PX(lcu, x_local - 1, y_local);
          }
          else {
            left_cu = uvg_cu_array_at_const(state->tile->frame->chroma_cu_array, x - 1, y);
          }
        }
        if (y) {
          if (y_local || tree_type != UVG_CHROMA_T) {
            above_cu = LCU_GET_CU_AT_PX(lcu, x_local, y_local - 1);
          }
          else {
            above_cu = uvg_cu_array_at_const(state->tile->frame->chroma_cu_array, x, y - 1);
          }
        }
        split_tree_t count_tree = split_tree;
        count_tree.split_tree = split_tree.split_tree | split_type << (split_tree.current_depth * 3);
        count_tree.mode_type_tree = split_tree.mode_type_tree; // | mode_type << (split_tree.current_depth * 2);
        const enum mode_type count_mode_type = mode_type_cond == MODE_TYPE_SIGNAL ? (cur_cu->type == CU_INTER ? MODE_TYPE_INTER : MODE_TYPE_INTRA) : split_mode_type;
        uvg_write_split_flag(
          state,
          &state->search_cabac,
          left_cu,
          above_cu, 
          cu_loc,
          count_tree,
          tree_type,
          count_mode_type,
          &is_implicit,
          &split_bits
          );
      }

      // 3.9
      const double factor    = state->qp > 30 ? 1.1 : 1.075;
      if (split_bits * state->lambda + cost / factor > cost) {
        can_split[split_type] = false;
        continue;
      }


      state->search_cabac.update = 0;
      split_cost += split_bits * state->lambda;

      // 3.7
      bool stop_to_qt = false;

      cu_loc_t new_cu_loc[4];
      uint8_t separate_chroma = 0;
      const int splits = uvg_get_split_locs(cu_loc, split_type, new_cu_loc, &separate_chroma);

      
      //TODO: remove if not necessary after inclusion of mode type 
      bool is_scipu = separate_chroma || split_tree.scipu_cb_depth > 0;
      if (is_scipu) new_split.scipu_cb_depth += 1;
      //Limit split when forced inter. Don't allow split if inter can't be used. TODO: Move to uvg_get_possible_splits?
      if (split_mode_type == MODE_TYPE_INTER &&
          ((new_cu_loc[0].width == 4 && new_cu_loc[0].height == 4) ||
           ((state->encoder_control->chroma_format != UVG_CSP_400) ? (new_cu_loc[0].chroma_height * new_cu_loc[0].chroma_width) < 16 : false))
         ) {
        can_split[split_type] = false;
        continue;
      }

      //Reset HMVP incase it has been modified while checking previous split types
      if (state->frame->slicetype != UVG_SLICE_I) {
        memcpy(&state->tile->frame->hmvp_lut[ctu_row_mul_five], hmvp_lut, sizeof(cu_info_t) * MAX_NUM_HMVP_CANDS);
        state->tile->frame->hmvp_size[ctu_row] = hmvp_lut_size;
      }
      if (state->encoder_control->cfg.ibc) {
        memcpy(&state->tile->frame->hmvp_lut_ibc[ctu_row_mul_five], hmvp_lut_ibc, sizeof(cu_info_t) * MAX_NUM_HMVP_CANDS);
        state->tile->frame->hmvp_size_ibc[ctu_row] = hmvp_lut_size_ibc;
      }

      separate_chroma |= !has_chroma;
      initialize_partial_work_tree(state, lcu, &split_lcu[split_type - 1], cu_loc , separate_chroma ? chroma_loc : cu_loc, tree_type);
      for (int split = 0; split < splits; ++split) {
        //if (split_mode_type == MODE_TYPE_ALL && is_scipu && split != 0) { //TODO: remove when proper search is added for scipu mode
        //  split_mode_type = (LCU_GET_CU_AT_PX( &split_lcu[split_type - 1],
        //                                           new_cu_loc[0].local_x,
        //                                           new_cu_loc[0].local_y
        //                                          )->type == CU_INTER) ? MODE_TYPE_INTER : MODE_TYPE_INTRA;
        //}
        new_split.mode_type_tree = new_split.mode_type_tree | split_mode_type << (split_tree.current_depth * 2);
        new_split.part_index = split;
        split_cost += search_cu(state, 
          &new_cu_loc[split], separate_chroma ? chroma_loc : &new_cu_loc[split],
          &split_lcu[split_type -1], 
          tree_type, new_split,
          !separate_chroma || (split == splits - 1 && has_chroma),
          split_mode_type);
        // If there is no separate chroma the block will always have chroma, otherwise it is the last block of the split that has the chroma

        // Set mode type for first split block
        if (split_mode_type == MODE_TYPE_ALL && is_scipu && split == 0) { //TODO: remove when proper search is added for scipu mode
          cu_info_t* const first_cu = LCU_GET_CU_AT_PX(&split_lcu[split_type - 1],
                                                        new_cu_loc[0].local_x,
                                                        new_cu_loc[0].local_y);
          split_mode_type = (first_cu->type == CU_INTER) ? MODE_TYPE_INTER : MODE_TYPE_INTRA;
          first_cu->mode_type_tree = first_cu->mode_type_tree | split_mode_type << (split_tree.current_depth * 2);
          first_cu->mode_type_tree = first_cu->mode_type_tree | split_mode_type << (new_split.current_depth * 2);
        }

        if (split_type == QT_SPLIT && completely_inside) {
          const cu_info_t * const t = LCU_GET_CU_AT_PX(
            &split_lcu[0],
            new_cu_loc[split].local_x,
            new_cu_loc[split].local_y);
          stop_to_qt |= GET_SPLITDATA(t, depth + 1) == QT_SPLIT;
        }

        if (split_cost > cost || split_cost > best_split_cost) {
          can_split[split_type] = false;
          break;
        }
      }

      //Check that sub-blocks for SCIPU are either all inter or all intra
      if (is_scipu /*&& can_split[split_type]*/) {
        const bool is_inter = LCU_GET_CU_AT_PX(&split_lcu[split_type - 1], new_cu_loc[0].local_x, new_cu_loc[0].local_y)->type == CU_INTER;
        for (int y_tmp = cu_loc->local_y; y_tmp < cu_loc->local_y + cu_loc->height; y_tmp += 4) {
          for (int x_tmp = cu_loc->local_x; x_tmp < cu_loc->local_x + cu_loc->width; x_tmp += 4) {
          //assert(is_inter ? LCU_GET_CU_AT_PX(&split_lcu[split_type - 1], new_cu_loc[split].local_x, new_cu_loc[split].local_y)->type == CU_INTER 
          //                : LCU_GET_CU_AT_PX(&split_lcu[split_type - 1], new_cu_loc[split].local_x, new_cu_loc[split].local_y)->type != CU_INTER);
            const cu_info_t * const t = LCU_GET_CU_AT_PX(&split_lcu[split_type - 1], x_tmp, y_tmp);
            if (t->type != CU_NOTSET) {
              assert(is_inter ? t->type == CU_INTER : t->type != CU_INTER);
            }
          }
        }
      }

      improved[split_type] = cost > split_cost;
      
      if (split_cost < best_split_cost) {
        best_split_cost = split_cost;
        best_split = split_type;
        memcpy(&best_split_cabac, &state->search_cabac, sizeof(cabac_data_t));

        //Store HMVP lut of best split
        if (state->frame->slicetype != UVG_SLICE_I) {
          memcpy(best_split_hmvp_lut, &state->tile->frame->hmvp_lut[ctu_row_mul_five], sizeof(cu_info_t) * MAX_NUM_HMVP_CANDS);
          best_split_hmvp_lut_size = state->tile->frame->hmvp_size[ctu_row];
        }
        if (state->encoder_control->cfg.ibc) {
          best_split_hmvp_lut_size_ibc = state->tile->frame->hmvp_size_ibc[ctu_row];
          memcpy(best_split_hmvp_lut_ibc, &state->tile->frame->hmvp_lut_ibc[ctu_row_mul_five], sizeof(cu_info_t) * MAX_NUM_HMVP_CANDS);
        }
      }
      if (stop_to_qt) break;
    }

    // If no search is not performed for this depth, try just the best mode
    // of the top left CU from the next depth. This should ensure that 64x64
    // gets used, at least in the most obvious cases, while avoiding any
    // searching.

    // TODO: Dual tree
    if (cur_cu->type == CU_NOTSET && depth < MAX_PU_DEPTH
        && x + cu_width <= frame_width && y + cu_width <= frame_height 
        && state->encoder_control->cfg.combine_intra_cus 
      && tree_type == UVG_BOTH_T)
    {

      cu_info_t *cu_d1 = LCU_GET_CU_AT_PX(&split_lcu[best_split - 1], x_local, y_local);

      // If the best CU in depth+1 is intra and the biggest it can be, try it.
      if (cu_d1->type == CU_INTRA && (cu_d1->log2_height + 1 == cur_cu->log2_height || cu_d1->log2_width + 1 == cur_cu->log2_width)) {
        cabac_data_t temp_cabac;
        memcpy(&temp_cabac, &state->search_cabac, sizeof(temp_cabac));
        memcpy(&state->search_cabac, &pre_search_cabac, sizeof(pre_search_cabac));
        cost = 0;
        double bits = 0;
        bool   is_implicit = false;
        uvg_write_split_flag(state, &state->search_cabac,
                             x > 0 ? LCU_GET_CU_AT_PX(lcu, SUB_SCU(x) - 1, SUB_SCU(y)) : NULL,
                             y > 0 ? LCU_GET_CU_AT_PX(lcu, SUB_SCU(x), SUB_SCU(y) - 1) : NULL, cu_loc, split_tree, tree_type, mode_type, &is_implicit,
                             &bits);

        cur_cu->intra = cu_d1->intra;
        cur_cu->type = CU_INTRA;
        if (cur_cu->intra.mode_chroma > 79) {
          cur_cu->intra.mode_chroma = cur_cu->intra.mode;
        }

        // Disable MRL in this case
        cur_cu->intra.multi_ref_idx = 0;
        cur_cu->lfnst_idx = 0;
        cur_cu->cr_lfnst_idx = 0;
        
        lcu_fill_cu_info(lcu, x_local, y_local, cu_width, cu_height, cur_cu);
        
        intra_search_data_t proxy;
        FILL(proxy, 0);
        proxy.pred_cu = *cur_cu;

        uvg_intra_recon_cu(state,
                           &proxy, cu_loc,
                           NULL,
                           lcu,
                           tree_type,
                           true,
                           state->encoder_control->chroma_format != UVG_CSP_400);

        double mode_bits = calc_mode_bits(state, lcu, cur_cu, cu_loc) + bits;
        cost += mode_bits * state->lambda;

        cost += cu_rd_cost_tr_split_accurate(state, cur_cu, lcu, tree_type, 0, cu_loc, chroma_loc, has_chroma);

        mark_deblocking(cu_loc, chroma_loc, lcu, tree_type, has_chroma, is_separate_tree, x_local, y_local);

        memcpy(&post_seach_cabac, &state->search_cabac, sizeof(post_seach_cabac));
        memcpy(&state->search_cabac, &temp_cabac, sizeof(temp_cabac));
      }
    }

    if (best_split_cost < cost) {
      // Copy split modes to this depth.
      cost = best_split_cost;
      memcpy(&state->search_cabac, &best_split_cabac, sizeof(best_split_cabac));
      work_tree_copy_up(&split_lcu[best_split -1], lcu, state->encoder_control->cfg.jccr, tree_type, cu_loc, is_separate_tree && !has_chroma ? NULL : chroma_loc);
      downsample_cclm_rec(
        state, x, y, cu_width / 2, cu_height / 2, lcu->rec.y, lcu->left_ref.y[64]
      );

      //Need to restore best split HMVP
      if (state->frame->slicetype != UVG_SLICE_I) {
        memcpy(&state->tile->frame->hmvp_lut[ctu_row_mul_five], best_split_hmvp_lut, sizeof(cu_info_t) * MAX_NUM_HMVP_CANDS);
        state->tile->frame->hmvp_size[ctu_row] = best_split_hmvp_lut_size;
      }
      if (state->encoder_control->cfg.ibc) {
        memcpy(&state->tile->frame->hmvp_lut_ibc[ctu_row_mul_five], best_split_hmvp_lut_ibc, sizeof(cu_info_t) * MAX_NUM_HMVP_CANDS);
        state->tile->frame->hmvp_size_ibc[ctu_row] = best_split_hmvp_lut_size_ibc;
      }
#if UVG_DEBUG
      //debug_split = 1;
#endif
    } else if (depth > 0) {
      // Copy this CU's mode all the way down for use in adjacent CUs mode
      // search.
      memcpy(&state->search_cabac, &post_seach_cabac, sizeof(post_seach_cabac));
      downsample_cclm_rec(
        state, x, y, cu_width / 2, cu_height / 2, lcu->rec.y, lcu->left_ref.y[64]
      );

      if (state->frame->slicetype != UVG_SLICE_I) {
        // Reset HMVP to the beginning of this CU level search and add this CU as the mvp
        memcpy(&state->tile->frame->hmvp_lut[ctu_row_mul_five], hmvp_lut, sizeof(cu_info_t) * MAX_NUM_HMVP_CANDS);
        state->tile->frame->hmvp_size[ctu_row] = hmvp_lut_size;        
      }
      if (state->encoder_control->cfg.ibc) {        
        memcpy(&state->tile->frame->hmvp_lut_ibc[ctu_row_mul_five], hmvp_lut_ibc, sizeof(cu_info_t) * MAX_NUM_HMVP_CANDS);
        state->tile->frame->hmvp_size_ibc[ctu_row] = hmvp_lut_size_ibc;        
      }
      // Add candidate when in inter slice or ibc is enabled
      if(state->frame->slicetype != UVG_SLICE_I || state->encoder_control->cfg.ibc) {
        uvg_hmvp_add_mv(state, x, y, cu_width, cu_height, cur_cu);
      }
    }
    else {
      downsample_cclm_rec(
        state, x, y, cu_width / 2, cu_height / 2, lcu->rec.y, lcu->left_ref.y[64]
      );      
    }
    FREE_POINTER(split_lcu);
  } else if (cur_cu->log2_height + cur_cu->log2_width > 4) {
    // Need to copy modes down since the lower level of the work tree is used
    // when searching SMP and AMP blocks.
    if(tree_type != UVG_CHROMA_T) {
      downsample_cclm_rec(
        state, x, y, cu_width / 2, cu_height / 2, lcu->rec.y, lcu->left_ref.y[64]
      );
    }

    if (state->frame->slicetype != UVG_SLICE_I) {
      // Reset HMVP to the beginning of this CU level search and add this CU as the mvp
      memcpy(&state->tile->frame->hmvp_lut[ctu_row_mul_five], hmvp_lut, sizeof(cu_info_t) * MAX_NUM_HMVP_CANDS);
      state->tile->frame->hmvp_size[ctu_row] = hmvp_lut_size;        
    }
    if (state->encoder_control->cfg.ibc) {        
      memcpy(&state->tile->frame->hmvp_lut_ibc[ctu_row_mul_five], hmvp_lut_ibc, sizeof(cu_info_t) * MAX_NUM_HMVP_CANDS);
      state->tile->frame->hmvp_size_ibc[ctu_row] = hmvp_lut_size_ibc;        
    }
    // Add candidate when in inter slice or ibc is enabled
    if(state->frame->slicetype != UVG_SLICE_I || state->encoder_control->cfg.ibc) {
      uvg_hmvp_add_mv(state, x, y, cu_width, cu_height, cur_cu);
    }
  }

  assert(cur_cu->type != CU_NOTSET);

  return cost;
}


/**
 * Initialize lcu_t for search.
 * - Copy reference CUs.
 * - Copy reference pixels from neighbouring LCUs.
 * - Copy reference pixels from this LCU.
 */
static void init_lcu_t(const encoder_state_t * const state, const int x, const int y, lcu_t *lcu, const yuv_t *hor_buf, const yuv_t *ver_buf)
{
  const videoframe_t * const frame = state->tile->frame;

  FILL(*lcu, 0);
  
  lcu->rec.chroma_format = state->encoder_control->chroma_format;
  lcu->ref.chroma_format = state->encoder_control->chroma_format;

  // Copy reference cu_info structs from neighbouring LCUs.

  // Copy top CU row.
  if (y > 0) {
    for (int i = 0; i < LCU_WIDTH; i += SCU_WIDTH) {
      const cu_info_t *from_cu = uvg_cu_array_at_const(frame->cu_array, x + i, y - 1);
      cu_info_t *to_cu = LCU_GET_CU_AT_PX(lcu, i, -1);
      memcpy(to_cu, from_cu, sizeof(*to_cu));
    }
  }
  // Copy left CU column.
  if (x > 0) {
    for (int i = 0; i < LCU_WIDTH; i += SCU_WIDTH) {
      const cu_info_t *from_cu = uvg_cu_array_at_const(frame->cu_array, x - 1, y + i);
      cu_info_t *to_cu = LCU_GET_CU_AT_PX(lcu, -1, i);
      memcpy(to_cu, from_cu, sizeof(*to_cu));
    }
  }
  // Copy top-left CU.
  if (x > 0 && y > 0) {
    const cu_info_t *from_cu = uvg_cu_array_at_const(frame->cu_array, x - 1, y - 1);
    cu_info_t *to_cu = LCU_GET_CU_AT_PX(lcu, -1, -1);
    memcpy(to_cu, from_cu, sizeof(*to_cu));
  }

  // Copy top-right CU, available only without WPP
  if (y > 0 && x + LCU_WIDTH < frame->width && !state->encoder_control->cfg.wpp) {
    const cu_info_t *from_cu = uvg_cu_array_at_const(frame->cu_array, x + LCU_WIDTH, y - 1);
    cu_info_t *to_cu = LCU_GET_TOP_RIGHT_CU(lcu);
    memcpy(to_cu, from_cu, sizeof(*to_cu));
  }

  // Copy reference pixels.
  {
    const int pic_width = frame->width;
    // Copy top reference pixels.
    if (y > 0) {
      // hor_buf is of size pic_width so there might not be LCU_REF_PX_WIDTH
      // number of allocated pixels left.
      int x_max = MIN(LCU_REF_PX_WIDTH, pic_width - x);
      int x_min_in_lcu = (x>0) ? 0 : 1;
      int luma_offset = OFFSET_HOR_BUF(x, y, frame, x_min_in_lcu - 1);
      int chroma_offset = OFFSET_HOR_BUF_C(x, y, frame, x_min_in_lcu - 1);
      int luma_bytes = (x_max + (1 - x_min_in_lcu))*sizeof(uvg_pixel);
      int chroma_bytes = (x_max / 2 + (1 - x_min_in_lcu))*sizeof(uvg_pixel);

      memcpy(&lcu->top_ref.y[x_min_in_lcu], &hor_buf->y[luma_offset], luma_bytes);

      if (state->encoder_control->chroma_format != UVG_CSP_400) {
        memcpy(&lcu->top_ref.u[x_min_in_lcu], &hor_buf->u[chroma_offset], chroma_bytes);
        memcpy(&lcu->top_ref.v[x_min_in_lcu], &hor_buf->v[chroma_offset], chroma_bytes);
      }
    }
    // Copy left reference pixels.
    if (x > 0) {
      int y_min_in_lcu = (y>0) ? 0 : 1;
      int luma_offset = OFFSET_VER_BUF(x, y, frame, y_min_in_lcu - 1);
      int chroma_offset = OFFSET_VER_BUF_C(x, y, frame, y_min_in_lcu - 1);
      int luma_bytes = (LCU_WIDTH + (1 - y_min_in_lcu)) * sizeof(uvg_pixel);
      int chroma_bytes = (LCU_WIDTH / 2 + (1 - y_min_in_lcu)) * sizeof(uvg_pixel);

      memcpy(&lcu->left_ref.y[y_min_in_lcu], &ver_buf->y[luma_offset], luma_bytes);

      if (state->encoder_control->chroma_format != UVG_CSP_400) {
        memcpy(&lcu->left_ref.u[y_min_in_lcu], &ver_buf->u[chroma_offset], chroma_bytes);
        memcpy(&lcu->left_ref.v[y_min_in_lcu], &ver_buf->v[chroma_offset], chroma_bytes);
      }
    }
  }

  // Copy LCU pixels.
  {
    const videoframe_t * const frame = state->tile->frame;
    int x_max = MIN(x + LCU_WIDTH, frame->width) - x;
    int y_max = MIN(y + LCU_WIDTH, frame->height) - y;

    int x_c = x / 2;
    int y_c = y / 2;
    int x_max_c = x_max / 2;
    int y_max_c = y_max / 2;

    uvg_pixel* source = NULL;
    if (state->tile->frame->lmcs_aps->m_sliceReshapeInfo.sliceReshaperEnableFlag) {
      source = frame->source_lmcs->y;
    } else {
      source = frame->source->y;
    }

    // Use LMCS pixels for luma if they are available, otherwise source_lmcs is mapped to normal source
    uvg_pixels_blit(&source[x + y * frame->source->stride], lcu->ref.y,
                        x_max, y_max, frame->source->stride, LCU_WIDTH);
    if (state->encoder_control->chroma_format != UVG_CSP_400) {
      uvg_pixels_blit(&frame->source->u[x_c + y_c * frame->source->stride / 2], lcu->ref.u,
                      x_max_c, y_max_c, frame->source->stride / 2, LCU_WIDTH / 2);
      uvg_pixels_blit(&frame->source->v[x_c + y_c * frame->source->stride / 2], lcu->ref.v,
                      x_max_c, y_max_c, frame->source->stride / 2, LCU_WIDTH / 2);
    }
  }
}


/**
 * Copy CU and pixel data to it's place in picture datastructure.
 */
static void copy_lcu_to_cu_data(const encoder_state_t * const state, int x_px, int y_px, const lcu_t *lcu, enum
                                uvg_tree_type tree_type)
{
  // Copy non-reference CUs to picture.
  uvg_cu_array_copy_from_lcu(
    tree_type != UVG_CHROMA_T ? state->tile->frame->cu_array : state->tile->frame->chroma_cu_array, 
    x_px,
    y_px,
    lcu);

  // Copy pixels to picture.
  {
    videoframe_t * const pic = state->tile->frame;
    const int pic_width = pic->width;
    const int x_max = MIN(x_px + LCU_WIDTH, pic_width) - x_px;
    const int y_max = MIN(y_px + LCU_WIDTH, pic->height) - y_px;

    if(tree_type != UVG_CHROMA_T) {
      uvg_pixels_blit(lcu->rec.y, &pic->rec->y[x_px + y_px * pic->rec->stride],
                          x_max, y_max, LCU_WIDTH, pic->rec->stride);
    }

    if (state->tile->frame->lmcs_aps->m_sliceReshapeInfo.sliceReshaperEnableFlag) {
      uvg_pixels_blit(lcu->rec.y, &pic->rec_lmcs->y[x_px + y_px * pic->rec->stride],
        x_max, y_max, LCU_WIDTH, pic->rec->stride);
    }

    if (state->encoder_control->chroma_format != UVG_CSP_400 && tree_type != UVG_LUMA_T) {
      uvg_pixels_blit(lcu->rec.u, &pic->rec->u[(x_px / 2) + (y_px / 2) * (pic->rec->stride / 2)],
                      x_max / 2, y_max / 2, LCU_WIDTH / 2, pic->rec->stride / 2);
      uvg_pixels_blit(lcu->rec.v, &pic->rec->v[(x_px / 2) + (y_px / 2) * (pic->rec->stride / 2)],
                      x_max / 2, y_max / 2, LCU_WIDTH / 2, pic->rec->stride / 2);
    }
  }
}


/**
 * Search LCU for modes.
 * - Best mode gets copied to current picture.
 */
void uvg_search_lcu(encoder_state_t * const state, const int x, const int y, const yuv_t * const hor_buf, const yuv_t * const ver_buf, lcu_coeff_t *coeff)
{
  memcpy(&state->search_cabac, &state->cabac, sizeof(cabac_data_t));
  state->search_cabac.only_count = 1;
  assert(x % LCU_WIDTH == 0);
  assert(y % LCU_WIDTH == 0);

  // Initialize the same starting state to every depth. The search process
  // will use these as temporary storage for predictions before making
  // a decision on which to use, and they get updated during the search
  // process.
  lcu_t work_tree;
  init_lcu_t(state, x, y, &work_tree, hor_buf, ver_buf);

  // If the ML depth prediction is enabled, 
  // generate the depth prediction interval 
  // for the current lcu
  constraint_t* constr = state->constraint;
  if (constr->ml_intra_depth_ctu) {
    uvg_lcu_luma_depth_pred(constr->ml_intra_depth_ctu, work_tree.ref.y, state->qp);
  }

  int tree_type = state->frame->slicetype == UVG_SLICE_I
                  && state->encoder_control->cfg.dual_tree
                    ? UVG_LUMA_T
                    : UVG_BOTH_T;

  cu_loc_t start;
  uvg_cu_loc_ctor(&start, x, y, LCU_WIDTH, LCU_WIDTH);
  split_tree_t split_tree = { 0, 0, 0, 0, 0, 0, 0 };
  double cost = 0; //TODO: save cost when resuming

#ifdef UVG_ENCODING_RESUME
  if (uvg_can_resume_encoding(state, x, y, false)) {
    uvg_process_resume_encoding(state, x, y, false, &cost, &work_tree, true);
  } else {
#endif // UVG_ENCODING_RESUME

  // Start search from depth 0.
  cost = search_cu(
    state,
    &start,
    &start,
    &work_tree,
    tree_type,
    split_tree,
    tree_type == UVG_BOTH_T,
    MODE_TYPE_ALL);

#ifdef UVG_ENCODING_RESUME
    uvg_process_resume_encoding(state, x, y, false, &cost, &work_tree, false);
  }
#endif // UVG_ENCODING_RESUME

  // Save squared cost for rate control.
  if(state->encoder_control->cfg.rc_algorithm == UVG_LAMBDA) {
    uvg_get_lcu_stats(state, x / LCU_WIDTH, y / LCU_WIDTH)->weight = cost * cost;
  }

  // The best decisions through out the LCU got propagated back to depth 0,
  // so copy those back to the frame.
  copy_lcu_to_cu_data(state, x, y, &work_tree, tree_type);

  // Copy coeffs to encoder state.
  copy_coeffs(work_tree.coeff.y, coeff->y, LCU_WIDTH, LCU_WIDTH, LCU_WIDTH);

  if(state->frame->slicetype == UVG_SLICE_I && state->encoder_control->cfg.dual_tree) {
#ifdef UVG_ENCODING_RESUME
    if (uvg_can_resume_encoding(state, x, y, true)) {
      uvg_process_resume_encoding(state, x, y, true, &cost, &work_tree, true);
    } else {
#endif // UVG_ENCODING_RESUME

    cost = search_cu(
      state, &start,
      &start,
      &work_tree, UVG_CHROMA_T,
      split_tree,
      true,
      MODE_TYPE_ALL);

#ifdef UVG_ENCODING_RESUME
       uvg_process_resume_encoding(state, x, y, true, &cost, &work_tree, false);
    }
#endif // UVG_ENCODING_RESUME

    if (state->encoder_control->cfg.rc_algorithm == UVG_LAMBDA) {
      uvg_get_lcu_stats(state, x / LCU_WIDTH, y / LCU_WIDTH)->weight += cost * cost;
    }
    copy_lcu_to_cu_data(state, x, y, &work_tree, UVG_CHROMA_T);
  }

  copy_coeffs(work_tree.coeff.u, coeff->u, LCU_WIDTH_C, LCU_WIDTH_C, LCU_WIDTH_C);
  copy_coeffs(work_tree.coeff.v, coeff->v, LCU_WIDTH_C, LCU_WIDTH_C, LCU_WIDTH_C);
  if (state->encoder_control->cfg.jccr) {
    copy_coeffs(work_tree.coeff.joint_uv, coeff->joint_uv, LCU_WIDTH_C, LCU_WIDTH_C, LCU_WIDTH_C);
  }
}
