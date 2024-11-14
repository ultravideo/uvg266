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

#include <string.h>
#include <stdlib.h>

#include "cu.h"

#include "alf.h"
#include "encoderstate.h"
#include "threads.h"



/**
 * \brief Number of PUs in a CU.
 *
 * Indexed by part_mode_t values.
 */
const uint8_t uvg_part_mode_num_parts[] = {
  1, // 2Nx2N
  2, // 2NxN
  2, // Nx2N
  4, // NxN
  2, // 2NxnU
  2, // 2NxnD
  2, // nLx2N
  2, // nRx2N
};

/**
 * \brief PU offsets.
 *
 * Indexed by [part mode][PU number][axis].
 *
 * Units are 1/4 of the width of the CU.
 */
const uint8_t uvg_part_mode_offsets[][4][2] = {
  { {0, 0}                         }, // 2Nx2N
  { {0, 0}, {0, 2}                 }, // 2NxN
  { {0, 0}, {2, 0}                 }, // Nx2N
  { {0, 0}, {2, 0}, {0, 2}, {2, 2} }, // NxN
  { {0, 0}, {0, 1}                 }, // 2NxnU
  { {0, 0}, {0, 3}                 }, // 2NxnD
  { {0, 0}, {1, 0}                 }, // nLx2N
  { {0, 0}, {3, 0}                 }, // nRx2N
};

/**
 * \brief PU sizes.
 *
 * Indexed by [part mode][PU number][axis].
 *
 * Units are 1/4 of the width of the CU.
 */
const uint8_t uvg_part_mode_sizes[][4][2] = {
  { {4, 4}                         }, // 2Nx2N
  { {4, 2}, {4, 2}                 }, // 2NxN
  { {2, 4}, {2, 4}                 }, // Nx2N
  { {2, 2}, {2, 2}, {2, 2}, {2, 2} }, // NxN
  { {4, 1}, {4, 3}                 }, // 2NxnU
  { {4, 3}, {4, 1}                 }, // 2NxnD
  { {1, 4}, {3, 4}                 }, // nLx2N
  { {3, 4}, {1, 4}                 }, // nRx2N
};


cu_info_t* uvg_cu_array_at(cu_array_t *cua, unsigned x_px, unsigned y_px)
{
  return (cu_info_t*) uvg_cu_array_at_const(cua, x_px, y_px);
}


void uvg_get_isp_cu_arr_coords(int *x, int *y, int dim)
{
  // Do nothing if dimensions are divisible by 4
  if (*y % 4 == 0 && *x % 4 == 0) return;
  const int remainder_y = *y % 4;
  const int remainder_x = *x % 4;

  if (remainder_y != 0) {
    // Horizontal ISP split
    if (remainder_y % 2 == 0 && dim == 8) {
      // 8x2 block
      *y -= 2;
      *x += 4;
    }
    else {
      // 16x1 block
      *y -= remainder_y;
      *x += remainder_y * 4;
    }
  }
  else {
    // Vertical ISP split
    if (*x % 2 == 0 && dim == 8) {
      // 2x8 block
      *y += 4;
      *x -= 2;
    }
    else {
      // 1x16 block
      *y += remainder_x * 4;
      *x -= remainder_x;
    }
  }
}


const cu_info_t* uvg_cu_array_at_const(const cu_array_t *cua, unsigned x_px, unsigned y_px)
{
  assert(x_px < cua->width);
  assert(y_px < cua->height);
  return &(cua)->data[(x_px >> 2) + (y_px >> 2) * ((cua)->stride >> 2)];
}


/**
 * \brief Allocate a CU array.
 *
 * \param width   width of the array in luma pixels
 * \param height  height of the array in luma pixels
 */
cu_array_t * uvg_cu_array_alloc(const int width, const int height)
{
  cu_array_t *cua = MALLOC(cu_array_t, 1);
  if (cua == NULL) return NULL;

  // Round up to a multiple of LCU width and divide by cell width.
  const int width_scu  = CEILDIV(width,  LCU_WIDTH) * LCU_WIDTH / SCU_WIDTH;
  const int height_scu = CEILDIV(height, LCU_WIDTH) * LCU_WIDTH / SCU_WIDTH;
  const unsigned cu_array_size = width_scu * height_scu;

  cua->base     = NULL;
  cua->data     = calloc(cu_array_size, sizeof(cu_info_t));
  cua->width    = width_scu  * SCU_WIDTH;
  cua->height   = height_scu * SCU_WIDTH;
  cua->stride   = cua->width;
  cua->refcount = 1;

  return cua;
}
cu_array_t * uvg_cu_array_chroma_alloc(const int width, const int height, enum uvg_chroma_format chroma)
{
  cu_array_t *cua = MALLOC(cu_array_t, 1);
  if (cua == NULL) return NULL;

  // Round up to a multiple of LCU width and divide by cell width.
  const int chroma_height = chroma == UVG_CSP_444 ? LCU_WIDTH : LCU_WIDTH_C;
  const int width_scu  = CEILDIV(width,  LCU_WIDTH_C) * LCU_WIDTH_C / SCU_WIDTH;
  const int height_scu = CEILDIV(height, chroma_height) * chroma_height / SCU_WIDTH;
  const unsigned cu_array_size = width_scu * height_scu;

  cua->base     = NULL;
  cua->data     = calloc(cu_array_size, sizeof(cu_info_t));
  cua->width    = width_scu  * SCU_WIDTH;
  cua->height   = height_scu * SCU_WIDTH;
  cua->stride   = cua->width;
  cua->refcount = 1;

  return cua;
}


cu_array_t * uvg_cu_subarray(cu_array_t *base,
                             const unsigned x_offset,
                             const unsigned y_offset,
                             const unsigned width,
                             const unsigned height)
{
  assert(x_offset + width <= base->width);
  assert(y_offset + height <= base->height);

  if (x_offset == 0 &&
      y_offset == 0 &&
      width == base->width &&
      height == base->height)
  {
    return uvg_cu_array_copy_ref(base);
  }

  cu_array_t *cua = MALLOC(cu_array_t, 1);
  if (cua == NULL) return NULL;

  // Find the real base array.
  cu_array_t *real_base = base;
  while (real_base->base) {
    real_base = real_base->base;
  }
  cua->base     = uvg_cu_array_copy_ref(real_base);
  cua->data     = uvg_cu_array_at(base, x_offset, y_offset);
  cua->width    = width;
  cua->height   = height;
  cua->stride   = base->stride;
  cua->refcount = 1;

  return cua;
}

void uvg_cu_array_free(cu_array_t **cua_ptr)
{
  cu_array_t *cua = *cua_ptr;
  if (cua == NULL) return;
  *cua_ptr = NULL;

  int new_refcount = UVG_ATOMIC_DEC(&cua->refcount);
  if (new_refcount > 0) {
    // Still we have some references, do nothing.
    return;
  }

  assert(new_refcount == 0);

  if (!cua->base) {
    FREE_POINTER(cua->data);
  } else {
    uvg_cu_array_free(&cua->base);
    cua->data = NULL;
  }

  FREE_POINTER(cua);
}


/**
 * \brief Get a new pointer to a cu array.
 *
 * Increment reference count and return the cu array.
 */
cu_array_t * uvg_cu_array_copy_ref(cu_array_t* cua)
{
  int32_t new_refcount = UVG_ATOMIC_INC(&cua->refcount);
  // The caller should have had another reference and we added one
  // reference so refcount should be at least 2.
  assert(new_refcount >= 2);
  return cua;
}


/**
 * \brief Copy an lcu to a cu array.
 *
 * All values are in luma pixels.
 *
 * \param dst     destination array
 * \param dst_x   x-coordinate of the left edge of the copied area in dst
 * \param dst_y   y-coordinate of the top edge of the copied area in dst
 * \param src     source lcu
 */
void uvg_cu_array_copy_from_lcu(cu_array_t* dst, int dst_x, int dst_y, const lcu_t *src)
{
  const int dst_stride = dst->stride >> 2;
  const int width = LCU_WIDTH;
  for (int y = 0; y < width; y += SCU_WIDTH) {
    for (int x = 0; x < width; x += SCU_WIDTH) {
      const cu_info_t *from_cu = LCU_GET_CU_AT_PX(src, x, y);
      const int x_scu = (dst_x + x) >> 2;
      const int y_scu = (dst_y + y) >> 2;
      cu_info_t *to_cu = &dst->data[x_scu + y_scu * dst_stride];
      memcpy(to_cu,                  from_cu, sizeof(*to_cu));
    }
  }
}

/*
 * \brief Constructs cu_loc_t based on given parameters. Calculates chroma dimensions automatically.
 *
 * \param loc     Destination cu_loc.
 * \param x       Block top left x coordinate.
 * \param y       Block top left y coordinate.
 * \param width   Block width.
 * \param height  Block height.
*/
void uvg_cu_loc_ctor(cu_loc_t* loc, int x, int y, int width, int height)
{
  assert(x >= 0 && y >= 0 && width >= 0 && height >= 0 && "Cannot give negative coordinates or block dimensions.");
  assert(!(width > LCU_WIDTH || height > LCU_WIDTH) && "Luma CU dimension exceeds maximum (dim > LCU_WIDTH).");
  // This check is no longer valid. With non-square blocks and ISP enabled, even 1x16 and 16x1 (ISP needs at least 16 samples) blocks are valid
  //assert(!(width < 4 || height < 4) && "Luma CU dimension smaller than 4.");
  
  loc->x = x;
  loc->y = y;
  loc->local_x = x % LCU_WIDTH;
  loc->local_y = y % LCU_WIDTH;
  loc->width = width;
  loc->height = height;
  // TODO: when MTT is implemented, chroma dimensions can be minimum 2.
  // Chroma width is half of luma width, when not at maximum depth.
  loc->chroma_width = width >> 1;
  loc->chroma_height = height >> 1;
}


int uvg_get_split_locs(
  const cu_loc_t* const origin,
  enum split_type split,
  cu_loc_t out[4],
  uint8_t* separate_chroma)
{
  const int half_width = origin->width >> 1;
  const int half_height = origin->height >> 1;
  const int quarter_width = origin->width >> 2;
  const int quarter_height = origin->height >> 2;
  if (origin->width == 4 && separate_chroma) *separate_chroma = 1;

  switch (split) {
    case NO_SPLIT:
      assert(0 && "trying to get split from no split");
    break;
    case QT_SPLIT:
      uvg_cu_loc_ctor(&out[0], origin->x, origin->y, half_width, half_height);
      uvg_cu_loc_ctor(&out[1], origin->x + half_width, origin->y, half_width, half_height);
      uvg_cu_loc_ctor(&out[2], origin->x, origin->y + half_height, half_width, half_height);
      uvg_cu_loc_ctor(&out[3], origin->x + half_width, origin->y + half_height, half_width, half_height);
      if (half_height == 4 && separate_chroma) *separate_chroma = 1;
      return 4;
    case BT_HOR_SPLIT:
      uvg_cu_loc_ctor(&out[0], origin->x, origin->y, origin->width, half_height);
      uvg_cu_loc_ctor(&out[1], origin->x, origin->y + half_height, origin->width, half_height);
      if (half_height * origin->width < 64 && separate_chroma) *separate_chroma = 1;
      return 2;
    case BT_VER_SPLIT:
      uvg_cu_loc_ctor(&out[0], origin->x, origin->y, half_width, origin->height);
      uvg_cu_loc_ctor(&out[1], origin->x + half_width, origin->y, half_width, origin->height);
      if ((half_width == 4 || half_width * origin->height < 64) && separate_chroma) *separate_chroma = 1;
      return 2;
    case TT_HOR_SPLIT:
      uvg_cu_loc_ctor(&out[0], origin->x, origin->y, origin->width, quarter_height);
      uvg_cu_loc_ctor(&out[1], origin->x, origin->y + quarter_height, origin->width, half_height);
      uvg_cu_loc_ctor(&out[2], origin->x, origin->y + quarter_height + half_height, origin->width, quarter_height);
      if (quarter_height * origin->width < 64 && separate_chroma) *separate_chroma = 1;
      return 3;
    case TT_VER_SPLIT:
      uvg_cu_loc_ctor(&out[0], origin->x, origin->y, quarter_width, origin->height);
      uvg_cu_loc_ctor(&out[1], origin->x + quarter_width, origin->y, half_width, origin->height);
      uvg_cu_loc_ctor(&out[2], origin->x + quarter_width + half_width, origin->y, quarter_width, origin->height);
      if ((quarter_width == 4 || quarter_width * origin->height < 64) && separate_chroma) *separate_chroma = 1;
      return 3;
  }
  return 0;
}


int uvg_get_implicit_split(
  const encoder_state_t* const state,
  const cu_loc_t* const cu_loc,
  uint8_t max_mtt_depth)
{
  bool right_ok = (state->tile->frame->width) >= cu_loc->x + cu_loc->width;
  bool bottom_ok = (state->tile->frame->height) >= cu_loc->y + cu_loc->height;

  if (right_ok && bottom_ok) return NO_SPLIT;
  if (right_ok && max_mtt_depth != 0) return BT_HOR_SPLIT;
  if (bottom_ok && max_mtt_depth != 0) return BT_VER_SPLIT;
  return QT_SPLIT;
}


enum mode_type_condition uvg_derive_mode_type_cond(const cu_loc_t* const cu_loc, const enum uvg_slice_type slice_type, const enum uvg_tree_type tree_type,
                                                   const enum uvg_chroma_format chroma_format, const enum split_type split_type, const enum mode_type mode_type)
{
  const bool is_dual_tree = slice_type == UVG_SLICE_I && tree_type != UVG_BOTH_T;
  if (is_dual_tree || mode_type != MODE_TYPE_ALL || chroma_format == UVG_CSP_444 || chroma_format == UVG_CSP_400) {
    return MODE_TYPE_INHERIT;
  }

  const unsigned luma_area = cu_loc->width * cu_loc->height;
  if ((luma_area == 64 && (split_type == QT_SPLIT || split_type == TT_HOR_SPLIT || split_type == TT_VER_SPLIT)) ||
      (luma_area == 32 && (split_type == BT_HOR_SPLIT || split_type == BT_VER_SPLIT))) {
    return MODE_TYPE_INFER;
  }

  if ((luma_area == 64 && (split_type == BT_HOR_SPLIT || split_type == BT_VER_SPLIT) && chroma_format == UVG_CSP_420) ||
      (luma_area == 128 && (split_type == TT_HOR_SPLIT || split_type == TT_VER_SPLIT) && chroma_format == UVG_CSP_420) ||
      (cu_loc->width == 8 && split_type == BT_VER_SPLIT) ||
      (cu_loc->width == 16 && split_type == TT_VER_SPLIT)) {
    return slice_type != UVG_SLICE_I ? MODE_TYPE_SIGNAL : MODE_TYPE_INFER;
  }

  return MODE_TYPE_INHERIT;
}

int uvg_get_possible_splits(const encoder_state_t * const state,
                            const cu_loc_t * const cu_loc, split_tree_t split_tree, enum uvg_tree_type tree_type,
                            bool splits[6])
{
  const unsigned width = cu_loc->width;
  const unsigned height = cu_loc->height;
  const int slice_type = state->frame->is_irap ? (tree_type == UVG_CHROMA_T ? 2 : 0) : 1;

  enum mode_type mode_type_parent = GET_MODETYPEDATA(&split_tree, split_tree.current_depth-1);

  const unsigned max_btd =
    state->encoder_control->cfg.max_btt_depth[slice_type] + split_tree.implicit_mtt_depth;
  const unsigned max_bt_size = state->encoder_control->cfg.max_bt_size[slice_type];
  const unsigned min_bt_size = 1 << MIN_SIZE;
  const unsigned max_tt_size = state->encoder_control->cfg.max_tt_size[slice_type];
  const unsigned min_tt_size = 1 << MIN_SIZE;
  const unsigned min_qt_size = state->encoder_control->cfg.min_qt_size[slice_type];

  const enum split_type implicitSplit = uvg_get_implicit_split(state, cu_loc, max_btd);
  
  splits[NO_SPLIT] = splits[QT_SPLIT] = splits[BT_HOR_SPLIT] = splits[TT_HOR_SPLIT] = splits[BT_VER_SPLIT] = splits[TT_VER_SPLIT] = true;
  bool can_btt = split_tree.mtt_depth < max_btd;
  
  const enum split_type last_split = GET_SPLITDATA(&split_tree, split_tree.current_depth - 1);
  const enum split_type parl_split = last_split == TT_HOR_SPLIT ? BT_HOR_SPLIT : BT_VER_SPLIT;

  if (tree_type == UVG_CHROMA_T && mode_type_parent == MODE_TYPE_INTRA) {
    splits[QT_SPLIT] = splits[BT_HOR_SPLIT] = splits[TT_HOR_SPLIT] = splits[BT_VER_SPLIT] = splits[TT_VER_SPLIT] = false;
    return implicitSplit != NO_SPLIT;
  }

  // don't allow QT-splitting below a BT split
  if (split_tree.current_depth != 0 && last_split != QT_SPLIT /* && !(width > 64 || height > 64)*/) splits[QT_SPLIT] = false;
  if (width <= min_qt_size)                              splits[QT_SPLIT] = false;

  if (tree_type == UVG_CHROMA_T && width <= 8) splits[QT_SPLIT] = false;

  if (implicitSplit != NO_SPLIT)
  {
    splits[NO_SPLIT] = splits[TT_HOR_SPLIT] = splits[TT_VER_SPLIT] = false;

    splits[BT_HOR_SPLIT] = implicitSplit == BT_HOR_SPLIT && height <= max_bt_size;
    splits[BT_VER_SPLIT] = implicitSplit == BT_VER_SPLIT && width <= max_bt_size;
    if (tree_type == UVG_CHROMA_T && width <= 8) splits[BT_VER_SPLIT] = false;
    if (!splits[BT_HOR_SPLIT] && !splits[BT_VER_SPLIT] && !splits[QT_SPLIT]) splits[QT_SPLIT] = true;
    return 1;
  }

  if ((last_split == TT_HOR_SPLIT || last_split == TT_VER_SPLIT) && split_tree.part_index == 1)
  {
    splits[BT_HOR_SPLIT] = parl_split != BT_HOR_SPLIT;
    splits[BT_VER_SPLIT] = parl_split != BT_VER_SPLIT;
  }

  if (can_btt && (width <= min_bt_size && height <= min_bt_size)
    && ((width <= min_tt_size && height <= min_tt_size)))
  {
    can_btt = false;
  }
  if (can_btt && (width > max_bt_size || height > max_bt_size)
    && ((width > max_tt_size || height > max_tt_size)))
  {
    can_btt = false;
  }

  if (!can_btt)
  {
    splits[BT_HOR_SPLIT] = splits[TT_HOR_SPLIT] = splits[BT_VER_SPLIT] = splits[TT_VER_SPLIT] = false;

    return 0;
  }

  if (width > max_bt_size || height > max_bt_size)
  {
    splits[BT_HOR_SPLIT] = splits[BT_VER_SPLIT] = false;
  }

  // specific check for BT splits
  if (height <= min_bt_size)                            splits[BT_HOR_SPLIT] = false;
  if (width > 64 && height <= 64) splits[BT_HOR_SPLIT] = false;
  if (tree_type == UVG_CHROMA_T && width * height <= 64)     splits[BT_HOR_SPLIT] = false;

  if (width <= min_bt_size)                              splits[BT_VER_SPLIT] = false;
  if (width <= 64 && height > 64) splits[BT_VER_SPLIT] = false;
  if (tree_type == UVG_CHROMA_T && (width * height <= 64 || width <= 8))     splits[BT_VER_SPLIT] = false;

  if (mode_type_parent == MODE_TYPE_INTER && width * height == 32)  splits[BT_VER_SPLIT] = splits[BT_HOR_SPLIT] = false;

  if (height <= 2 * min_tt_size || height > max_tt_size || width > max_tt_size)
    splits[TT_HOR_SPLIT] = false;
  if (width > 64 || height > 64)  splits[TT_HOR_SPLIT] = false;
  if (tree_type == UVG_CHROMA_T && width * height <= 64 * 2)     splits[TT_HOR_SPLIT] = false;

  if (width <= 2 * min_tt_size || width > max_tt_size || height > max_tt_size)
    splits[TT_VER_SPLIT] = false;
  if (width > 64 || height > 64)  splits[TT_VER_SPLIT] = false;
  if (tree_type == UVG_CHROMA_T && (width * height <= 64 * 2 || width <= 16))     splits[TT_VER_SPLIT] = false;

  if (mode_type_parent == MODE_TYPE_INTER && width * height == 64)  splits[TT_VER_SPLIT] = splits[TT_HOR_SPLIT] = false;

  return 0;
}


int uvg_count_available_edge_cus(const cu_loc_t* const cu_loc, const lcu_t* const lcu, bool left)
{
  if ((left && cu_loc->x == 0) || (!left && cu_loc->y == 0)) {
    return 0;
  }
  if (left && cu_loc->local_x == 0) return (LCU_WIDTH - cu_loc->local_y) / 4;
  if (!left && cu_loc->local_y == 0) return (cu_loc->width) / 2;

  int amount = left ? cu_loc->height & ~3 : cu_loc->width & ~3;
  if(left) {
    const cu_info_t* cu = LCU_GET_CU_AT_PX(lcu, cu_loc->local_x, cu_loc->local_y);
    if (cu_loc->local_y == 0 && cu_loc->local_x == 32 && cu->log2_height == 6 && cu->log2_width == 6) return 8;
    while (cu_loc->local_y + amount < LCU_WIDTH && LCU_GET_CU_AT_PX(lcu, cu_loc->local_x - TR_MIN_WIDTH, cu_loc->local_y + amount)->type != CU_NOTSET) {
      amount += TR_MIN_WIDTH;
    }
    return MAX(amount / TR_MIN_WIDTH, cu_loc->height / TR_MIN_WIDTH);
  }
  while (cu_loc->local_x + amount < LCU_WIDTH && LCU_GET_CU_AT_PX(lcu, cu_loc->local_x + amount, cu_loc->local_y - TR_MIN_WIDTH)->type != CU_NOTSET) {
    amount += TR_MIN_WIDTH;
  }
  return MAX(amount / TR_MIN_WIDTH, cu_loc->width / TR_MIN_WIDTH);
}
