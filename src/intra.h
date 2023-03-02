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
* \ingroup Reconstruction
* \file
* Intra prediction.
*/

#include "cu.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "uvg266.h"

// Maximum possible reference line length for intra blocks
#define INTRA_REF_LENGTH (2 * 128 + 3 + 33 * MAX_REF_LINE_IDX)

typedef struct {
  uvg_pixel left[INTRA_REF_LENGTH];
  uvg_pixel top[INTRA_REF_LENGTH];
} uvg_intra_ref;
typedef struct
{
  uvg_intra_ref ref;
  uvg_intra_ref filtered_ref;
  bool filtered_initialized;
} uvg_intra_references;

typedef struct
{
  int16_t a;
  int16_t shift;
  int16_t b;
} cclm_parameters_t;

typedef struct {
  cu_info_t pred_cu;
  cclm_parameters_t cclm_parameters[2];
  double cost;
  double bits;
  double coeff_bits;
  double distortion;
  double lfnst_costs[3];
  uint8_t best_isp_cbfs;
} intra_search_data_t ;


#define UVG_NUM_INTRA_MODES 67

/**
* \brief Function for deriving intra luma predictions
* \param x          x-coordinate of the PU in pixels
* \param y          y-coordinate of the PU in pixels
* \param preds      output buffer for 3 predictions
* \param cur_pu     PU to check
* \param left_pu    PU to the left of cur_pu
* \param above_pu   PU above cur_pu
* \returns          1 if predictions are found, otherwise 0
*/
int8_t uvg_intra_get_dir_luma_predictor(
  const uint32_t x,
  const uint32_t y,
  int8_t *preds,
  const cu_info_t *const cur_pu,
  const cu_info_t *const left_pu,
  const cu_info_t *const above_pu);

/**
* \brief Build intra prediction reference buffers.
* \param log2_width    Log2 of width, range 2..5.
* \param color         What color pixels to use.
* \param luma_px       Luma coordinates of the prediction block.
* \param pic_px        Picture dimensions in luma pixels.
* \param lcu           LCU struct.
* \param refs          Pointer to top and left references.
* \param entropy_sync  Indicate that top right is not available if WPP is enabled.
* \param extra_refs    Additional left edge reference lines for use with MRL.
* \param multi_ref_idx Multi reference line index for the prediction block.
*/
void uvg_intra_build_reference(
  const encoder_state_t* const state,
  const cu_loc_t* const pu_loc,
  const cu_loc_t* const cu_loc,
  const color_t color,
  const vector2d_t *const luma_px,
  const vector2d_t *const pic_px,
  const lcu_t *const lcu,
  uvg_intra_references *const refs,
  bool entropy_sync,
  uvg_pixel *extra_refs,
  uint8_t multi_ref_idx,
  const uint8_t isp_mode);

/**
 * \brief Generate intra predictions.
 * \param refs            Reference pixels used for the prediction.
 * \param log2_width      Width of the predicted block.
 * \param mode            Intra mode used for the prediction.
 * \param color           Color of the prediction.
 * \param dst             Buffer for the predicted pixels.
 * \param filter_boundary Whether to filter the boundary on modes 10 and 26.
 */
void uvg_intra_predict(
  const encoder_state_t* const state,
  uvg_intra_references* const refs,
  const cu_loc_t* const cu_loc,
  const cu_loc_t* const pu_loc,
  const color_t color,
  uvg_pixel* dst,
  const intra_search_data_t* data,
  const lcu_t* lcu
);

void uvg_intra_recon_cu(
  encoder_state_t* const state,
  intra_search_data_t* search_data,
  const cu_loc_t* cu_loc,
  cu_info_t *cur_cu,
  lcu_t *lcu,
  enum uvg_tree_type tree_type,
  bool recon_luma,
  bool recon_chroma);

double uvg_recon_and_estimate_cost_isp(encoder_state_t* const state,
                                       const cu_loc_t* const cu_loc,
                                       double cost_treshold,
                                       intra_search_data_t* const search_data,
                                       lcu_t* const lcu, bool* violates_lfnst);

int8_t uvg_get_co_located_luma_mode(
  const cu_loc_t* const chroma_loc,
  const cu_loc_t* const cu_loc,
  const cu_info_t* luma_cu,
  const lcu_t* const lcu,
  const cu_array_t* const cu_array,
  enum uvg_tree_type tree_type);
bool uvg_cclm_is_allowed(const encoder_state_t* const state, const cu_loc_t* const luma_loc, cu_info_t const* const cur_cu, enum
                         uvg_tree_type tree_type);

uint8_t uvg_get_mip_flag_context(
  const cu_loc_t* const cu_loc,
  const lcu_t* lcu,
  cu_array_t* const cu_a);

int8_t uvg_wide_angle_correction(
  int_fast8_t mode,
  const int log2_width,
  const int log2_height,
  const bool account_for_dc_planar);

// ISP related defines
#define NUM_ISP_MODES 3
#define ISP_MODE_NO_ISP 0
#define ISP_MODE_HOR 1
#define ISP_MODE_VER 2
#define SPLIT_TYPE_HOR 1
#define SPLIT_TYPE_VER 2

int uvg_get_isp_split_dim(const int width, const int height, const int split_type, const bool is_transform_block);
int uvg_get_isp_split_num(const int width, const int height, const int split_type, const bool is_transform_block);
void uvg_get_isp_split_loc(cu_loc_t *loc, const int x, const int y, const int block_w, const int block_h, int split_idx, const int split_type, const bool is_transform_block);
bool uvg_can_use_isp(const int width, const int height);
bool uvg_can_use_isp_with_lfnst(const int width, const int height, const int isp_mode, const enum uvg_tree_type tree_type);
