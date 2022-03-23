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
#include "kvazaar.h"

// Maximum possible reference line length for intra blocks
#define INTRA_REF_LENGTH (2 * 128 + 3 + 33 * MAX_REF_LINE_IDX)

typedef struct {
  kvz_pixel left[INTRA_REF_LENGTH];
  kvz_pixel top[INTRA_REF_LENGTH];
} kvz_intra_ref;
typedef struct
{
  kvz_intra_ref ref;
  kvz_intra_ref filtered_ref;
  bool filtered_initialized;
} kvz_intra_references;

typedef struct
{
  int16_t a;
  int16_t shift;
  int16_t b;
} cclm_parameters_t;

typedef struct {
  int8_t luma_mode;
  int8_t chroma_mode;
  cclm_parameters_t cclm_parameters[2];
  uint8_t multi_ref_idx;
  bool mip_flag;
  bool mip_transp;
  int8_t mts_idx;
  int8_t jccr;
} intra_parameters_t;

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
int8_t kvz_intra_get_dir_luma_predictor(
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
void kvz_intra_build_reference(
  const int_fast8_t log2_width,
  const color_t color,
  const vector2d_t *const luma_px,
  const vector2d_t *const pic_px,
  const lcu_t *const lcu,
  kvz_intra_references *const refs,
  bool entropy_sync,
  kvz_pixel *extra_refs,
  uint8_t multi_ref_idx);

/**
 * \brief Generate intra predictions.
 * \param refs            Reference pixels used for the prediction.
 * \param log2_width      Width of the predicted block.
 * \param mode            Intra mode used for the prediction.
 * \param color           Color of the prediction.
 * \param dst             Buffer for the predicted pixels.
 * \param filter_boundary Whether to filter the boundary on modes 10 and 26.
 */
void kvz_intra_predict(
  encoder_state_t *const state,
  kvz_intra_references *refs,
  int_fast8_t log2_width,
  int_fast8_t mode,
  color_t color,
  kvz_pixel *dst,
  bool filter_boundary,
  const uint8_t multi_ref_idx);

void kvz_intra_recon_cu(
  encoder_state_t *const state,
  int x,
  int y,
  int depth,
  const intra_parameters_t * intra_parameters,
  cu_info_t *cur_cu,
  lcu_t *lcu);


void kvz_predict_cclm(
  encoder_state_t const* const state,
  const color_t color,
  const int8_t width,
  const int8_t height,
  const int16_t x0,
  const int16_t y0,
  const int16_t stride,
  const int8_t mode,
  lcu_t* const lcu,
  kvz_intra_references* chroma_ref,
  kvz_pixel* dst,
  cclm_parameters_t* cclm_params
);

int kvz_get_mip_flag_context(int x, int y, int width, int height, const lcu_t* lcu, cu_array_t* const cu_a);

void kvz_mip_predict(
  encoder_state_t const * const state,
  kvz_intra_references * refs,
  const uint16_t width,
  const uint16_t height,
  kvz_pixel* dst,
  const int mip_mode,
  const bool mip_transp
);