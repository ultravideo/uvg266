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

#include "filter.h"

#include <stdlib.h>

#include "cu.h"
#include "encoder.h"
#include "uvg266.h"
#include "transform.h"
#include "videoframe.h"


//////////////////////////////////////////////////////////////////////////
// INITIALIZATIONS
const uint16_t uvg_g_tc_table_8x8[66] =
{
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   3,   4,   4,   4,
  4,  5,  5,  5,  5,  7,  7,  8,  9,  10,  10,  11,  13,  14,  15,  17,  19,  21,  24,  25,  29,  33,
  36, 41, 45, 51, 57, 64, 71, 80, 89, 100, 112, 125, 141, 157, 177, 198, 222, 250, 280, 314, 352, 395
};

const uint8_t uvg_g_beta_table_8x8[64] =
{
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 22, 24,
  26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56,
  58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88
};

const int8_t uvg_g_luma_filter[16][8] =
{
  {  0, 0,   0, 64,  0,   0,  0,  0 },
  {  0, 1,  -3, 63,  4,  -2,  1,  0 },
  { -1, 2,  -5, 62,  8,  -3,  1,  0 },
  { -1, 3,  -8, 60, 13,  -4,  1,  0 },
  { -1, 4, -10, 58, 17,  -5,  1,  0 }, //1/4
  { -1, 4, -11, 52, 26,  -8,  3, -1 },
  { -1, 3,  -9, 47, 31, -10,  4, -1 },
  { -1, 4, -11, 45, 34, -10,  4, -1 },
  { -1, 4, -11, 40, 40, -11,  4, -1 }, //1/2
  { -1, 4, -10, 34, 45, -11,  4, -1 },
  { -1, 4, -10, 31, 47,  -9,  3, -1 },
  { -1, 3,  -8, 26, 52, -11,  4, -1 },
  {  0, 1,  -5, 17, 58, -10,  4, -1 }, //3/4
  {  0, 1,  -4, 13, 60,  -8,  3, -1 },
  {  0, 1,  -3,  8, 62,  -5,  2, -1 },
  {  0, 1,  -2,  4, 63,  -3,  1,  0 }
};

const int8_t uvg_g_chroma_filter[32][4] =
{
  {  0, 64,  0,  0 },
  { -1, 63,  2,  0 },
  { -2, 62,  4,  0 },
  { -2, 60,  7, -1 },
  { -2, 58, 10, -2 },
  { -3, 57, 12, -2 },
  { -4, 56, 14, -2 },
  { -4, 55, 15, -2 },
  { -4, 54, 16, -2 },
  { -5, 53, 18, -2 },
  { -6, 52, 20, -2 },
  { -6, 49, 24, -3 },
  { -6, 46, 28, -4 },
  { -5, 44, 29, -4 },
  { -4, 42, 30, -4 },
  { -4, 39, 33, -4 },
  { -4, 36, 36, -4 },
  { -4, 33, 39, -4 },
  { -4, 30, 42, -4 },
  { -4, 29, 44, -5 },
  { -4, 28, 46, -6 },
  { -3, 24, 49, -6 },
  { -2, 20, 52, -6 },
  { -2, 18, 53, -5 },
  { -2, 16, 54, -4 },
  { -2, 15, 55, -4 },
  { -2, 14, 56, -4 },
  { -2, 12, 57, -3 },
  { -2, 10, 58, -2 },
  { -1,  7, 60, -2 },
  {  0,  4, 62, -2 },
  {  0,  2, 63, -1 },
};

//////////////////////////////////////////////////////////////////////////
// FUNCTIONS

/**
 * \brief Perform in strong luma filtering in place.
 * \param line  line of 8 pixels, with center at index 4
 * \param tc  tc treshold
 * \return  Reach of the filter starting from center.
 */
static INLINE int uvg_filter_deblock_luma_strong(
    uvg_pixel *line,
    int32_t tc)
{
  const uvg_pixel m0 = line[0];
  const uvg_pixel m1 = line[1];
  const uvg_pixel m2 = line[2];
  const uvg_pixel m3 = line[3];
  const uvg_pixel m4 = line[4];
  const uvg_pixel m5 = line[5];
  const uvg_pixel m6 = line[6];
  const uvg_pixel m7 = line[7];
  const uint8_t tcW[3] = { 3, 2, 1 }; //Wheights for tc

  line[1] = CLIP(m1 - tcW[2]*tc, m1 + tcW[2]*tc, (2*m0 + 3*m1 +   m2 +   m3 +   m4 + 4) >> 3);
  line[2] = CLIP(m2 - tcW[1]*tc, m2 + tcW[1]*tc, (  m1 +   m2 +   m3 +   m4        + 2) >> 2);
  line[3] = CLIP(m3 - tcW[0]*tc, m3 + tcW[0]*tc, (  m1 + 2*m2 + 2*m3 + 2*m4 +   m5 + 4) >> 3);
  line[4] = CLIP(m4 - tcW[0]*tc, m4 + tcW[0]*tc, (  m2 + 2*m3 + 2*m4 + 2*m5 +   m6 + 4) >> 3);
  line[5] = CLIP(m5 - tcW[1]*tc, m5 + tcW[1]*tc, (  m3 +   m4 +   m5 +   m6        + 2) >> 2);
  line[6] = CLIP(m6 - tcW[2]*tc, m6 + tcW[2]*tc, (  m3 +   m4 +   m5 + 3*m6 + 2*m7 + 4) >> 3);

  return 3;
}

/**
 * \brief Perform in weak luma filtering in place.
 * \param encoder  Encoder
 * \param line  Line of 8 pixels, with center at index 4
 * \param tc  The tc treshold
 * \param p_2nd  Whether to filter the 2nd line of P
 * \param q_2nd  Whether to filter the 2nd line of Q
 */
static INLINE int uvg_filter_deblock_luma_weak(
    const encoder_control_t * const encoder,
    uvg_pixel *line,
    int32_t tc,
    bool p_2nd,
    bool q_2nd)
{
  const uvg_pixel m1 = line[1];
  const uvg_pixel m2 = line[2];
  const uvg_pixel m3 = line[3];
  const uvg_pixel m4 = line[4];
  const uvg_pixel m5 = line[5];
  const uvg_pixel m6 = line[6];

  int32_t delta = (9 * (m4 - m3) - 3 * (m5 - m2) + 8) >> 4;

  if (abs(delta) >= tc * 10) {
    return 0;
  } else {
    int32_t tc2 = tc >> 1;
    delta = CLIP(-tc, tc, delta);
    line[3] = CLIP(0, (1 << encoder->bitdepth) - 1, (m3 + delta));
    line[4] = CLIP(0, (1 << encoder->bitdepth) - 1, (m4 - delta));

    if (p_2nd) {
      int32_t delta1 = CLIP(-tc2, tc2, (((m1 + m3 + 1) >> 1) - m2 + delta) >> 1);
      line[2] = CLIP(0, (1 << encoder->bitdepth) - 1, m2 + delta1);
    }
    if (q_2nd) {
      int32_t delta2 = CLIP(-tc2, tc2, (((m6 + m4 + 1) >> 1) - m5 - delta) >> 1);
      line[5] = CLIP(0, (1 << encoder->bitdepth) - 1, m5 + delta2);
    }
    
    if (p_2nd || q_2nd) {
      return 2;
    } else {
      return 1;
    }
  }
}

/**
 * \brief Performe strong/weak filtering for chroma
 */
static INLINE void uvg_filter_deblock_chroma(const encoder_control_t * const encoder,
  uvg_pixel *src,
  int32_t offset,
  int32_t tc,
  int8_t part_P_nofilter,
  int8_t part_Q_nofilter,
  bool sw,
  bool large_boundary,
  bool is_chroma_hor_CTB_boundary)
{
  int32_t delta;
  int16_t m0 = src[-offset * 4];
  int16_t m1 = src[-offset * 3];
  int16_t m2 = src[-offset * 2];
  int16_t m3 = src[-offset];
  int16_t m4 = src[0];
  int16_t m5 = src[offset];
  int16_t m6 = src[offset * 2];
  int16_t m7 = src[offset * 3];

  if (sw) {
    if (is_chroma_hor_CTB_boundary) {
      src[-offset * 1] = CLIP(m3 - tc, m3 + tc, (3 * m2 + 2 * m3 + m4 + m5 + m6 + 4) >> 3);
      src[0] = CLIP(m4 - tc, m4 + tc, (2 * m2 + m3 + 2 * m4 + m5 + m6 + m7 + 4) >> 3);
    } else {
      src[-offset * 3] = CLIP(m1 - tc, m1 + tc, (3 * m0 + 2 * m1 + m2 + m3 + m4 + 4) >> 3);
      src[-offset * 2] = CLIP(m2 - tc, m2 + tc, (2 * m0 + m1 + 2 * m2 + m3 + m4 + m5 + 4) >> 3);
      src[-offset * 1] = CLIP(m3 - tc, m3 + tc, (m0 + m1 + m2 + 2 * m3 + m4 + m5 + m6 + 4) >> 3);
      src[0] = CLIP(m4 - tc, m4 + tc, (m1 + m2 + m3 + 2 * m4 + m5 + m6 + m7 + 4) >> 3);

    }

    src[offset * 1] = CLIP(m5 - tc, m5 + tc, (m2 + m3 + m4 + 2 * m5 + m6 + 2 * m7 + 4) >> 3);
    src[offset * 2] = CLIP(m6 - tc, m6 + tc, (m3 + m4 + m5 + 2 * m6 + 3 * m7 + 4) >> 3);
  } else {
    delta = CLIP(-tc, tc, (((m4 - m3) * 4) + m2 - m5 + 4) >> 3);
    src[-offset] = CLIP(0, (1 << encoder->bitdepth) - 1, m3 + delta);
    src[0] = CLIP(0, (1 << encoder->bitdepth) - 1, m4 - delta);
  }

  if (part_P_nofilter) {
    if (large_boundary) {
      src[-offset * 3] = (uvg_pixel)m1;
      src[-offset * 2] = (uvg_pixel)m2;
    }
    src[-offset * 1] = (uvg_pixel)m3;
  }
  if (part_Q_nofilter) {
    if (large_boundary) {
      src[offset * 1] = (uvg_pixel)m5;
      src[offset * 2] = (uvg_pixel)m6;
    }
    src[0] = (uvg_pixel)m4;
  }
}

/**
 * \brief Check whether an edge is a TU boundary.
 *
 * \param state   encoder state
 * \param x       x-coordinate of the scu in pixels
 * \param y       y-coordinate of the scu in pixels
 * \param dir     direction of the edge to check
 * \return        true, if the edge is a TU boundary, otherwise false
 */
static bool is_tu_boundary(
  const encoder_state_t *const state,
  int32_t x,
  int32_t y,
  edge_dir dir,
  enum uvg_tree_type tree_type)
{
  x >>= tree_type == UVG_CHROMA_T;
  y >>= tree_type == UVG_CHROMA_T;
  // if (x & 3 || y & 3) return false;
  const cu_info_t *const scu =
    uvg_cu_array_at_const(tree_type != UVG_CHROMA_T ? state->tile->frame->cu_array : state->tile->frame->chroma_cu_array, x, y);
  const int tu_width = LCU_WIDTH >> (scu->tr_depth + (tree_type == UVG_CHROMA_T));

  if (dir == EDGE_HOR) {
    return (y & (tu_width - 1)) == 0;
  } else {
    return (x & (tu_width - 1)) == 0;
  }
}


/**
 * \brief Check whether an edge is a PU boundary.
 *
 * \param state   encoder state
 * \param x       x-coordinate of the scu in pixels
 * \param y       y-coordinate of the scu in pixels
 * \param dir     direction of the edge to check
 * \return        true, if the edge is a TU boundary, otherwise false
 */
static bool is_pu_boundary(const encoder_state_t *const state,
                           int32_t x,
                           int32_t y,
                           edge_dir dir)
{
  /*
   TODO: it appears that this function can never be true when is_tu_boundary
   is false. Therefore it should be safe to remove this function but let's keep
   it for now, in case some other tool requires it.
  */
  return false;
  //const cu_info_t *const scu =
  //  uvg_cu_array_at_const(state->tile->frame->cu_array, x, y);
  //// Get the containing CU.
  //const int32_t cu_width = LCU_WIDTH >> scu->depth;
  //const int32_t x_cu = x & ~(cu_width - 1);
  //const int32_t y_cu = y & ~(cu_width - 1);
  //const cu_info_t *const cu =
  //  uvg_cu_array_at_const(state->tile->frame->cu_array, x_cu, y_cu);

  //const int num_pu = uvg_part_mode_num_parts[cu->part_size];
  //for (int i = 0; i < num_pu; i++) {
  //  if (dir == EDGE_HOR) {
  //    int y_pu = PU_GET_Y(cu->part_size, cu_width, y_cu, i);
  //    if (y_pu == y) {
  //      return true;
  //    }

  //  } else {
  //    int x_pu = PU_GET_X(cu->part_size, cu_width, x_cu, i);
  //    if (x_pu == x) {
  //      return true;
  //    }
  //  }
  //}

  //return false;
}


/**
 * \brief Check whether an edge is aligned on a 8x8 grid.
 *
 * \param x     x-coordinate of the edge
 * \param y     y-coordinate of the edge
 * \param dir   direction of the edge
 * \return      true, if the edge is aligned on a 8x8 grid, otherwise false
 */
static bool is_on_8x8_grid(int x, int y, edge_dir dir)
{
  if (dir == EDGE_HOR) {
    return (y & 7) == 0 && (x & 2) == 0;
  } else {
    return (x & 7) == 0 && (y & 2) == 0;
  }
}

static int8_t get_qp_y_pred(const encoder_state_t* state, int x, int y, edge_dir dir)
{
  if (state->frame->max_qp_delta_depth < 0) {
    return state->qp;
  }

  int32_t qp_p;
  if (dir == EDGE_HOR && y > 0) {
    qp_p = uvg_cu_array_at_const(state->tile->frame->cu_array, x, y - 1)->qp;
  } else if (dir == EDGE_VER && x > 0) {
    qp_p = uvg_cu_array_at_const(state->tile->frame->cu_array, x - 1, y)->qp;
  } else {
    // TODO: This seems to be dead code. Investigate.
    qp_p = state->encoder_control->cfg.set_qp_in_cu ? 26 : state->frame->QP;
  }

  const int32_t qp_q =
    uvg_cu_array_at_const(state->tile->frame->cu_array, x, y)->qp;

  return (qp_p + qp_q + 1) >> 1;
}

/**
 * \brief Gather pixels needed for deblocking
 */
static INLINE void gather_deblock_pixels(
    const uvg_pixel *src,
    int step, 
    int stride,
    int reach,
    uvg_pixel *dst)
{
  for (int i = -reach; i < +reach; ++i) {
    dst[i + 4] = src[i * step + stride];
  }
}

/**
* \brief Gather pixels from src to dst using a custom stride and step for src
*/
static INLINE void gather_pixels(
    const uvg_pixel *src,
    int step,
    int stride,
    int numel,
    uvg_pixel *dst)
{
  for (int i = 0; i < numel; ++i) {
    dst[i] = src[i * step + stride];
  }
}

/**
* \brief Scatter pixels
*/
static INLINE void scatter_deblock_pixels(
    const uvg_pixel *src,
    int step, 
    int stride,
    int reach,
    uvg_pixel *dst)
{
  for (int i = -reach; i < +reach; ++i) {
    dst[i * step + stride] = src[i + 4];
  }
}

/**
 * \brief Perform large block strong luma filtering in place.
 * \param line  line of 8 pixels, with center at index 4
 * \param lineL extended pixels with P pixels in [0,3] and Q pixels in [4,7]
 * \param tc  tc treshold
 * \param filter_length_P filter length in the P block
 * \param filter_length_Q filter length in the Q block
 * \return  Reach of the filter starting from center.
 */
static INLINE int uvg_filter_deblock_large_block(uvg_pixel *line, uvg_pixel *lineL, const int32_t tc,
                                                 const uint8_t filter_length_P, const uint8_t filter_length_Q)
{
  int ref_P = 0;
  int ref_Q = 0;
  int ref_middle = 0;

  const int coeffs7[7] = { 59, 50, 41, 32, 23, 14, 5 };
  const int coeffs5[5] = { 58, 45, 32, 19, 6 };
  const int coeffs3[3] = { 53, 32, 11 };

  const int *coeffs_P = NULL;
  const int *coeffs_Q = NULL;

  //Form P/Q arrays that contain all of the samples to make things simpler later
  const uvg_pixel lineP[8] = { line[3], line[2], line[1], line[0],
                               lineL[3], lineL[2], lineL[1], lineL[0] };
  const uvg_pixel lineQ[8] = { line[4], line[5], line[6], line[7],
                               lineL[4], lineL[5], lineL[6], lineL[7] };
  //Separate destination arrays with only six output pixels going in line and  rest to lineL to simplify things later
  uvg_pixel* dstP[7] = { line + 3, line + 2, line + 1,
                         lineL + 3, lineL + 2, lineL + 1, lineL + 0 };
  uvg_pixel* dstQ[7] = { line + 4, line + 5, line + 6,
                         lineL + 4, lineL + 5, lineL + 6, lineL + 7 };

  //Get correct filter coeffs and Q/P end samples
  switch (filter_length_P)
  {
  case 7:
    ref_P = (lineP[6] + lineP[7] + 1) >> 1;
    coeffs_P = coeffs7;
    break;

  case 5:
    ref_P = (lineP[4] + lineP[5] + 1) >> 1;
    coeffs_P = coeffs5;
    break;

  case 3:
    ref_P = (lineP[2] + lineP[3] + 1) >> 1;
    coeffs_P = coeffs3;
    break;
  }

  switch (filter_length_Q)
  {
  case 7:
    ref_Q = (lineQ[6] + lineQ[7] + 1) >> 1;
    coeffs_Q = coeffs7;
    break;

  case 5:
    ref_Q = (lineQ[4] + lineQ[5] + 1) >> 1;
    coeffs_Q = coeffs5;
    break;

  case 3:
    ref_Q = (lineQ[2] + lineQ[3] + 1) >> 1;
    coeffs_Q = coeffs3;
    break;
  }

  //Get middle samples
  if (filter_length_P == filter_length_Q) {
    if (filter_length_P == 7) {
      ref_middle = (lineP[6] + lineP[5] + lineP[4] + lineP[3] + lineP[2] + lineP[1]
                    + 2 * (lineP[0] + lineQ[0])
                    + lineQ[1] + lineQ[2] + lineQ[3] + lineQ[4] + lineQ[5] + lineQ[6] + 8) >> 4;
    }
    else { //filter_length_P == 5
      ref_middle = (lineP[4] + lineP[3]
                    + 2 * (lineP[2] + lineP[1] + lineP[0] + lineQ[0] + lineQ[1] + lineQ[2])
                    + lineQ[3] + lineQ[4] + 8) >> 4;
    }
  }
  else {
    const uint8_t lenS = MIN(filter_length_P, filter_length_Q);
    const uint8_t lenL = MAX(filter_length_P, filter_length_Q);
    const uvg_pixel *refS = filter_length_P < filter_length_Q ? lineP : lineQ;
    const uvg_pixel *refL = filter_length_P < filter_length_Q ? lineQ : lineP;

    if (lenL == 7 && lenS == 5) {
      ref_middle = (lineP[5] + lineP[4] + lineP[3] + lineP[2]
                    + 2 * (lineP[1] + lineP[0] + lineQ[0] + lineQ[1])
                    + lineQ[2] + lineQ[3] + lineQ[4] + lineQ[5] + 8) >> 4;
    }
    else if (lenL == 7 && lenS == 3) {
      ref_middle = (3 * refS[0] + 2 * refL[0] + 3 * refS[1] + refL[1] + 2 * refS[2]
                    + refL[2] + refL[3] + refL[4] + refL[5] + refL[6] + 8) >> 4;
    }
    else { //lenL == 5 && lenS == 3
    ref_middle = (lineP[3] + lineP[2] + lineP[1] + lineP[0]
                  + lineQ[0] + lineQ[1] + lineQ[2] + lineQ[3] + 4) >> 3;

    }
  }

  //Filter pixels in the line

  const uint8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1 };
  const uint8_t tc3[3] = { 6, 4, 2 };

  const uint8_t *tc_coeff_P = (filter_length_P == 3) ? tc3 : tc7;
  const uint8_t *tc_coeff_Q = (filter_length_Q == 3) ? tc3 : tc7;

  for (size_t i = 0; i < filter_length_P; i++)
  {
    int range = (tc * tc_coeff_P[i]) >> 1;
    *dstP[i] = CLIP(lineP[i] - range, lineP[i] + range, (ref_middle * coeffs_P[i] + ref_P * (64 - coeffs_P[i]) + 32) >> 6);
  }

  for (size_t i = 0; i < filter_length_Q; i++)
  {
    int range = (tc * tc_coeff_Q[i]) >> 1;
    *dstQ[i] = CLIP(lineQ[i] - range, lineQ[i] + range, (ref_middle * coeffs_Q[i] + ref_Q * (64 - coeffs_Q[i]) + 32) >> 6);
  }

  return 3;
}

/**
* \brief Determine if strong or weak filtering should be used
*/
static INLINE bool use_strong_filtering(const uvg_pixel * const b0, const uvg_pixel * const b3,
                                        const uvg_pixel * const b0L, const uvg_pixel * const b3L,
                                        const int_fast32_t dp0, const int_fast32_t dq0,
                                        const int_fast32_t dp3, const int_fast32_t dq3,
                                        const int32_t tc, const int32_t beta,
                                        const bool is_side_P_large, const bool is_side_Q_large,
                                        const uint8_t max_filter_length_P, const uint8_t max_filter_length_Q,
                                        const bool is_chroma_CTB_boundary)
{
  int_fast32_t sp0 = is_chroma_CTB_boundary ? abs(b0[2] - b0[3]) : abs(b0[0] - b0[3]);
  int_fast32_t sp3 = is_chroma_CTB_boundary ? abs(b3[2] - b3[3]) : abs(b3[0] - b3[3]);

  if (is_side_P_large || is_side_Q_large) { //Large block decision
    int_fast32_t sq0 = abs(b0[4] - b0[7]);
    int_fast32_t sq3 = abs(b3[4] - b3[7]);
    uvg_pixel tmp0, tmp3;
    if (is_side_P_large) {
      if (max_filter_length_P == 7) {
        tmp0 = b0L[0];
        tmp3 = b3L[0];
        sp0 = sp0 + abs(b0L[3] - b0L[2] - b0L[1] + tmp0);
        sp3 = sp3 + abs(b3L[3] - b3L[2] - b3L[1] + tmp3);
      } else {
        tmp0 = b0L[2];
        tmp3 = b3L[2];
      }
      sp0 = (sp0 + abs(b0[0] - tmp0) + 1) >> 1;
      sp3 = (sp3 + abs(b3[0] - tmp3) + 1) >> 1;
    }
    if (is_side_Q_large) {
      if (max_filter_length_Q == 7) {
        tmp0 = b0L[7];
        tmp3 = b3L[7];
        sq0 = sq0 + abs(b0L[4] - b0L[5] - b0L[6] + tmp0);
        sq3 = sq3 + abs(b3L[4] - b3L[5] - b3L[6] + tmp3);
      } else {
        tmp0 = b0L[5];
        tmp3 = b3L[5];
      }
      sq0 = (sq0 + abs(tmp0 - b0[7]) + 1) >> 1;
      sq3 = (sq3 + abs(tmp3 - b3[7]) + 1) >> 1;
    }
    return 2 * (dp0 + dq0) < beta >> 4 &&
      2 * (dp3 + dq3) < beta >> 4 &&
      abs(b0[3] - b0[4]) < (5 * tc + 1) >> 1 &&
      abs(b3[3] - b3[4]) < (5 * tc + 1) >> 1 &&
      sp0 + sq0 < (beta * 3 >> 5) &&
      sp3 + sq3 < (beta * 3 >> 5);
  } else { //Normal decision
    return 2 * (dp0 + dq0) < beta >> 2 &&
      2 * (dp3 + dq3) < beta >> 2 &&
      abs(b0[3] - b0[4]) < (5 * tc + 1) >> 1 &&
      abs(b3[3] - b3[4]) < (5 * tc + 1) >> 1 &&
      sp0 + abs(b0[4] - b0[7]) < beta >> 3 &&
      sp3 + abs(b3[4] - b3[7]) < beta >> 3;
  }
}

static INLINE void get_max_filter_length(uint8_t *filt_len_P, uint8_t *filt_len_Q,
                                         const encoder_state_t * const state, const uint32_t x, const uint32_t y,
                                         const edge_dir dir, const bool transform_edge,
                                         const int tu_size_P_side, const int tu_size_Q_side,
                                         const int pu_pos, const int pu_size, 
                                         const bool merge_flag, const color_t comp,
                                         enum uvg_tree_type tree_type)
{
  //const int tu_size_P_side = 0;
  //const int tu_size_Q_side = 0;
  //const int size = 0;
  const int x_mul = dir == EDGE_HOR ? 0 : 1;
  const int y_mul = dir == EDGE_HOR ? 1 : 0;
  const int pos = dir == EDGE_HOR ? y : x;
  const int len = EDGE_HOR ? state->tile->frame->height : state->tile->frame->width;
  //const bool transform_edge = is_tu_boundary(state, x, y, dir);
  bool transform_edge_4x4[2] = { false, false };
  bool transform_edge_8x8[2] = { false, false };
  
  if (pos >= 4) transform_edge_4x4[0] = is_tu_boundary(state, x - x_mul * 4, y - y_mul * 4, dir, tree_type);
  if (pos >= 8) transform_edge_8x8[0] = is_tu_boundary(state, x - x_mul * 8, y - y_mul * 8, dir, tree_type);
  if (pos + 4 < len) transform_edge_4x4[1] = is_tu_boundary(state, x + x_mul * 4, y + y_mul * 4, dir, tree_type);
  if (pos + 8 < len) transform_edge_8x8[1] = is_tu_boundary(state, x + x_mul * 8, y + y_mul * 8, dir, tree_type);

  if (comp == COLOR_Y) {
    if (tu_size_P_side <= 4 || tu_size_Q_side <= 4){
      *filt_len_P = 1;
      *filt_len_Q = 1;
    }
    else {
      *filt_len_P = tu_size_P_side >= 32 ? 7 : 3;
      *filt_len_Q = tu_size_Q_side >= 32 ? 7 : 3;
    }

    if ((merge_flag && false) || false) //TODO: Add merge_mode == SUBPU_ATMVP and cu.affine
    {
      if (transform_edge) {
        *filt_len_Q = MIN(*filt_len_Q, 5);
        if (pu_pos > 0) {
          *filt_len_P = MIN(*filt_len_P, 5);
        }
      } else if (pu_pos > 0 && (transform_edge_4x4[0] || (pu_pos + 4) >= pu_size || transform_edge_4x4[1])) { //adjacent to transform edge (4x4 grid)
        *filt_len_P = 1;
        *filt_len_Q = 1;
      } else if (pu_pos > 0 && (pu_pos == 8 || transform_edge_8x8[0] || (pu_pos + 8) >= pu_size || transform_edge_8x8[1])) { //adjacent to transform edge (8x8 grid)
        *filt_len_P = 2;
        *filt_len_Q = 2;
      } else {
        *filt_len_P = 3;
        *filt_len_Q = 3;
      }
    }
  }
  else {
    *filt_len_P = (tu_size_P_side >= 8 && tu_size_Q_side >= 8) ? 3 : 1;
    *filt_len_Q = (tu_size_P_side >= 8 && tu_size_Q_side >= 8) ? 3 : 1;
  }
}

/**
 * \brief Apply the deblocking filter to luma pixels on a single edge.
 *
 * The caller should check that the edge is a TU boundary or a PU boundary.
 *
 \verbatim

         .-- filter this edge if dir == EDGE_HOR
         v
     +--------+
     |o <-- pixel at (x, y)
     |        |
     |<-- filter this edge if dir == EDGE_VER
     |        |
     +--------+

 \endverbatim
 *
 * \param state     encoder state
 * \param x         x-coordinate in pixels (see above)
 * \param y         y-coordinate in pixels (see above)
 * \param length    length of the edge in pixels
 * \param dir       direction of the edge to filter
 * \param tu_boundary   whether the edge is a TU boundary
 */
static void filter_deblock_edge_luma(encoder_state_t * const state,
                                     int32_t x,
                                     int32_t y,
                                     int32_t length,
                                     edge_dir dir,
                                     bool tu_boundary)
{
  videoframe_t * const frame = state->tile->frame;
  const encoder_control_t * const encoder = state->encoder_control;
  
  {
    int32_t stride = frame->rec->stride;
    int32_t beta_offset_div2 = encoder->cfg.deblock_beta;
    int32_t tc_offset_div2   = encoder->cfg.deblock_tc;
    // TODO: support 10+bits
    uvg_pixel *orig_src = &frame->rec->y[x + y*stride];
    uvg_pixel *src = orig_src;

    const int32_t qp = get_qp_y_pred(state, x, y, dir);

    const int MAX_QP = 63; //TODO: Make DEFAULT_INTRA_TC_OFFSET(=2) a define?
    const int8_t lumaBitdepth = encoder->bitdepth;

    int8_t strength = 0;
    int32_t bitdepth_scale  = 1 << (lumaBitdepth - 8);
    int32_t b_index         = CLIP(0, MAX_QP, qp + (beta_offset_div2 << 1));
    int32_t beta            = uvg_g_beta_table_8x8[b_index] * bitdepth_scale;
    int32_t side_threshold  = (beta + (beta >>1 )) >> 3;
    int32_t tc_index;
    int32_t tc;

    //Deblock adapted to halve pixel mvd.
    const int16_t mvdThreashold = 1 << (INTERNAL_MV_PREC - 1);

    uint32_t num_4px_parts  = length / 4;

    // Transpose the image by swapping x and y strides when doing horizontal
    // edges.
    const int32_t x_stride = (dir == EDGE_VER) ? 1 : stride;
    const int32_t y_stride = (dir == EDGE_VER) ? stride : 1;

    // TODO: add CU based QP calculation

    // For each 4-pixel part in the edge
    for (uint32_t block_idx = 0; block_idx < num_4px_parts; ++block_idx) {
      
      // CUs on both sides of the edge
      cu_info_t *cu_p;
      cu_info_t *cu_q;
      int32_t y_coord = y;
      int32_t x_coord = x;
      {
        if (dir == EDGE_VER) {
          y_coord = y + 4 * block_idx;
          cu_p = uvg_cu_array_at(frame->cu_array, x - 1, y_coord);
          cu_q = uvg_cu_array_at(frame->cu_array, x, y_coord);

        } else {
          x_coord = x + 4 * block_idx;
          cu_p = uvg_cu_array_at(frame->cu_array, x_coord, y - 1);
          cu_q = uvg_cu_array_at(frame->cu_array, x_coord, y);
        }

        bool nonzero_coeffs = cbf_is_set(cu_q->cbf, cu_q->tr_depth, COLOR_Y)
          || cbf_is_set(cu_p->cbf, cu_p->tr_depth, COLOR_Y);

        // Filter strength
        strength = 0;
        if (cu_q->type == CU_INTRA || cu_p->type == CU_INTRA) { // Intra is used
          strength = 2;
        }
        else if (tu_boundary && nonzero_coeffs) {
          // Non-zero residual/coeffs and transform boundary
          // Neither CU is intra so tr_depth <= MAX_DEPTH.
          strength = 1;
        }
        else if(cu_p->inter.mv_dir == 3 || cu_q->inter.mv_dir == 3 || state->frame->slicetype == UVG_SLICE_B) { // B-slice related checks. TODO: Need to account for cu_p being in another slice?

          // Zero all undefined motion vectors for easier usage
          if(!(cu_q->inter.mv_dir & 1)) {
            cu_q->inter.mv[0][0] = 0;
            cu_q->inter.mv[0][1] = 0;
          }
          if(!(cu_q->inter.mv_dir & 2)) {
            cu_q->inter.mv[1][0] = 0;
            cu_q->inter.mv[1][1] = 0;
          }

          if(!(cu_p->inter.mv_dir & 1)) {
            cu_p->inter.mv[0][0] = 0;
            cu_p->inter.mv[0][1] = 0;
          }
          if(!(cu_p->inter.mv_dir & 2)) {
            cu_p->inter.mv[1][0] = 0;
            cu_p->inter.mv[1][1] = 0;
          }
          const int refP0 = (cu_p->inter.mv_dir & 1) ? state->frame->ref_LX[0][cu_p->inter.mv_ref[0]] : -1;
          const int refP1 = (cu_p->inter.mv_dir & 2) ? state->frame->ref_LX[1][cu_p->inter.mv_ref[1]] : -1;
          const int refQ0 = (cu_q->inter.mv_dir & 1) ? state->frame->ref_LX[0][cu_q->inter.mv_ref[0]] : -1;
          const int refQ1 = (cu_q->inter.mv_dir & 2) ? state->frame->ref_LX[1][cu_q->inter.mv_ref[1]] : -1;
          const mv_t* mvQ0 = cu_q->inter.mv[0];
          const mv_t* mvQ1 = cu_q->inter.mv[1];

          const mv_t* mvP0 = cu_p->inter.mv[0];
          const mv_t* mvP1 = cu_p->inter.mv[1];

          if(( refP0 == refQ0 &&  refP1 == refQ1 ) || ( refP0 == refQ1 && refP1==refQ0 ))
          {
            // Different L0 & L1
            if ( refP0 != refP1 ) {          
              if ( refP0 == refQ0 ) {
                strength  = ((abs(mvQ0[0] - mvP0[0]) >= mvdThreashold) ||
                             (abs(mvQ0[1] - mvP0[1]) >= mvdThreashold) ||
                             (abs(mvQ1[0] - mvP1[0]) >= mvdThreashold) ||
                             (abs(mvQ1[1] - mvP1[1]) >= mvdThreashold)) ? 1 : 0;
              } else {
                strength  = ((abs(mvQ1[0] - mvP0[0]) >= mvdThreashold) ||
                             (abs(mvQ1[1] - mvP0[1]) >= mvdThreashold) ||
                             (abs(mvQ0[0] - mvP1[0]) >= mvdThreashold) ||
                             (abs(mvQ0[1] - mvP1[1]) >= mvdThreashold)) ? 1 : 0;
              }
            // Same L0 & L1
            } else {  
              strength  = ((abs(mvQ0[0] - mvP0[0]) >= mvdThreashold) ||
                           (abs(mvQ0[1] - mvP0[1]) >= mvdThreashold) ||
                           (abs(mvQ1[0] - mvP1[0]) >= mvdThreashold) ||
                           (abs(mvQ1[1] - mvP1[1]) >= mvdThreashold)) &&
                          ((abs(mvQ1[0] - mvP0[0]) >= mvdThreashold) ||
                           (abs(mvQ1[1] - mvP0[1]) >= mvdThreashold) ||
                           (abs(mvQ0[0] - mvP1[0]) >= mvdThreashold) ||
                           (abs(mvQ0[1] - mvP1[1]) >= mvdThreashold)) ? 1 : 0;
            }
          } else {
            strength = 1;
          }
        }
        else /*if (cu_p->inter.mv_dir != 3 && cu_q->inter.mv_dir != 3)*/ { //is P-slice
          if (cu_q->inter.mv_ref[cu_q->inter.mv_dir - 1] != cu_p->inter.mv_ref[cu_p->inter.mv_dir - 1]) {
            // Reference pictures are different
            strength = 1;
          } else if (
            ((abs(cu_q->inter.mv[cu_q->inter.mv_dir - 1][0] - cu_p->inter.mv[cu_p->inter.mv_dir - 1][0]) >= mvdThreashold) ||
            (abs(cu_q->inter.mv[cu_q->inter.mv_dir - 1][1] - cu_p->inter.mv[cu_p->inter.mv_dir - 1][1]) >= mvdThreashold))) {
            // Absolute motion vector diff between blocks >= 0.5 (Integer pixel)
            strength = 1;
          }
        }
      
        tc_index        = CLIP(0, MAX_QP + 2, (int32_t)(qp + 2*(strength - 1) + (tc_offset_div2 << 1)));
        tc              = lumaBitdepth < 10 ? ((uvg_g_tc_table_8x8[tc_index] + (1 << (9 - lumaBitdepth))) >> (10 - lumaBitdepth))
                                            : ((uvg_g_tc_table_8x8[tc_index] << (lumaBitdepth - 10)));
      }

      if (strength == 0) continue;

      bool is_side_P_large = false;
      bool is_side_Q_large = false;
      uint8_t max_filter_length_P = 0;
      uint8_t max_filter_length_Q = 0;
      const int cu_size = LCU_WIDTH >> cu_q->depth;
      const int pu_part_idx = (y + PU_GET_H(cu_q->part_size, cu_size, 0) <= y_coord ? 
                               1 + (uvg_part_mode_num_parts[cu_q->part_size] >> 2) : 0)
                            + (x + PU_GET_W(cu_q->part_size, cu_size, 0) <= x_coord ? 1 : 0);
      const int pu_size = dir == EDGE_HOR ? PU_GET_H(cu_q->part_size, cu_size, pu_part_idx)
                                          : PU_GET_W(cu_q->part_size, cu_size, pu_part_idx);
      const int pu_pos = dir == EDGE_HOR ? y_coord - PU_GET_Y(cu_q->part_size, cu_size, 0, pu_part_idx) 
                                         : x_coord - PU_GET_X(cu_q->part_size, cu_size, 0, pu_part_idx);
      get_max_filter_length(&max_filter_length_P, &max_filter_length_Q, state, x_coord, y_coord,
                            dir, tu_boundary,
                            LCU_WIDTH >> cu_p->tr_depth,
                            LCU_WIDTH >> cu_q->tr_depth,
                            pu_pos, pu_size, cu_q->merged, COLOR_Y,
                            UVG_LUMA_T);

      if (max_filter_length_P > 3) {
        is_side_P_large = dir == EDGE_HOR && y % LCU_WIDTH == 0 ? false : true;
        //TODO: Add affine/ATMVP related stuff
        /*if (max_filter_length_P > 5 && cu_p->affine) {
          max_filter_length_P = MIN(max_filter_length_P, 5);
        }*/
      }
      if (max_filter_length_Q > 3) {
        is_side_Q_large = true;
      }
      //                               +-- edge_src
      //                               v
      // line0 p7 p6 p5 p4 p3 p2 p1 p0 q0 q1 q2 q3 q4 q5 q6 q7
      uvg_pixel *edge_src = &src[block_idx * 4 * y_stride];

      // Gather the lines of pixels required for the filter on/off decision.
      //TODO: May need to limit reach in small blocks?
      uvg_pixel b[4][8];
      gather_deblock_pixels(edge_src, x_stride, 0 * y_stride, 4, &b[0][0]);
      gather_deblock_pixels(edge_src, x_stride, 3 * y_stride, 4, &b[3][0]);

      int_fast32_t dp0 = abs(b[0][1] - 2 * b[0][2] + b[0][3]);
      int_fast32_t dq0 = abs(b[0][4] - 2 * b[0][5] + b[0][6]);
      int_fast32_t dp3 = abs(b[3][1] - 2 * b[3][2] + b[3][3]);
      int_fast32_t dq3 = abs(b[3][4] - 2 * b[3][5] + b[3][6]);
      int_fast32_t dp = dp0 + dp3;
      int_fast32_t dq = dq0 + dq3;

      bool sw = false;

      if (is_side_P_large || is_side_Q_large) {
        int_fast32_t dp0L = dp0;
        int_fast32_t dq0L = dq0;
        int_fast32_t dp3L = dp3;
        int_fast32_t dq3L = dq3;
        
        //In case of large blocks, need to gather extra pixels
        //bL:
        //line0 p7 p6 p5 p4 q4 q5 q6 q7
        uvg_pixel bL[4][8];

        if (is_side_P_large) {
          gather_pixels(edge_src - 8 * x_stride, x_stride, 0 * y_stride, 4, &bL[0][0]);
          gather_pixels(edge_src - 8 * x_stride, x_stride, 3 * y_stride, 4, &bL[3][0]);
          dp0L = (dp0L + abs(bL[0][2] - 2 * bL[0][3] + b[0][0]) + 1) >> 1;
          dp3L = (dp3L + abs(bL[3][2] - 2 * bL[3][3] + b[3][0]) + 1) >> 1;
        }
        if (is_side_Q_large) {
          gather_pixels(edge_src + 4 * x_stride, x_stride, 0 * y_stride, 4, &bL[0][4]);
          gather_pixels(edge_src + 4 * x_stride, x_stride, 3 * y_stride, 4, &bL[3][4]);
          dq0L = (dq0L + abs(b[0][7] - 2 * bL[0][4] + bL[0][5]) + 1) >> 1;
          dq3L = (dq3L + abs(b[3][7] - 2 * bL[3][4] + bL[3][5]) + 1) >> 1;
        }
        
        int_fast32_t dpL = dp0L + dp3L;
        int_fast32_t dqL = dq0L + dq3L;

        if (dpL + dqL < beta) {
          sw = use_strong_filtering(&b[0][0], &b[3][0], &bL[0][0], &bL[3][0],
                                    dp0L, dq0L, dp3L, dq3L, tc, beta,
                                    is_side_P_large, is_side_Q_large,
                                    max_filter_length_P, max_filter_length_Q, false);
          if (sw) {
            gather_deblock_pixels(edge_src, x_stride, 1 * y_stride, 4, &b[1][0]);
            gather_deblock_pixels(edge_src, x_stride, 2 * y_stride, 4, &b[2][0]);
            if (is_side_P_large)
            {
              gather_pixels(edge_src - 8 * x_stride, x_stride, 1 * y_stride, 4, &bL[1][0]);
              gather_pixels(edge_src - 8 * x_stride, x_stride, 2 * y_stride, 4, &bL[2][0]);
            }
            if (is_side_Q_large)
            {
              gather_pixels(edge_src + 4 * x_stride, x_stride, 1 * y_stride, 4, &bL[1][4]);
              gather_pixels(edge_src + 4 * x_stride, x_stride, 2 * y_stride, 4, &bL[2][4]);
            }

            for (int i = 0; i < 4; ++i) {
              int filter_reach;
              filter_reach = uvg_filter_deblock_large_block(&b[i][0], &bL[i][0], tc,
                                                            is_side_P_large ? max_filter_length_P : 3, 
                                                            is_side_Q_large ? max_filter_length_Q : 3);
              scatter_deblock_pixels(&b[i][0], x_stride, i * y_stride, filter_reach, edge_src);
              if (is_side_P_large) {
                const int diff_reach = (max_filter_length_P - filter_reach) >> 1;
                const int dst_offset = (filter_reach + diff_reach) * x_stride;
                scatter_deblock_pixels(&bL[i][0] - diff_reach, x_stride, i * y_stride, diff_reach, edge_src - dst_offset);
              }
              if (is_side_Q_large) {
                const int diff_reach = (max_filter_length_Q - filter_reach) >> 1;
                const int dst_offset = (filter_reach + diff_reach) * x_stride;
                scatter_deblock_pixels(&bL[i][0] + diff_reach, x_stride, i * y_stride, diff_reach, edge_src + dst_offset);
              }
            }
          }
        }
      }

      if (!sw)
      {
        if (dp + dq < beta) {
          if (max_filter_length_P > 2 && max_filter_length_Q > 2) {
            // Strong filtering flag checking.
            sw = use_strong_filtering(b[0], b[3], NULL, NULL,
                                      dp0, dq0, dp3, dq3, tc, beta,
                                      false, false, 7, 7, false);
          }

          // Read lines 1 and 2. Weak filtering doesn't use the outermost pixels
          // but let's give them anyway to simplify control flow.
          gather_deblock_pixels(edge_src, x_stride, 1 * y_stride, 4, &b[1][0]);
          gather_deblock_pixels(edge_src, x_stride, 2 * y_stride, 4, &b[2][0]);

          for (int i = 0; i < 4; ++i) {
            int filter_reach;
            if (sw) {
              filter_reach = uvg_filter_deblock_luma_strong(&b[i][0], tc);
            } else {
              bool p_2nd = false;
              bool q_2nd = false;
              if (max_filter_length_P > 1 && max_filter_length_Q > 1) {
                p_2nd = dp < side_threshold;
                q_2nd = dq < side_threshold;
              }
              filter_reach = uvg_filter_deblock_luma_weak(encoder, &b[i][0], tc, p_2nd, q_2nd);
            }
            scatter_deblock_pixels(&b[i][0], x_stride, i * y_stride, filter_reach, edge_src);
          }
        }
      }
    }
  }
}

/**
 * \brief Apply the deblocking filter to chroma pixels on a single edge.
 *
 * The caller should check that the edge is a TU boundary or a PU boundary.
 *
 \verbatim

         .-- filter this edge if dir == EDGE_HOR
         v
     +--------+
     |o <-- pixel at (x, y)
     |        |
     |<-- filter this edge if dir == EDGE_VER
     |        |
     +--------+

 \endverbatim
 *
 * \param state         encoder state
 * \param x             x-coordinate in chroma pixels (see above)
 * \param y             y-coordinate in chroma pixels (see above)
 * \param length        length of the edge in chroma pixels
 * \param dir           direction of the edge to filter
 * \param tu_boundary   whether the edge is a TU boundary
 */
static void filter_deblock_edge_chroma(encoder_state_t * const state,
                                       int32_t x,
                                       int32_t y,
                                       int32_t length,
                                       edge_dir dir,
                                       bool tu_boundary,
                                       enum uvg_tree_type tree_type)
{
  const encoder_control_t * const encoder = state->encoder_control;
  const videoframe_t * const frame = state->tile->frame;
   
  // For each subpart
  {
    int32_t stride = frame->rec->stride >> 1;
    int32_t tc_offset_div2 = encoder->cfg.deblock_tc;
    int32_t beta_offset_div2 = encoder->cfg.deblock_beta;
    // TODO: support 10+bits
    uvg_pixel *src[] = {
      &frame->rec->u[x + y*stride],
      &frame->rec->v[x + y*stride],
    };
    
    const uint8_t MAX_QP = 63;

    const int32_t luma_qp  = get_qp_y_pred(state, x << 1, y << 1, dir);
    int32_t QP = uvg_get_scaled_qp(1, luma_qp, 0, state->encoder_control->qp_map[0]);//uvg_g_chroma_scale[luma_qp]; //TODO: Add BDOffset?
    int32_t bitdepth_scale = 1 << (encoder->bitdepth - 8);
    
    //TU size should be in chroma samples (?)
    const int chroma_shift = dir == EDGE_HOR ? (encoder->chroma_format == UVG_CSP_420 ? 1 : 0)
                                             : (encoder->chroma_format != UVG_CSP_444 ? 1 : 0);
    //TODO: Replace two (2) with min CU log2 size when its updated to the correct value 
    const int min_chroma_width_log2 = 2-(encoder->chroma_format == UVG_CSP_420 ? 1 : 0);
    const int min_chroma_height_log2 = 2 -(encoder->chroma_format != UVG_CSP_444 ? 1 : 0);
    const int min_chroma_size_log2 = dir == EDGE_HOR ? min_chroma_height_log2 : min_chroma_width_log2;
    const int min_chroma_length = 1 << min_chroma_size_log2;
    const uint32_t num_parts = (length >> min_chroma_size_log2);
    
    const int32_t offset = (dir == EDGE_HOR) ? stride :      1;
    const int32_t step   = (dir == EDGE_HOR) ?      1 : stride;

    for (uint32_t blk_idx = 0; blk_idx < num_parts; ++blk_idx)
    {
      // CUs on both sides of the edge
      cu_info_t *cu_p;
      cu_info_t *cu_q;
      int32_t x_coord = x << (tree_type != UVG_CHROMA_T);
      int32_t y_coord = y << (tree_type != UVG_CHROMA_T);
      cu_array_t* cua = tree_type != UVG_CHROMA_T ? frame->cu_array : frame->chroma_cu_array;
      if (dir == EDGE_VER) {
        y_coord = (y + min_chroma_length * blk_idx) << (tree_type != UVG_CHROMA_T);
        cu_p = uvg_cu_array_at(cua, x_coord - 1, y_coord);
        cu_q = uvg_cu_array_at(cua, x_coord    , y_coord);

      } else {
        x_coord = (x + min_chroma_length * blk_idx) << (tree_type != UVG_CHROMA_T);
        cu_p = uvg_cu_array_at(cua, x_coord, y_coord - 1);
        cu_q = uvg_cu_array_at(cua, x_coord, y_coord    );
      }

      const int cu_size = LCU_WIDTH >> (cu_q->depth + (tree_type == UVG_CHROMA_T));
      const int pu_part_idx = ((y << (tree_type != UVG_CHROMA_T)) + PU_GET_H(cu_q->part_size, cu_size, 0) <= y_coord ?
                               1 + (uvg_part_mode_num_parts[cu_q->part_size] >> 2) : 0)
                              + ((x << (tree_type != UVG_CHROMA_T)) + PU_GET_W(cu_q->part_size, cu_size, 0) <= x_coord ? 1 : 0);
      const int pu_size = dir == EDGE_HOR ? PU_GET_H(cu_q->part_size, cu_size, pu_part_idx)
                                          : PU_GET_W(cu_q->part_size, cu_size, pu_part_idx);
      const int pu_pos = dir == EDGE_HOR ? y_coord - PU_GET_Y(cu_q->part_size, cu_size, 0, pu_part_idx)
                                         : x_coord - PU_GET_X(cu_q->part_size, cu_size, 0, pu_part_idx);
      uint8_t max_filter_length_P = 0;
      uint8_t max_filter_length_Q = 0;
      
      const int tu_p_size = LCU_WIDTH >> (cu_p->tr_depth + (chroma_shift));
      const int tu_q_size = LCU_WIDTH >> (cu_q->tr_depth + (chroma_shift));
      get_max_filter_length(&max_filter_length_P, &max_filter_length_Q, state, x_coord, y_coord,
                            dir, tu_boundary, tu_p_size, tu_q_size,
                            pu_pos, pu_size, cu_q->merged, COLOR_U,
                            tree_type);


      const bool large_boundary = (max_filter_length_P >= 3 && max_filter_length_Q >= 3);
      const bool is_chroma_hor_CTB_boundary = (dir == EDGE_HOR && y_coord % (LCU_WIDTH >> (tree_type == UVG_CHROMA_T)) == 0);
      uint8_t c_strength[2] = { 0, 0 };
      

      if (cu_q->type == CU_INTRA || cu_p->type == CU_INTRA) {
        c_strength[0] = 2;
        c_strength[1] = 2;
      }
      else if (tu_boundary){ //TODO: Add ciip/IBC related stuff
        bool nonzero_coeffs_U = cbf_is_set(cu_q->cbf, cu_q->tr_depth, COLOR_U)
                                || cbf_is_set(cu_p->cbf, cu_p->tr_depth, COLOR_U);
        bool nonzero_coeffs_V = cbf_is_set(cu_q->cbf, cu_q->tr_depth, COLOR_V)
                                || cbf_is_set(cu_p->cbf, cu_p->tr_depth, COLOR_V);
        c_strength[0] = nonzero_coeffs_U ? 1 : 0;
        c_strength[1] = nonzero_coeffs_V ? 1 : 0;
      }

      for (int component = 0; component < 2; component++) {
        if (c_strength[component] == 2 || (large_boundary && c_strength[component] == 1)){

          int32_t TC_index = CLIP(0, MAX_QP + 2, (int32_t)(QP + 2 * (c_strength[component] - 1) + (tc_offset_div2 << 1)));
          int32_t Tc = encoder->bitdepth < 10 ? ((uvg_g_tc_table_8x8[TC_index] + (1 << (9 - encoder->bitdepth))) >> (10 - encoder->bitdepth))
                                              : (uvg_g_tc_table_8x8[TC_index] << (encoder->bitdepth - 10));

          bool use_long_filter = false;

          //                   +-- edge_src
          //                   v
          // line0 p3 p2 p1 p0 q0 q1 q2 q3
          uvg_pixel *edge_src = &src[component][blk_idx * min_chroma_length * step];
            
          if (large_boundary) {
            const int beta_index = CLIP(0, MAX_QP, QP + (beta_offset_div2 << 1));
            const int beta = uvg_g_beta_table_8x8[beta_index] * bitdepth_scale;


            const uint8_t sss = chroma_shift == 1 ? 1 : 3;
            // Gather the lines of pixels required for the filter on/off decision.
            //TODO: May need to limit reach in small blocks?
            uvg_pixel b[2][8];
            gather_deblock_pixels(edge_src, offset, 0 * step, 4, &b[0][0]);
            gather_deblock_pixels(edge_src, offset, sss * step, 4, &b[1][0]);

            const uint8_t p_ind = is_chroma_hor_CTB_boundary ? 2 : 1;

            int_fast32_t dp0 = abs(b[0][p_ind] - 2 * b[0][2] + b[0][3]);
            int_fast32_t dq0 = abs(b[0][4] - 2 * b[0][5] + b[0][6]);
            int_fast32_t dp3 = abs(b[1][p_ind] - 2 * b[1][2] + b[1][3]);
            int_fast32_t dq3 = abs(b[1][4] - 2 * b[1][5] + b[1][6]);
            int_fast32_t dp = dp0 + dp3;
            int_fast32_t dq = dq0 + dq3;

            if (dp + dq < beta) {
              use_long_filter = true;
              const bool sw = use_strong_filtering(b[0], b[1], NULL, NULL,
                                                   dp0, dq0, dp3, dq3, Tc, beta,
                                                   false, false, 7, 7, is_chroma_hor_CTB_boundary);
              for (int i = 0; i < min_chroma_length; i++) {
                uvg_filter_deblock_chroma(encoder, edge_src + step * i, offset, Tc, 0, 0,
                                          sw, large_boundary, is_chroma_hor_CTB_boundary);
              }
            }
          }
          if (!use_long_filter)
          {
            for (int i = 0; i < min_chroma_length; i++) {
              uvg_filter_deblock_chroma(encoder, edge_src + step * i, offset, Tc, 0, 0,
                                        false, large_boundary, is_chroma_hor_CTB_boundary);
            }
          }
        }
      }
    }
  }
}


/**
 * \brief Filter edge of a single PU or TU
 *
 * \param state     encoder state
 * \param x         block x-position in pixels
 * \param y         block y-position in pixels
 * \param width     block width in pixels
 * \param height    block height in pixels
 * \param dir       direction of the edges to filter
 * \param tu_boundary   whether the edge is a TU boundary
 */
static void filter_deblock_unit(
  encoder_state_t * const state,
  int x,
  int y,
  int width,
  int height,
  edge_dir dir,
  bool tu_boundary,
  bool previous_ctu,
  enum uvg_tree_type tree_type)
{
  // no filtering on borders (where filter would use pixels outside the picture)
  if (x == 0 && dir == EDGE_VER) return;
  if (y == 0 && dir == EDGE_HOR) return;

  // Length of luma and chroma edges.
  int32_t length;
  if (dir == EDGE_HOR) {
    const videoframe_t* const frame = state->tile->frame;
    const int32_t x_right = x + width;
    const bool rightmost_8px_of_lcu = x_right % LCU_WIDTH == 0 || x_right % LCU_WIDTH == LCU_WIDTH - width;
    const bool rightmost_8px_of_frame = x_right == frame->width || x_right + width == frame->width;

    if (rightmost_8px_of_lcu && !rightmost_8px_of_frame && !previous_ctu) {
      // The last 8 pixels will be deblocked when processing the next LCU.
      length   = width - 4;
      if (length == 0) return;

    } else {
      length   = width;
    }
  } else {
    length   = height;
  }

  if(tree_type != UVG_CHROMA_T) {
    filter_deblock_edge_luma(state, x, y, length, dir, tu_boundary);
  }

  // Chroma pixel coordinates.
  const int32_t x_c = x >> 1;
  const int32_t y_c = y >> 1;
  if (state->encoder_control->chroma_format != UVG_CSP_400 && 
    (is_on_8x8_grid(x_c, y_c, dir && (x_c + 4) % 32)
     || (x == state->tile->frame->width - 8 && dir == 1 && y_c % 8 == 0)) 
    && tree_type != UVG_LUMA_T) {
    filter_deblock_edge_chroma(state, x_c, y_c, length, dir, tu_boundary, tree_type);
  }
}


/**
 * \brief Deblock PU and TU boundaries inside an LCU.
 *
 * \param state     encoder state
 * \param x_px      block x-position in pixels
 * \param y_px      block y-position in pixels
 * \param dir       direction of the edges to filter
 *
 * Recursively traverse the CU/TU quadtree. At the lowest level, apply the
 * deblocking filter to the left edge (when dir == EDGE_VER) or the top edge
 * (when dir == EDGE_HOR) as needed. Both luma and chroma are filtered.
 */
static void filter_deblock_lcu_inside(encoder_state_t * const state,
                                      int32_t x,
                                      int32_t y,
                                      edge_dir dir)
{
  const int end_x = MIN(x + LCU_WIDTH, state->tile->frame->width);
  const int end_y = MIN(y + LCU_WIDTH, state->tile->frame->height);

  const enum uvg_tree_type luma_tree = state->frame->is_irap && state->encoder_control->cfg.dual_tree ? UVG_LUMA_T : UVG_BOTH_T;
  const enum uvg_tree_type chroma_tree = state->frame->is_irap && state->encoder_control->cfg.dual_tree ? UVG_CHROMA_T : UVG_BOTH_T;

  for (int edge_y = y; edge_y < end_y; edge_y += 4) {
    for (int edge_x = x; edge_x < end_x; edge_x += 4) {
      bool tu_boundary = is_tu_boundary(state, edge_x, edge_y, dir, luma_tree);
      if (tu_boundary || is_pu_boundary(state, edge_x, edge_y, dir)) {
        filter_deblock_unit(state, edge_x, edge_y, 4, 4, dir, tu_boundary, edge_x < x, luma_tree);
      }
      if(chroma_tree == UVG_CHROMA_T && is_tu_boundary(state, edge_x, edge_y, dir, chroma_tree)) {
        filter_deblock_unit(state, edge_x, edge_y, 4, 4, dir, tu_boundary, edge_x < x, chroma_tree);        
      }
    }
  }
}


/**
 * \brief Filter rightmost 8 pixels of the horizontal egdes of an LCU.
 *
 * \param state     encoder state
 * \param x_px      x-coordinate of the *right* edge of the LCU in pixels
 * \param y_px      y-coordinate of the top edge of the LCU in pixels
 */
static void filter_deblock_lcu_rightmost(encoder_state_t * const state,
                                         int32_t x_px,
                                         int32_t y_px)
{
  // Luma
  const enum uvg_tree_type luma_tree = state->frame->is_irap && state->encoder_control->cfg.dual_tree ? UVG_LUMA_T : UVG_BOTH_T;
  const enum uvg_tree_type chroma_tree = state->frame->is_irap && state->encoder_control->cfg.dual_tree ? UVG_CHROMA_T : UVG_BOTH_T;

  const int end = MIN(y_px + LCU_WIDTH, state->tile->frame->height);
  for (int x = x_px - 8; x < x_px; x += 4) {
    for (int y = y_px; y < end; y += 4) {
      // The top edge of the whole frame is not filtered.
      bool tu_boundary = is_tu_boundary(state, x, y, EDGE_HOR, luma_tree);
      if (y > 0 && (tu_boundary || is_pu_boundary(state, x, y, EDGE_HOR))) {
        filter_deblock_edge_luma(state, x, y, 4, EDGE_HOR, tu_boundary);
      }
    }
  }

  // Chroma
  if (state->encoder_control->chroma_format != UVG_CSP_400) {
    const int x_px_c = x_px >> 1;
    const int y_px_c = y_px >> 1;
    const int x_c = x_px_c - 4;
    const int end_c = MIN(y_px_c + LCU_WIDTH_C, state->tile->frame->height >> 1);
    for (int y_c = y_px_c; y_c < end_c; y_c += 8) {
      // The top edge of the whole frame is not filtered.
      bool tu_boundary = is_tu_boundary(state, x_c << 1, y_c << 1, EDGE_HOR, chroma_tree);
      if (y_c > 0 && (tu_boundary || is_pu_boundary(state, x_c << 1, y_c << 1, EDGE_HOR))) {
        filter_deblock_edge_chroma(state, x_c , y_c, 4, EDGE_HOR, tu_boundary, chroma_tree);
      }
    }
  }
}


/**
 * \brief Deblock a single LCU without using data from right or down.
 *
 * Filter the following vertical edges (horizontal filtering):
 *  1. The left edge of the LCU.
 *  2. All vertical edges within the LCU.
 *
 * Filter the following horizontal edges (vertical filtering):
 *  1. The rightmost 8 pixels of the top edge of the LCU to the left.
 *  2. The rightmost 8 pixels of all horizontal edges within the LCU to the
 *     left.
 *  3. The top edge and all horizontal edges within the LCU, excluding the
 *     rightmost 8 pixels. If the LCU is the rightmost LCU of the frame, the
 *     last 8 pixels are also filtered.
 *
 * What is not filtered:
 *  - The rightmost 8 pixels of the top edge and all horizontal edges within
 *    the LCU, unless the LCU is the rightmost LCU of the frame.
 *  - The bottom edge of the LCU.
 *  - The right edge of the LCU.
 *
 * \param state   encoder state
 * \param x_px    x-coordinate of the left edge of the LCU in pixels
 * \param y_px    y-coordinate of the top edge of the LCU in pixels
 */
 //TODO: Things to check/fix for VVC:
// - Strength calculation to include average Luma level (Luma Adaptive Deblocing Filter LADF) (optional)
// - Deblocking strength for CIIP and IBC modes (CIIP/IBC not currently used)
// - Handle new prediction modes (i.e. PLT) (PLT not currently used)
// - Deblocking filter for subblock boundaries
// - Allow loop filtering across slice/tile boundaries?
void uvg_filter_deblock_lcu(encoder_state_t * const state, int x_px, int y_px)
{
  assert(!state->encoder_control->cfg.lossless);
  filter_deblock_lcu_inside(state, x_px, y_px, EDGE_VER);
  if (x_px > 0) {
    filter_deblock_lcu_rightmost(state, x_px, y_px);
  }
  filter_deblock_lcu_inside(state, x_px, y_px, EDGE_HOR);
}
