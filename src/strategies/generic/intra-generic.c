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

#include "strategies/generic/intra-generic.h"

#include <stdlib.h>

#include "kvazaar.h"
#include "strategyselector.h"
#include "kvz_math.h"


 /**
 * \brief Generage angular predictions.
 * \param log2_width    Log2 of width, range 2..5.
 * \param intra_mode    Angular mode in range 2..34.
 * \param in_ref_above  Pointer to -1 index of above reference, length=width*2+1.
 * \param in_ref_left   Pointer to -1 index of left reference, length=width*2+1.
 * \param dst           Buffer of size width*width.
 */
static void kvz_angular_pred_generic(
  const int_fast8_t log2_width,
  const int_fast8_t intra_mode,
  const kvz_pixel *const in_ref_above,
  const kvz_pixel *const in_ref_left,
  kvz_pixel *const dst)
{
  
  assert(log2_width >= 2 && log2_width <= 5);
  assert(intra_mode >= 2 && intra_mode <= 66);

  static const int8_t modedisp2sampledisp[27] = { 0, 1, 2, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 26, 29, 32, 35, 39, 45, 49, 54, 60, 68, 79, 93, 114 };
  static const int16_t modedisp2invsampledisp[27] = { 0, 8192, 4096, 2731, 1638, 1170, 910, 745, 630, 546, 482, 431, 390, 356, 315, 282, 256, 234, 210, 182, 167, 152, 137, 120, 104, 88, 72 }; // (256 * 32) / sampledisp

                                                    // Temporary buffer for modes 11-25.
                                                    // It only needs to be big enough to hold indices from -width to width-1.
  kvz_pixel tmp_ref[2 * 32];
  const int_fast8_t width = 1 << log2_width;

  // Whether to swap references to always project on the left reference row.
  const bool vertical_mode = intra_mode >= 34;
  // Modes distance to horizontal or vertical mode.
  const int_fast8_t mode_disp = vertical_mode ? intra_mode - 50 : 18 - intra_mode;
  //const int_fast8_t mode_disp = vertical_mode ? intra_mode - 26 : 10 - intra_mode;
  
  // Sample displacement per column in fractions of 32.
  const int_fast8_t sample_disp = (mode_disp < 0 ? -1 : 1) * modedisp2sampledisp[abs(mode_disp)];

  // Pointer for the reference we are interpolating from.
  const kvz_pixel *ref_main;
  // Pointer for the other reference.
  const kvz_pixel *ref_side;

  // Set ref_main and ref_side such that, when indexed with 0, they point to
  // index 0 in block coordinates.
  if (sample_disp < 0) {
    // Negative sample_disp means, we need to use both references.

    // TODO: update refs to take into account variating block size and shapes
    //       (height is not always equal to width)
    ref_side = (vertical_mode ? in_ref_left : in_ref_above) + 1;
    ref_main = (vertical_mode ? in_ref_above : in_ref_left) + 1;

    // Move the reference pixels to start from the middle to the later half of
    // the tmp_ref, so there is room for negative indices.
    for (int_fast8_t x = -1; x < width; ++x) {
      tmp_ref[x + width] = ref_main[x];
    }
    // Get a pointer to block index 0 in tmp_ref.
    ref_main = &tmp_ref[width];

    // Extend the side reference to the negative indices of main reference.
    int_fast32_t col_sample_disp = 128; // rounding for the ">> 8"
    int_fast16_t inv_abs_sample_disp = modedisp2invsampledisp[abs(mode_disp)];
    // TODO: add 'vertical_mode ? height : width' instead of 'width'
    int_fast8_t most_negative_index = (width * sample_disp) >> 5;
    for (int_fast8_t x = -2; x >= most_negative_index; --x) {
      col_sample_disp += inv_abs_sample_disp;
      int_fast8_t side_index = col_sample_disp >> 8;
      tmp_ref[x + width] = ref_side[side_index - 1];
    }
  }
  else {
    // sample_disp >= 0 means we don't need to refer to negative indices,
    // which means we can just use the references as is.
    ref_main = (vertical_mode ? in_ref_above : in_ref_left) + 1;
    ref_side = (vertical_mode ? in_ref_left : in_ref_above) + 1;
  }

  if (sample_disp != 0) {
    // The mode is not horizontal or vertical, we have to do interpolation.

    int_fast16_t delta_pos = 0;
    for (int_fast8_t y = 0; y < width; ++y) {
      delta_pos += sample_disp;
      int_fast8_t delta_int = delta_pos >> 5;
      int_fast8_t delta_fract = delta_pos & (32 - 1);

      if (delta_fract) {
        // Do linear filtering
        for (int_fast8_t x = 0; x < width; ++x) {
          // TODO: is the ref1 same for the first iterator as in VTM?
          kvz_pixel ref1 = ref_main[x + delta_int];
          kvz_pixel ref2 = ref_main[x + delta_int + 1];
          dst[y * width + x] = ((32 - delta_fract) * ref1 + delta_fract * ref2 + 16) >> 5;
        }
      }
      else {
        // Just copy the integer samples
        for (int_fast8_t x = 0; x < width; x++) {
          // TODO: check if [x + delta_int] or [x + delta_int + 1]
          dst[y * width + x] = ref_main[x + delta_int + 1];
        }
      }

      // TODO: replace latter width with height
      int scale = ((kvz_math_floor_log2(width) - 2 + kvz_math_floor_log2(width) - 2 + 2) >> 2);

      if (sample_disp == 2 || sample_disp == 66) {
        int wT = 16 >> MIN(31, ((y << 1) >> scale));
        for (int x = 0; x < width; x++) {
          int wL = 16 >> MIN(31, ((x << 1) >> scale));
          if (wT + wL == 0) break;
          int c = x + y + 1;
          const kvz_pixel left = (wL != 0) ? ref_side[c + 1] : 0;
          const kvz_pixel top  = (wT != 0) ? ref_main[c + 1] : 0;
          dst[y * width + x] = CLIP_TO_PIXEL((wL * left + wT * top + (64 - wL - wT) * dst[y * width + x] + 32) >> 6);
        }
      }
      else if ((sample_disp >= 58) || sample_disp <= (10)) {
        int inv_angle_sum_0 = 2;
        for (int x = 0; x < width; x++) {
          inv_angle_sum_0 += modedisp2invsampledisp[abs(mode_disp)];
          int delta_pos_0 = inv_angle_sum_0 >> 2;
          int delta_frac_0 = delta_pos_0 & 63;
          int delta_int_0 = delta_pos_0 >> 6;
          int delta_y = y + delta_int_0 + 1;
          // TODO: convert to JVET_K0500_WAIP
          if (delta_y > width + width - 1) break;

          int wL = 32 >> MIN(31, ((x << 1) >> scale));
          if (wL == 0) break;
          kvz_pixel *p = ref_side + delta_y;
          kvz_pixel left = (((64 - delta_frac_0) * p[0] + delta_frac_0 * p[1] + 32) >> 6);
          dst[y * width + x] = CLIP_TO_PIXEL((wL * left + (64 - wL) * dst[y * width + x] + 32) >> 6);
        }
      }
    }
  }
  else {
    // Mode is horizontal or vertical, just copy the pixels.

    // TODO: update outer loop to use height instead of width
    for (int_fast8_t y = 0; y < width; ++y) {
      for (int_fast8_t x = 0; x < width; ++x) {
        dst[y * width + x] = ref_main[x+1];
      }
    }
  }

  // Flip the block if this is was a horizontal mode.
  if (!vertical_mode) {
    for (int_fast8_t y = 0; y < width - 1; ++y) {
      for (int_fast8_t x = y + 1; x < width; ++x) {
        SWAP(dst[y * width + x], dst[x * width + y], kvz_pixel);
      }
    }
  }
}


/**
 * \brief Generate planar prediction.
 * \param log2_width    Log2 of width, range 2..5.
 * \param in_ref_above  Pointer to -1 index of above reference, length=width*2+1.
 * \param in_ref_left   Pointer to -1 index of left reference, length=width*2+1.
 * \param dst           Buffer of size width*width.
 */
static void kvz_intra_pred_planar_generic(
  const int_fast8_t log2_width,
  const kvz_pixel *const ref_top,
  const kvz_pixel *const ref_left,
  kvz_pixel *const dst)
{
  assert(log2_width >= 2 && log2_width <= 5);

  const int_fast8_t width = 1 << log2_width;
  const kvz_pixel top_right = ref_top[width + 1];
  const kvz_pixel bottom_left = ref_left[width + 1];

#if 0
  // Unoptimized version for reference.
  for (int y = 0; y < width; ++y) {
    for (int x = 0; x < width; ++x) {
      int_fast16_t hor = (width - 1 - x) * ref_left[y + 1] + (x + 1) * top_right;
      int_fast16_t ver = (width - 1 - y) * ref_top[x + 1] + (y + 1) * bottom_left;
      dst[y * width + x] = (ver + hor + width) >> (log2_width + 1);
    }
  }
#else
  int_fast16_t top[32];
  for (int i = 0; i < width; ++i) {
    top[i] = ref_top[i + 1] << log2_width;
  }

  for (int y = 0; y < width; ++y) {
    int_fast16_t hor = (ref_left[y + 1] << log2_width) + width;
    for (int x = 0; x < width; ++x) {
      hor += top_right - ref_left[y + 1];
      top[x] += bottom_left - ref_top[x + 1];
      dst[y * width + x] = (hor + top[x]) >> (log2_width + 1);
    }
  }
#endif
}

int kvz_strategy_register_intra_generic(void* opaque, uint8_t bitdepth)
{
  bool success = true;

  success &= kvz_strategyselector_register(opaque, "angular_pred", "generic", 0, &kvz_angular_pred_generic);
  success &= kvz_strategyselector_register(opaque, "intra_pred_planar", "generic", 0, &kvz_intra_pred_planar_generic);

  return success;
}
