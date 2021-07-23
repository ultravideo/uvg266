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
  const int_fast8_t channel_type,
  const kvz_pixel *const in_ref_above,
  const kvz_pixel *const in_ref_left,
  kvz_pixel *const dst)
{
  
  assert(log2_width >= 2 && log2_width <= 5);
  assert(intra_mode >= 2 && intra_mode <= 66);

  static const int16_t modedisp2sampledisp[32] = { 0,    1,    2,    3,    4,    6,     8,   10,   12,   14,   16,   18,   20,   23,   26,   29,   32,   35,   39,  45,  51,  57,  64,  73,  86, 102, 128, 171, 256, 341, 512, 1024 };
  static const int16_t modedisp2invsampledisp[32] = { 0, 16384, 8192, 5461, 4096, 2731, 2048, 1638, 1365, 1170, 1024, 910, 819, 712, 630, 565, 512, 468, 420, 364, 321, 287, 256, 224, 191, 161, 128, 96, 64, 48, 32, 16 }; // (512 * 32) / sampledisp
  static const int32_t pre_scale[] = { 8, 7, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, -1, -1, -2, -3 };
  static const int16_t intraGaussFilter[32][4] = {
  { 16, 32, 16, 0 },
  { 15, 29, 17, 3 },
  { 15, 29, 17, 3 },
  { 14, 29, 18, 3 },
  { 13, 29, 18, 4 },
  { 13, 28, 19, 4 },
  { 13, 28, 19, 4 },
  { 12, 28, 20, 4 },
  { 11, 28, 20, 5 },
  { 11, 27, 21, 5 },
  { 10, 27, 22, 5 },
  { 9, 27, 22, 6 },
  { 9, 26, 23, 6 },
  { 9, 26, 23, 6 },
  { 8, 25, 24, 7 },
  { 8, 25, 24, 7 },
  { 8, 24, 24, 8 },
  { 7, 24, 25, 8 },
  { 7, 24, 25, 8 },
  { 6, 23, 26, 9 },
  { 6, 23, 26, 9 },
  { 6, 22, 27, 9 },
  { 5, 22, 27, 10 },
  { 5, 21, 27, 11 },
  { 5, 20, 28, 11 },
  { 4, 20, 28, 12 },
  { 4, 19, 28, 13 },
  { 4, 19, 28, 13 },
  { 4, 18, 29, 13 },
  { 3, 18, 29, 14 },
  { 3, 17, 29, 15 },
  { 3, 17, 29, 15 }
  };
  static const int16_t cubic_filter[32][4] =
  {
    { 0, 64,  0,  0 },
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
    { 0,  4, 62, -2 },
    { 0,  2, 63, -1 },
  };

                                                    // Temporary buffer for modes 11-25.
                                                    // It only needs to be big enough to hold indices from -width to width-1.
  kvz_pixel temp_main[2 * 128] = { 0 };
  kvz_pixel temp_side[2 * 128] = { 0 };
  const int_fast32_t width = 1 << log2_width;

  uint32_t pred_mode = intra_mode; // ToDo: handle WAIP

  // Whether to swap references to always project on the left reference row.
  const bool vertical_mode = intra_mode >= 34;
  // Modes distance to horizontal or vertical mode.
  const int_fast8_t mode_disp = vertical_mode ? pred_mode - 50 : -(pred_mode - 18);
  //const int_fast8_t mode_disp = vertical_mode ? intra_mode - 26 : 10 - intra_mode;
  
  // Sample displacement per column in fractions of 32.
  const int_fast8_t sample_disp = (mode_disp < 0 ? -1 : 1) * modedisp2sampledisp[abs(mode_disp)];
  
  // TODO: replace latter width with height
  int scale = MIN(2, log2_width - pre_scale[abs(mode_disp)]);

  // Pointer for the reference we are interpolating from.
  kvz_pixel *ref_main;
  // Pointer for the other reference.
  const kvz_pixel *ref_side;

  // Set ref_main and ref_side such that, when indexed with 0, they point to
  // index 0 in block coordinates.
  if (sample_disp < 0) {
    for (int i = 0; i <= width + 1; i++) {
      temp_main[width + i] = (vertical_mode ? in_ref_above[i] : in_ref_left[i]);
      temp_side[width + i] = (vertical_mode ? in_ref_left[i] : in_ref_above[i]);
    }

    ref_main = temp_main + width;
    ref_side = temp_side + width;

    for (int i = -width; i <= -1; i++) {
      ref_main[i] = ref_side[MIN((-i * modedisp2invsampledisp[abs(mode_disp)] + 256) >> 9, width)];
    }

    

    //const uint32_t index_offset = width + 1;
    //const int32_t last_index = width;
    //const int_fast32_t most_negative_index = (width * sample_disp) >> 5;
    //// Negative sample_disp means, we need to use both references.

    //// TODO: update refs to take into account variating block size and shapes
    ////       (height is not always equal to width)
    //ref_side = (vertical_mode ? in_ref_left : in_ref_above) + 1;
    //ref_main = (vertical_mode ? in_ref_above : in_ref_left) + 1;

    //// Move the reference pixels to start from the middle to the later half of
    //// the tmp_ref, so there is room for negative indices.
    //for (int_fast32_t x = -1; x < width; ++x) {
    //  tmp_ref[x + index_offset] = ref_main[x];
    //}
    //// Get a pointer to block index 0 in tmp_ref.
    //ref_main = &tmp_ref[index_offset];
    //tmp_ref[index_offset -1] = tmp_ref[index_offset];

    //// Extend the side reference to the negative indices of main reference.
    //int_fast32_t col_sample_disp = 128; // rounding for the ">> 8"
    //int_fast16_t inv_abs_sample_disp = modedisp2invsampledisp[abs(mode_disp)];
    //// TODO: add 'vertical_mode ? height : width' instead of 'width'
    //
    //for (int_fast32_t x = -1; x > most_negative_index; x--) {
    //  col_sample_disp += inv_abs_sample_disp;
    //  int_fast32_t side_index = col_sample_disp >> 8;
    //  tmp_ref[x + index_offset - 1] = ref_side[side_index - 1];
    //}
    //tmp_ref[last_index + index_offset] = tmp_ref[last_index + index_offset - 1];
    //tmp_ref[most_negative_index + index_offset - 1] = tmp_ref[most_negative_index + index_offset];
  }
  else {
    
    for (int i = 0; i <= (width << 1); i++) {
      temp_main[i] = (vertical_mode ? in_ref_above[i] : in_ref_left[i]);
      temp_side[i] = (vertical_mode ? in_ref_left[i] : in_ref_above[i]);
    }

    const int s = 0;
    const int max_index = (0 << s) + 2;
    const int ref_length = width << 1;
    const kvz_pixel val = temp_main[ref_length];
    for (int j = 0; j <= max_index; j++) {
      temp_main[ref_length + j] = val;
    }

    ref_main = temp_main;
    ref_side = temp_side;
    //// sample_disp >= 0 means we don't need to refer to negative indices,
    //// which means we can just use the references as is.
    //ref_main = (vertical_mode ? in_ref_above : in_ref_left) + 1;
    //ref_side = (vertical_mode ? in_ref_left : in_ref_above) + 1;

    //memcpy(tmp_ref + width, ref_main, (width*2) * sizeof(kvz_pixel));
    //ref_main = &tmp_ref[width];
    //tmp_ref[width-1] = tmp_ref[width];
    //int8_t last_index = 1 + width*2;
    //tmp_ref[width + last_index] = tmp_ref[width + last_index - 1];
  }

  if (sample_disp != 0) {
    // The mode is not horizontal or vertical, we have to do interpolation.

    int_fast32_t delta_pos = 0;
    for (int_fast32_t y = 0; y < width; ++y) {
      delta_pos += sample_disp;
      int_fast32_t delta_int = delta_pos >> 5;
      int_fast32_t delta_fract = delta_pos & (32 - 1);

      if ((abs(sample_disp) & 0x1F) != 0) {
        
        // Luma Channel
        if (channel_type == 0) {
          int32_t ref_main_index = delta_int;
          kvz_pixel p[4];
          bool use_cubic = true; // Default to cubic filter
          static const int kvz_intra_hor_ver_dist_thres[8] = { 24, 24, 24, 14, 2, 0, 0, 0 };
          int filter_threshold = kvz_intra_hor_ver_dist_thres[log2_width];
          int dist_from_vert_or_hor = MIN(abs(pred_mode - 50), abs(pred_mode - 18));
          if (dist_from_vert_or_hor > filter_threshold) {
            static const int16_t modedisp2sampledisp[32] = { 0,    1,    2,    3,    4,    6,     8,   10,   12,   14,   16,   18,   20,   23,   26,   29,   32,   35,   39,  45,  51,  57,  64,  73,  86, 102, 128, 171, 256, 341, 512, 1024 };
            const int_fast8_t mode_disp = (pred_mode >= 34) ? pred_mode - 50 : 18 - pred_mode;
            const int_fast8_t sample_disp = (mode_disp < 0 ? -1 : 1) * modedisp2sampledisp[abs(mode_disp)];
            if ((abs(sample_disp) & 0x1F) != 0)
            {
              use_cubic = false;
            }
          }
          const int16_t filter_coeff[4] = { 16 - (delta_fract >> 1), 32 - (delta_fract >> 1), 16 + (delta_fract >> 1), delta_fract >> 1 };
          int16_t const * const f = use_cubic ? cubic_filter[delta_fract] : filter_coeff;
          // Do 4-tap intra interpolation filtering
          for (int_fast32_t x = 0; x < width; x++, ref_main_index++) {
            p[0] = ref_main[ref_main_index];
            p[1] = ref_main[ref_main_index + 1];
            p[2] = ref_main[ref_main_index + 2];
            p[3] = ref_main[ref_main_index + 3];
         
            dst[y * width + x] = CLIP_TO_PIXEL(((int32_t)(f[0] * p[0]) + (int32_t)(f[1] * p[1]) + (int32_t)(f[2] * p[2]) + (int32_t)(f[3] * p[3]) + 32) >> 6);

          }
        }
        else {
        
          // Do linear filtering
          for (int_fast32_t x = 0; x < width; ++x) {
            kvz_pixel ref1 = ref_main[x + delta_int + 1];
            kvz_pixel ref2 = ref_main[x + delta_int + 2];
            dst[y * width + x] = ref1 + ((delta_fract * (ref2-ref1) + 16) >> 5);
          }
        }
      }
      else {
        // Just copy the integer samples
        for (int_fast32_t x = 0; x < width; x++) {
          dst[y * width + x] = ref_main[x + delta_int + 1];
        }
      }

     
      // PDPC
      bool PDPC_filter = (width >= 4 || channel_type != 0);
      if (pred_mode > 1 && pred_mode < 67) {
        if (mode_disp < 0) {
          PDPC_filter = false;
        }
        else if (mode_disp > 0) {
          PDPC_filter = (scale >= 0);
        }
      }
      if(PDPC_filter) {
        int       inv_angle_sum = 256;
        for (int x = 0; x < MIN(3 << scale, width); x++) {
          inv_angle_sum += modedisp2invsampledisp[abs(mode_disp)];

          int wL = 32 >> (2 * x >> scale);
          const kvz_pixel left = ref_side[y + (inv_angle_sum >> 9) + 1];
          dst[y * width + x] = dst[y * width + x] + ((wL * (left - dst[y * width + x]) + 32) >> 6);
        }
      }

        /*
      if (pred_mode == 2 || pred_mode == 66) {
        int wT = 16 >> MIN(31, ((y << 1) >> scale));
        for (int x = 0; x < width; x++) {
          int wL = 16 >> MIN(31, ((x << 1) >> scale));
          if (wT + wL == 0) break;
          int c = x + y + 1;
          if (c >= 2 * width) { wL = 0; }
          if (c >= 2 * width) { wT = 0; }
          const kvz_pixel left = (wL != 0) ? ref_side[c] : 0;
          const kvz_pixel top  = (wT != 0) ? ref_main[c] : 0;
          dst[y * width + x] = CLIP_TO_PIXEL((wL * left + wT * top + (64 - wL - wT) * dst[y * width + x] + 32) >> 6);
        }
      } else if (sample_disp == 0 || sample_disp >= 12) {
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
          const kvz_pixel *p = ref_side + delta_y - 1;
          kvz_pixel left = p[delta_frac_0 >> 5];
          dst[y * width + x] = CLIP_TO_PIXEL((wL * left + (64 - wL) * dst[y * width + x] + 32) >> 6);
        }
      }*/
    }
  }
  else {
    // Mode is horizontal or vertical, just copy the pixels.

    // TODO: update outer loop to use height instead of width
    for (int_fast32_t y = 0; y < width; ++y) {
      for (int_fast32_t x = 0; x < width; ++x) {
        dst[y * width + x] = ref_main[x + 1];
      }
      if ((width >= 4 || channel_type != 0) && sample_disp >= 0) {
        int scale = (log2_width + log2_width - 2) >> 2;
        const kvz_pixel top_left = ref_main[0];
        const kvz_pixel left = ref_side[1 + y];
        for (int i = 0; i < MIN(3 << scale, width); i++) {
          const int wL = 32 >> (2 * i >> scale);
          const kvz_pixel val = dst[y * width + i];
          dst[y * width + i] = CLIP_TO_PIXEL(val + ((wL * (left - top_left) + 32) >> 6));
        }
      }
    }
  }

  // Flip the block if this is was a horizontal mode.
  if (!vertical_mode) {
    for (int_fast32_t y = 0; y < width - 1; ++y) {
      for (int_fast32_t x = y + 1; x < width; ++x) {
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
  // TODO: Add height
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
      //
    }
  }
#endif
}

/**
* \brief Generage intra DC prediction with post filtering applied.
* \param log2_width    Log2 of width, range 2..5.
* \param in_ref_above  Pointer to -1 index of above reference, length=width*2+1.
* \param in_ref_left   Pointer to -1 index of left reference, length=width*2+1.
* \param dst           Buffer of size width*width.
*/
static void kvz_intra_pred_filtered_dc_generic(
  const int_fast8_t log2_width,
  const kvz_pixel *const ref_top,
  const kvz_pixel *const ref_left,
  kvz_pixel *const out_block)
{
  assert(log2_width >= 2 && log2_width <= 5);

  const int_fast8_t width = 1 << log2_width;

  int_fast16_t sum = 0;
  for (int_fast8_t i = 0; i < width; ++i) {
    sum += ref_top[i + 1];
    sum += ref_left[i + 1];
  }

  const kvz_pixel dc_val = (sum + width) >> (log2_width + 1);

  // Filter top-left with ([1 2 1] / 4)
  out_block[0] = (ref_left[1] + 2 * dc_val + ref_top[1] + 2) / 4;

  // Filter rest of the boundary with ([1 3] / 4)
  for (int_fast8_t x = 1; x < width; ++x) {
    out_block[x] = (ref_top[x + 1] + 3 * dc_val + 2) / 4;
  }
  for (int_fast8_t y = 1; y < width; ++y) {
    out_block[y * width] = (ref_left[y + 1] + 3 * dc_val + 2) / 4;
    for (int_fast8_t x = 1; x < width; ++x) {
      out_block[y * width + x] = dc_val;
    }
  }
}


int kvz_strategy_register_intra_generic(void* opaque, uint8_t bitdepth)
{
  bool success = true;

  success &= kvz_strategyselector_register(opaque, "angular_pred", "generic", 0, &kvz_angular_pred_generic);
  success &= kvz_strategyselector_register(opaque, "intra_pred_planar", "generic", 0, &kvz_intra_pred_planar_generic);
  success &= kvz_strategyselector_register(opaque, "intra_pred_filtered_dc", "generic", 0, &kvz_intra_pred_filtered_dc_generic);

  return success;
}
