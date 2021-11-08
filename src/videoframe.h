#ifndef VIDEOFRAME_H_
#define VIDEOFRAME_H_
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
 * \ingroup DataStructures
 * \file
 * \brief Container for the frame currently being encoded.
 */

#include "cu.h"
#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"


/**
 * \brief Struct which contains all picture data
 */
typedef struct videoframe
{
  kvz_picture *source;         //!< \brief Source image.
  kvz_picture *source_lmcs;    //!< \brief LMCS mapped source image if available, otherwise points to source.
  kvz_picture *rec;            //!< \brief Reconstructed image.
  kvz_picture *rec_lmcs;       //!< \brief LMCS mapped reconstructed image, if available, otherwise points to source.

  uint8_t* lmcs_avg_processed; //!< \brief For each LCU, indicates if already calculated average of border pixels is available
  int32_t* lmcs_avg;           //!< \brief Average of LCU border pixels

  int32_t width;          //!< \brief Luma pixel array width.
  int32_t height;         //!< \brief Luma pixel array height.
  int32_t height_in_lcu;  //!< \brief Picture width in number of LCU's.
  int32_t width_in_lcu;   //!< \brief Picture height in number of LCU's.

  cu_array_t* cu_array;     //!< \brief Info for each CU at each depth.
  struct lmcs_aps* lmcs_aps; //!< \brief LMCS parameters for both the current frame.
  struct sao_info_t *sao_luma;   //!< \brief Array of sao parameters for every LCU.
  struct sao_info_t *sao_chroma;   //!< \brief Array of sao parameters for every LCU.
  struct alf_info_t *alf_info;   //!< \brief Array of alf parameters for both luma and chroma.
  struct param_set_map* alf_param_set_map;

  int32_t poc;           //!< \brief Picture order count
  cu_info_t* hmvp_lut; //!< \brief Look-up table for HMVP, one for each LCU row

  uint8_t* hmvp_size; //!< \brief HMVP LUT size
  bool source_lmcs_mapped; //!< \brief Indicate if source_lmcs is available and mapped to LMCS
  bool lmcs_top_level; //!< \brief Indicate that in this level the LMCS images are allocated
  bool rec_lmcs_mapped; //!< \brief Indicate if rec_lmcs is available and mapped to LMCS
} videoframe_t;


videoframe_t *kvz_videoframe_alloc(int32_t width, int32_t height, enum kvz_chroma_format chroma_format, enum kvz_alf alf_type);
int kvz_videoframe_free(videoframe_t * const frame);

void kvz_videoframe_set_poc(videoframe_t * frame, int32_t poc);

#endif
