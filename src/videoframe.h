#ifndef VIDEOFRAME_H_
#define VIDEOFRAME_H_
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
 * \ingroup DataStructures
 * \file
 * \brief Container for the frame currently being encoded.
 */

#include "cu.h"
#include "global.h" // IWYU pragma: keep
#include "uvg266.h"
#include "hashmap.h"


/**
 * \brief Struct which contains all picture data
 */
typedef struct videoframe
{
  uvg_picture *source;         //!< \brief Source image.
  uvg_picture *source_lmcs;    //!< \brief LMCS mapped source image if available, otherwise points to source.
  uvg_picture *rec;            //!< \brief Reconstructed image.
  uvg_picture *rec_lmcs;       //!< \brief LMCS mapped reconstructed image, if available, otherwise points to source.

  uvg_pixel *cclm_luma_rec;    //!< \brief buffer for the downsampled luma reconstruction for cclm
  uvg_pixel *cclm_luma_rec_top_line;    //!< \brief buffer for the downsampled luma reconstruction for cclm

  uint8_t* lmcs_avg_processed; //!< \brief For each LCU, indicates if already calculated average of border pixels is available
  int32_t* lmcs_avg;           //!< \brief Average of LCU border pixels

  int32_t width;          //!< \brief Luma pixel array width.
  int32_t height;         //!< \brief Luma pixel array height.
  int32_t height_in_lcu;  //!< \brief Picture width in number of LCU's.
  int32_t width_in_lcu;   //!< \brief Picture height in number of LCU's.
  int32_t chroma_width;          //!< \brief Luma pixel array width.
  int32_t chroma_height;         //!< \brief Luma pixel array height.
  int32_t chroma_height_in_lcu;  //!< \brief Picture width in number of LCU's.
  int32_t chroma_width_in_lcu;   //!< \brief Picture height in number of LCU's.

  cu_array_t* cu_array;     //!< \brief Info for each CU at each depth.
  cu_array_t* chroma_cu_array;     //!< \brief Info for each CU at each depth.
  struct lmcs_aps* lmcs_aps; //!< \brief LMCS parameters for both the current frame.
  struct sao_info_t *sao_luma;   //!< \brief Array of sao parameters for every LCU.
  struct sao_info_t *sao_chroma;   //!< \brief Array of sao parameters for every LCU.
  struct alf_info_t *alf_info;   //!< \brief Array of alf parameters for both luma and chroma.
  struct param_set_map* alf_param_set_map;

  int32_t poc;           //!< \brief Picture order count
    
  uvg_pixel **ibc_buffer_y; //!< \brief Intra Block Copy buffer for each LCU row 
  uvg_pixel **ibc_buffer_u; //!< \brief Intra Block Copy buffer for each LCU row 
  uvg_pixel **ibc_buffer_v; //!< \brief Intra Block Copy buffer for each LCU row
  uvg_hashmap_t **ibc_hashmap_row; //!< \brief Hashmap for IBC hash search for each LCU row
  uint32_t *ibc_hashmap_pos_to_hash; //!< \brief Hashmap reverse search for position to hash
  uint32_t ibc_hashmap_pos_to_hash_stride; //!< \brief Hashmap position to hash stride
  cu_info_t* hmvp_lut_ibc; //!< \brief Look-up table for HMVP in IBC, one for each LCU row
  uint8_t* hmvp_size_ibc; //!< \brief HMVP IBC LUT size

  cu_info_t* hmvp_lut; //!< \brief Look-up table for HMVP, one for each LCU row
  uint8_t* hmvp_size; //!< \brief HMVP LUT size
  bool source_lmcs_mapped; //!< \brief Indicate if source_lmcs is available and mapped to LMCS
  bool lmcs_top_level; //!< \brief Indicate that in this level the LMCS images are allocated
  bool rec_lmcs_mapped; //!< \brief Indicate if rec_lmcs is available and mapped to LMCS

} videoframe_t;


videoframe_t *uvg_videoframe_alloc(int32_t width, int32_t height, enum uvg_chroma_format chroma_format, enum uvg_alf alf_type, bool cclm);
int uvg_videoframe_free(videoframe_t * const frame);

void uvg_videoframe_set_poc(videoframe_t *const frame, const int32_t poc);

#endif
