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

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "scalinglist.h"

#include "uvg_math.h"
#include "rdo.h"
#include "tables.h"


const uint8_t uvg_g_scaling_list_num[SCALING_LIST_SIZE_NUM] = {6, 6, 6, 2};
const uint16_t uvg_g_scaling_list_size[SCALING_LIST_SIZE_NUM][SCALING_LIST_SIZE_NUM] =
{
  {1, 2, 4, 8, 16, 32, 64, 128,},
  {2, 4, 8, 16, 32, 64, 128, 256,},
  {4, 8, 16, 32, 64, 128, 256, 512,},
  {8, 16, 32, 64, 128, 256, 512, 1024,},
  {16, 32, 64, 128, 256, 512, 1024, 2048,},
  {32, 64, 128, 256, 512, 1024, 2048, 4096,},
  {64, 128, 256, 512, 1024, 2048, 4096, 8192,},
  {128, 256, 512, 1024, 2048, 4096, 8192, 16384,},
};
static const uint8_t g_scaling_list_size_x[SCALING_LIST_SIZE_NUM] = {1, 2, 4, 8, 16, 32, 64, 128};

static const int32_t g_quant_default_4x4[16] =
{
  16, 16, 16, 16,
  16, 16, 16, 16,
  16, 16, 16, 16,
  16, 16, 16, 16
};

// ToDo: check these
static const int32_t g_quant_intra_default_8x8[64] =
{
  16, 16, 16, 16, 17, 18, 21, 24,
  16, 16, 16, 16, 17, 19, 22, 25,
  16, 16, 17, 18, 20, 22, 25, 29,
  16, 16, 18, 21, 24, 27, 31, 36,
  17, 17, 20, 24, 30, 35, 41, 47,
  18, 19, 22, 27, 35, 44, 54, 65,
  21, 22, 25, 31, 41, 54, 70, 88,
  24, 25, 29, 36, 47, 65, 88, 115
};

static const int32_t g_quant_inter_default_8x8[64] =
{
  16, 16, 16, 16, 17, 18, 20, 24,
  16, 16, 16, 17, 18, 20, 24, 25,
  16, 16, 17, 18, 20, 24, 25, 28,
  16, 17, 18, 20, 24, 25, 28, 33,
  17, 18, 20, 24, 25, 28, 33, 41,
  18, 20, 24, 25, 28, 33, 41, 54,
  20, 24, 25, 28, 33, 41, 54, 71,
  24, 25, 28, 33, 41, 54, 71, 91
};

const int16_t uvg_g_quant_scales[6] = {26214, 23302, 20560, 18396, 16384, 14564};
const int16_t uvg_g_inv_quant_scales[6] = {40, 45, 51, 57, 64, 72};


/**
 * \brief Initialize scaling lists
 *
 */
void uvg_scalinglist_init(scaling_list_t* const scaling_list) {
  uint32_t size_id_y, size_id_x, list_id, qp;
  for (size_id_x = 0; size_id_x < SCALING_LIST_SIZE_NUM; size_id_x++) {
    for (size_id_y = 0; size_id_y < SCALING_LIST_SIZE_NUM; size_id_y++) {
      for (list_id = 0; list_id < SCALING_LIST_NUM; list_id++) {
        for (qp = 0; qp < 6; qp++) {
          scaling_list->quant_coeff[size_id_x][size_id_y][list_id][qp] = (int32_t*)calloc(
            uvg_g_scaling_list_size[size_id_x][size_id_y], sizeof(int32_t));
          scaling_list->de_quant_coeff[size_id_x][size_id_y][list_id][qp] = (int32_t*)calloc(
            uvg_g_scaling_list_size[size_id_x][size_id_y], sizeof(int32_t));
          scaling_list->error_scale[size_id_x][size_id_y][list_id][qp] = (double*)calloc(
            uvg_g_scaling_list_size[size_id_x][size_id_y], sizeof(double));

        }
        scaling_list->scaling_list_coeff[size_id_x][size_id_y][list_id] = (int32_t*)calloc(
          MIN(MAX_MATRIX_COEF_NUM, uvg_g_scaling_list_size[size_id_x][size_id_y]), sizeof(int32_t));
      }
    }
  }

  //ToDo: Find out what this is about
  // alias, assign pointer to an existing array
  //for (qp = 0; qp < 6; qp++) {
  //  scaling_list->quant_coeff[3][3][qp]    = scaling_list->quant_coeff[3][1][qp];
  //  scaling_list->de_quant_coeff[3][3][qp] = scaling_list->de_quant_coeff[3][1][qp];
  //  scaling_list->error_scale[3][3][qp]    = scaling_list->error_scale[3][1][qp];
  //}

  //Initialize dc (otherwise we switch on undef in uvg_scalinglist_set)
  for (size_id_x = 0; size_id_x < SCALING_LIST_SIZE_NUM; size_id_x++) {
    for (size_id_y = 0; size_id_y < SCALING_LIST_SIZE_NUM; size_id_y++) {
      for (list_id = 0; list_id < SCALING_LIST_NUM; ++list_id) {
        scaling_list->scaling_list_dc[size_id_x][size_id_y][list_id] = 0;
      }
    }
  }

  scaling_list->enable = 0;
  scaling_list->use_default_list = 0;
}

/**
 * \brief Destroy scaling list allocated memory
 *
 */
void uvg_scalinglist_destroy(scaling_list_t* const scaling_list) {
  uint32_t size_id_y, size_id_x, list_id, qp;

  for (size_id_x = 0; size_id_x < SCALING_LIST_SIZE_NUM; size_id_x++) {
    for (size_id_y = 0; size_id_y < SCALING_LIST_SIZE_NUM; size_id_y++) {
      for (list_id = 0; list_id < SCALING_LIST_NUM; list_id++) {
        for (qp = 0; qp < 6; qp++) {
          FREE_POINTER(scaling_list->quant_coeff[size_id_x][size_id_y][list_id][qp]);
          FREE_POINTER(scaling_list->de_quant_coeff[size_id_x][size_id_y][list_id][qp]);
          FREE_POINTER(scaling_list->error_scale[size_id_x][size_id_y][list_id][qp]);

        }
        FREE_POINTER(scaling_list->scaling_list_coeff[size_id_x][size_id_y][list_id]);
      }
    }
  }
}

int uvg_scalinglist_parse(scaling_list_t* const scaling_list, FILE* fp) {
  // ToDo: fix
  return 0;
  //#define LINE_BUFSIZE 1024
  //static const char matrix_type[4][6][20] =
  //{
  //  {
  //    "INTRA4X4_LUMA",
  //    "INTRA4X4_CHROMAU",
  //    "INTRA4X4_CHROMAV",
  //    "INTER4X4_LUMA",
  //    "INTER4X4_CHROMAU",
  //    "INTER4X4_CHROMAV"
  //  },
  //  {
  //    "INTRA8X8_LUMA",
  //    "INTRA8X8_CHROMAU",
  //    "INTRA8X8_CHROMAV",
  //    "INTER8X8_LUMA",
  //    "INTER8X8_CHROMAU",
  //    "INTER8X8_CHROMAV"
  //  },
  //  {
  //    "INTRA16X16_LUMA",
  //    "INTRA16X16_CHROMAU",
  //    "INTRA16X16_CHROMAV",
  //    "INTER16X16_LUMA",
  //    "INTER16X16_CHROMAU",
  //    "INTER16X16_CHROMAV"
  //  },
  //  {
  //    "INTRA32X32_LUMA",
  //    "INTER32X32_LUMA",
  //  },
  //};
  //static const char matrix_type_dc[2][6][22] =
  //{
  //  {
  //    "INTRA16X16_LUMA_DC",
  //    "INTRA16X16_CHROMAU_DC",
  //    "INTRA16X16_CHROMAV_DC",
  //    "INTER16X16_LUMA_DC",
  //    "INTER16X16_CHROMAU_DC",
  //    "INTER16X16_CHROMAV_DC"
  //  },
  //  {
  //    "INTRA32X32_LUMA_DC",
  //    "INTER32X32_LUMA_DC",
  //  },
  //};

  //uint32_t size_id;
  //for (size_id = 0; size_id < SCALING_LIST_SIZE_NUM; size_id++) {
  //  uint32_t list_id;
  //  uint32_t size = MIN(MAX_MATRIX_COEF_NUM, (int32_t)uvg_g_scaling_list_size[size_id]);
  //  //const uint32_t * const scan = (size_id == 0) ? uvg_g_sig_last_scan[SCAN_DIAG][1] : g_sig_last_scan_32x32;

  //  for (list_id = 0; list_id < uvg_g_scaling_list_num[size_id]; list_id++) {
  //    int found;
  //    uint32_t i;
  //    int32_t data;
  //    //This IS valid (our pointer is dynamically allocated in uvg_scalinglist_init)
  //    int32_t *coeff = (int32_t*) scaling_list->scaling_list_coeff[size_id][list_id];
  //    char line[LINE_BUFSIZE + 1] = { 0 }; // +1 for null-terminator

  //    // Go back for each matrix.
  //    fseek(fp, 0, SEEK_SET);

  //    do {
  //      if (!fgets(line, LINE_BUFSIZE, fp) ||
  //          ((found = !!strstr(line, matrix_type[size_id][list_id])) == 0 && feof(fp)))
  //        return 0;
  //    } while (!found);

  //    for (i = 0; i < size;) {
  //      char *p;
  //      if (!fgets(line, LINE_BUFSIZE, fp))
  //        return 0;
  //      p = line;

  //      // Read coefficients per line.
  //      // The comma (,) character is used as a separator.
  //      // The coefficients are stored in up-right diagonal order.
  //      do {
  //        int ret = sscanf(p, "%d", &data);
  //        if (ret != 1)
  //          break;
  //        else if (data < 1 || data > 255)
  //          return 0;

  //        coeff[i++] = data;
  //        if (i == size)
  //          break;

  //        // Seek to the next newline, null-terminator or comma.
  //        while (*p != '\n' && *p != '\0' && *p != ',')
  //          ++p;
  //        if (*p == ',')
  //          ++p;
  //      } while (*p != '\n' && *p != '\0');
  //    }

  //    // Set DC value.
  //    if (size_id >= SCALING_LIST_16x16) {
  //      fseek(fp, 0, SEEK_SET);

  //      do {
  //        if (!fgets(line, LINE_BUFSIZE, fp) ||
  //            ((found = !!strstr(line, matrix_type_dc[size_id - SCALING_LIST_16x16][list_id])) == 0 && feof(fp)))
  //          return 0;
  //      } while (!found);
  //      if (1 != fscanf(fp, "%d", &data) || data < 1 || data > 255)
  //        return 0;

  //      scaling_list->scaling_list_dc[size_id][list_id] = data;
  //    } else
  //      scaling_list->scaling_list_dc[size_id][list_id] = coeff[0];
  //  }
  //}

  //scaling_list->enable = 1;
  //return 1;
  //#undef LINE_BUFSIZE
}

const int32_t* uvg_scalinglist_get_default(const uint32_t size_id, const uint32_t list_id) {
  const int32_t* list_ptr = g_quant_intra_default_8x8; // Default to "8x8" intra
  switch (size_id) {
  case SCALING_LIST_4x4:
    list_ptr = g_quant_default_4x4;
    break;
  case SCALING_LIST_8x8:
  case SCALING_LIST_16x16:
    if (list_id > 2) list_ptr = g_quant_inter_default_8x8;
    break;
  case SCALING_LIST_32x32:
    if (list_id > 0) list_ptr = g_quant_inter_default_8x8;
    break;
  }
  return list_ptr;
}


/**
 * \brief get scaling list for decoder
 *
 */
static void scalinglist_process_dec(const int32_t* const coeff, int32_t* dequantcoeff,
                                    int32_t inv_quant_scales, uint32_t height,
                                    uint32_t width, uint32_t ratio,
                                    int32_t size_num, uint32_t dc,
                                    uint8_t flat) {
  uint32_t j, i;

  // Flat scaling list
  if (flat) {
    for (j = 0; j < height * width; j++) {
      *dequantcoeff++ = inv_quant_scales << 4;
    }
  }
  else {
    for (j = 0; j < height; j++) {
      for (i = 0; i < width; i++) {
        dequantcoeff[j * width + i] = inv_quant_scales * coeff[size_num * (j / ratio) + i / ratio];
      }
    }
    if (ratio > 1) {
      dequantcoeff[0] = inv_quant_scales * dc;
    }
  }
}

/**
 * \brief get scaling list for encoder
 *
 */
void uvg_scalinglist_process_enc(const int32_t* const coeff, int32_t* quantcoeff, const int32_t quant_scales,
                                 const uint32_t height, const uint32_t width, const uint32_t ratio,
                                 const int32_t size_num, const uint32_t dc, const uint8_t flat) {
  uint32_t j, i;
  int32_t nsqth = (height < width) ? 4 : 1; //!< height ratio for NSQT
  int32_t nsqtw = (width < height) ? 4 : 1; //!< width ratio for NSQT

  // Flat scaling list
  if (flat) {
    for (j = 0; j < height * width; j++) {
      *quantcoeff++ = quant_scales >> 4;
    }
  }
  else {
    for (j = 0; j < height; j++) {
      for (i = 0; i < width; i++) {
        uint32_t coeffpos = size_num * (j * nsqth / ratio) + i * nsqtw / ratio;
        quantcoeff[j * width + i] = quant_scales / ((coeffpos > 63) ? 1 : coeff[coeffpos]);
      }
    }
    if (ratio > 1) {
      quantcoeff[0] = quant_scales / dc;
    }
  }
}


/** set error scale coefficients
 * \param list List ID
 * \param uiSize Size
 * \param uiQP Quantization parameter
 */
static void scalinglist_set_err_scale(uint8_t bitdepth, scaling_list_t* const scaling_list, uint32_t list,
                                      uint32_t size_x, uint32_t size_y, uint32_t qp) {
  const uint8_t width = g_scaling_list_size_x[size_x];
  const uint8_t height = g_scaling_list_size_x[size_y];
  int32_t transform_shift = MAX_TR_DYNAMIC_RANGE - bitdepth - ((uvg_math_floor_log2(width) + uvg_math_floor_log2(height)) >> 1); // Represents scaling through forward transform

  uint32_t i, max_num_coeff = width * height;
  const int32_t* quantcoeff = scaling_list->quant_coeff[size_x][size_y][list][qp];
  //This cast is allowed, since error_scale is a malloc'd pointer in uvg_scalinglist_init
  double* err_scale = (double*)scaling_list->error_scale[size_x][size_y][list][qp];

  const bool needsSqrt2 = ((uvg_math_floor_log2(width) + uvg_math_floor_log2(height)) & 1) == 1;
  double dTransShift = (double)transform_shift + (needsSqrt2 ? -0.5 : 0.0);
  // Compensate for scaling of bitcount in Lagrange cost function
  double scale = CTX_FRAC_ONE_BIT;
  // Compensate for scaling through forward transform
  scale = scale * pow(2.0, -2.0 * dTransShift);
  for (i = 0; i < max_num_coeff; i++) {
    err_scale[i] = scale / quantcoeff[i] / quantcoeff[i] / (1 << (2 * (bitdepth - 8)));
  }
}


/**
 * \brief set scaling lists
 *
 */
void uvg_scalinglist_set(scaling_list_t* const scaling_list, const int32_t* const coeff, uint32_t listId,
                         uint32_t size_id_x, uint32_t size_id_y, uint32_t qp) {
  const uint32_t width = g_scaling_list_size_x[size_id_x];
  const uint32_t height = g_scaling_list_size_x[size_id_y];
  const uint32_t ratio = g_scaling_list_size_x[size_id_x] / MIN(8, g_scaling_list_size_x[size_id_x]);
  const uint32_t dc = scaling_list->scaling_list_dc[size_id_x][size_id_y][listId] != 0
                        ? scaling_list->scaling_list_dc[size_id_x][size_id_y][listId]
                        : 16;
  //These cast are allowed, since these are pointer's to malloc'd area in uvg_scalinglist_init
  int32_t* quantcoeff = (int32_t*)scaling_list->quant_coeff[size_id_x][size_id_y][listId][qp];
  int32_t* dequantcoeff = (int32_t*)scaling_list->de_quant_coeff[size_id_x][size_id_y][listId][qp];

  // Encoder list
  uvg_scalinglist_process_enc(coeff, quantcoeff, uvg_g_quant_scales[qp] << 4, height, width, ratio,
                              MIN(8, g_scaling_list_size_x[size_id_x]), dc, !scaling_list->enable);
  // Decoder list
  scalinglist_process_dec(coeff, dequantcoeff, uvg_g_inv_quant_scales[qp], height, width, ratio,
                          MIN(8, g_scaling_list_size_x[size_id_x]), dc, !scaling_list->enable);


  // TODO: support NSQT
  // if(sizeId == /*SCALING_LIST_32x32*/3 || sizeId == /*SCALING_LIST_16x16*/2) { //for NSQT
  //   quantcoeff   = g_quant_coeff[listId][qp][sizeId-1][/*SCALING_LIST_VER*/1];
  //   uvg_scalinglist_process_enc(coeff,quantcoeff,g_quantScales[qp]<<4,height,width>>2,ratio,MIN(8,g_scalingListSizeX[sizeId]),/*scalingList->getScalingListDC(sizeId,listId)*/0);

  //   quantcoeff   = g_quant_coeff[listId][qp][sizeId-1][/*SCALING_LIST_HOR*/2];
  //   uvg_scalinglist_process_enc(coeff,quantcoeff,g_quantScales[qp]<<4,height>>2,width,ratio,MIN(8,g_scalingListSizeX[sizeId]),/*scalingList->getScalingListDC(sizeId,listId)*/0);
  // }
}

/**
 * \brief
 *
 */
void uvg_scalinglist_process(scaling_list_t* const scaling_list, uint8_t bitdepth) {
  uint32_t size_id_x, size_id_y, list, qp;

  for (size_id_x = 0; size_id_x < SCALING_LIST_SIZE_NUM; size_id_x++) {
    for (size_id_y = 0; size_id_y < SCALING_LIST_SIZE_NUM; size_id_y++) {
      for (list = 0; list < SCALING_LIST_NUM; list++) {
        const int32_t* const list_ptr = scaling_list->use_default_list
                                          ? uvg_scalinglist_get_default(size_id_x, list)
                                          : scaling_list->scaling_list_coeff[size_id_x][size_id_y][list];

        for (qp = 0; qp < SCALING_LIST_REM_NUM; qp++) {
          uvg_scalinglist_set(scaling_list, list_ptr, list, size_id_x, size_id_y, qp);
          scalinglist_set_err_scale(bitdepth, scaling_list, list, size_id_x, size_id_y, qp);
        }
      }
    }
  }
}
