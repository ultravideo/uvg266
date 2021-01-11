#ifndef CABAC_H_
#define CABAC_H_
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
 * \ingroup CABAC
 * \file
 * Coding bins using CABAC.
 */

#include "global.h" // IWYU pragma: keep

#include "bitstream.h"

struct encoder_state_t;

// Types
typedef struct
{
  uint16_t  state[2];
  uint8_t  rate;
} cabac_ctx_t;

typedef struct
{
  cabac_ctx_t *cur_ctx;
  uint32_t   low;
  uint32_t   range;
  uint32_t   buffered_byte;
  int32_t    num_buffered_bytes;
  int32_t    bits_left;
  int8_t     only_count;
  bitstream_t *stream;

  // CONTEXTS
  struct {
    cabac_ctx_t alf_ctb_flag_model[9];
    cabac_ctx_t alf_latest_filt;
    cabac_ctx_t alf_temporal_filt;
    cabac_ctx_t alf_ctb_alternatives[2];
    cabac_ctx_t alf_luma_coeff_delta_prediction_flag;
    cabac_ctx_t alf_cc_filter_control_flag[6];
    cabac_ctx_t sao_merge_flag_model;
    cabac_ctx_t sao_type_idx_model;
    cabac_ctx_t split_flag_model[9]; //!< \brief split flag context models
    cabac_ctx_t qt_split_flag_model[6]; //!< \brief qt split flag context models
    cabac_ctx_t intra_luma_mpm_flag_model;    //!< \brief intra mode context models
    cabac_ctx_t intra_subpart_model[2];    //!< \brief intra sub part context models
    cabac_ctx_t chroma_pred_model;
    cabac_ctx_t inter_dir[6];
    cabac_ctx_t qt_cbf_model_luma[4];
    cabac_ctx_t qt_cbf_model_cr[3];
    cabac_ctx_t qt_cbf_model_cb[2];
    cabac_ctx_t cu_qp_delta_abs[2];
    cabac_ctx_t part_size_model[4];
    cabac_ctx_t cu_sig_model_luma[3][12];
    cabac_ctx_t cu_sig_model_chroma[3][8];
    cabac_ctx_t cu_parity_flag_model_luma[21];
    cabac_ctx_t cu_parity_flag_model_chroma[11];
    cabac_ctx_t cu_gtx_flag_model_luma[2][21];
    cabac_ctx_t cu_gtx_flag_model_chroma[2][11];
    cabac_ctx_t cu_ctx_last_y_luma[20];
    cabac_ctx_t cu_ctx_last_y_chroma[3];
    cabac_ctx_t cu_ctx_last_x_luma[20];
    cabac_ctx_t cu_ctx_last_x_chroma[3];
    cabac_ctx_t cu_pred_mode_model[2];
    cabac_ctx_t cu_skip_flag_model[3];
    cabac_ctx_t cu_merge_idx_ext_model;
    cabac_ctx_t cu_merge_flag_ext_model;
    cabac_ctx_t cu_transquant_bypass;
    cabac_ctx_t cu_mvd_model[2];
    cabac_ctx_t cu_ref_pic_model[2];
    cabac_ctx_t mvp_idx_model;
    cabac_ctx_t cu_qt_root_cbf_model;
    cabac_ctx_t sig_coeff_group_model[4];
    cabac_ctx_t luma_planar_model[2];
    cabac_ctx_t multi_ref_line[2];
    cabac_ctx_t bdpcm_mode[4];
    cabac_ctx_t joint_bc_br[3];
  } ctx;
} cabac_data_t;


// Globals
extern const uint8_t kvz_g_auc_renorm_table[32];


// Functions
void kvz_cabac_start(cabac_data_t *data);
void kvz_cabac_encode_bin(cabac_data_t *data, uint32_t bin_value);
void kvz_cabac_encode_bin_ep(cabac_data_t *data, uint32_t bin_value);
void kvz_cabac_encode_trunc_bin(cabac_data_t *data, uint32_t bin_value, uint32_t max_value);
void kvz_cabac_encode_bins_ep(cabac_data_t *data, uint32_t bin_values, int num_bins);
void kvz_cabac_encode_bin_trm(cabac_data_t *data, uint8_t bin_value);
void kvz_cabac_write(cabac_data_t *data);
void kvz_cabac_finish(cabac_data_t *data);
void kvz_cabac_write_coeff_remain(cabac_data_t *cabac, uint32_t symbol,
                              uint32_t r_param, const unsigned int cutoff);
void kvz_cabac_write_ep_ex_golomb(struct encoder_state_t * const state, cabac_data_t *data,
                uint32_t symbol, uint32_t count);
void kvz_cabac_write_unary_max_symbol(cabac_data_t *data, cabac_ctx_t *ctx,
                                  uint32_t symbol, int32_t offset,
                                  uint32_t max_symbol);
void kvz_cabac_write_unary_max_symbol_ep(cabac_data_t *data, unsigned int symbol, unsigned int max_symbol);

#define CTX_PROB_BITS 15
#define CTX_PROB_BITS_0 10
#define CTX_PROB_BITS_1 14
#define CTX_MASK_0 (~(~0u << CTX_PROB_BITS_0) << (CTX_PROB_BITS - CTX_PROB_BITS_0))
#define CTX_MASK_1 (~(~0u << CTX_PROB_BITS_1) << (CTX_PROB_BITS - CTX_PROB_BITS_1))

// Macros
#define CTX_GET_STATE(ctx) ( (ctx)->state[0]+(ctx)->state[1] )
#define CTX_STATE(ctx) ( CTX_GET_STATE(ctx)>>8 )
#define CTX_SET_STATE(ctx, state) {\
  (ctx)->state[0]=(state >> 1) & (int)CTX_MASK_0;\
  (ctx)->state[1]=(state >> 1) & (int)CTX_MASK_1;\
}
#define CTX_MPS(ctx) (CTX_STATE(ctx)>>7)
#define CTX_LPS(ctx,range) ((uint8_t)( ((((CTX_STATE(ctx)&0x80) ? (CTX_STATE(ctx)^0xff) : (CTX_STATE(ctx))) >>2)*(range>>5)>>1)+4  ))
#define CTX_UPDATE(ctx,bin) { \
  int rate0 = (ctx)->rate >> 4;\
  int rate1 = (ctx)->rate & 15;\
\
  (ctx)->state[0] -= ((ctx)->state[0] >> rate0) & (int)CTX_MASK_0;\
  (ctx)->state[1] -= ((ctx)->state[1] >> rate1) & (int)CTX_MASK_1;\
  if (bin) {\
    (ctx)->state[0] += (0x7fffu >> rate0) & (int)CTX_MASK_0;\
    (ctx)->state[1] += (0x7fffu >> rate1) & (int)CTX_MASK_1;\
  }\
}

#define CTX_SET_LOG2_WIN(ctx, size) { \
  int rate0 = 2 + ((size >> 2) & 3); \
  int rate1 = 3 + rate0 + (size & 3);\
 (ctx)->rate = 16 * rate0 + rate1;\
}

#ifdef KVZ_DEBUG_PRINT_CABAC
extern uint32_t kvz_cabac_bins_count;
extern bool kvz_cabac_bins_verbose;
#define CABAC_BIN(data, value, name) { \
    uint32_t prev_state = CTX_STATE(data->cur_ctx); \
    if(kvz_cabac_bins_verbose && !data->only_count) {printf("%d %d  [%d:%d]  %s = %u, range = %u LPS = %u state = %u -> ", \
           kvz_cabac_bins_count++, (data)->range, (data)->range-CTX_LPS(data->cur_ctx,(data)->range), CTX_LPS(data->cur_ctx,(data)->range), (name), (uint32_t)(value), (data)->range, CTX_LPS(data->cur_ctx,(data)->range), prev_state); }\
    kvz_cabac_encode_bin((data), (value)); \
    if(kvz_cabac_bins_verbose && !data->only_count) printf("%u\n", CTX_STATE(data->cur_ctx)); }
    

  #define CABAC_BINS_EP(data, value, bins, name) { \
    uint32_t prev_state = CTX_STATE(data->cur_ctx); \
    kvz_cabac_encode_bins_ep((data), (value), (bins)); \
    if(kvz_cabac_bins_verbose && !data->only_count) { printf("%d %s = %u(%u bins), state = %u -> %u\n", \
           kvz_cabac_bins_count, (name), (uint32_t)(value), (bins), prev_state, CTX_STATE(data->cur_ctx));  kvz_cabac_bins_count+=bins;}}

  #define CABAC_BIN_EP(data, value, name) { \
    uint32_t prev_state = CTX_STATE(data->cur_ctx); \
    kvz_cabac_encode_bin_ep((data), (value)); \
    if(kvz_cabac_bins_verbose && !data->only_count) {printf("%d %s = %u, state = %u -> %u\n", \
           kvz_cabac_bins_count++, (name), (uint32_t)(value), prev_state, CTX_STATE(data->cur_ctx)); }}
#else
  #define CABAC_BIN(data, value, name) \
    kvz_cabac_encode_bin((data), (value));
  #define CABAC_BINS_EP(data, value, bins, name) \
    kvz_cabac_encode_bins_ep((data), (value), (bins));
  #define CABAC_BIN_EP(data, value, name) \
    kvz_cabac_encode_bin_ep((data), (value));
#endif

#endif
