#ifndef CABAC_H_
#define CABAC_H_
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
  int8_t     only_count : 4;
  int8_t     update : 4;
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
    cabac_ctx_t lfnst_idx_model[3];
    cabac_ctx_t mts_idx_model[4];
    cabac_ctx_t split_flag_model[9]; //!< \brief split flag context models
    cabac_ctx_t qt_split_flag_model[6]; //!< \brief qt split flag context models
    cabac_ctx_t mtt_vertical_model[5]; 
    cabac_ctx_t mtt_binary_model[4]; 
    cabac_ctx_t intra_luma_mpm_flag_model;    //!< \brief intra mode context models
    cabac_ctx_t intra_subpart_model[2];    //!< \brief intra sub part context models
    cabac_ctx_t chroma_pred_model;
    cabac_ctx_t inter_dir[6];
    cabac_ctx_t imv_flag[5];
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
    cabac_ctx_t mip_flag[4];
    cabac_ctx_t bdpcm_mode[4];
    cabac_ctx_t joint_cb_cr[3];
    cabac_ctx_t transform_skip_model_luma;
    cabac_ctx_t transform_skip_model_chroma;
    cabac_ctx_t transform_skip_sig_coeff_group[3];
    cabac_ctx_t transform_skip_sig[3];
    cabac_ctx_t transform_skip_res_sign[6];
    cabac_ctx_t transform_skip_gt1[4];
    cabac_ctx_t transform_skip_par;
    cabac_ctx_t transform_skip_gt2[5];
    cabac_ctx_t cclm_flag;
    cabac_ctx_t cclm_model;
    cabac_ctx_t ibc_flag[3];

  } ctx;
} cabac_data_t;


// Globals
extern const uint8_t uvg_g_auc_renorm_table[32];


// Functions
void uvg_cabac_start(cabac_data_t *const data);
void uvg_cabac_encode_bin(cabac_data_t *const data, const uint32_t bin_value);
void uvg_cabac_encode_bin_ep(cabac_data_t *const data, const uint32_t bin_value);
void uvg_cabac_encode_trunc_bin(cabac_data_t *const data, const uint32_t bin_value, const uint32_t max_value, double* bits_out);
void uvg_cabac_encode_bins_ep(cabac_data_t *const data, uint32_t bin_values, int num_bins);
void uvg_cabac_encode_bin_trm(cabac_data_t *const data, const uint8_t bin_value);
void uvg_cabac_write(cabac_data_t *const data);
void uvg_cabac_finish(cabac_data_t *const data);
int uvg_cabac_write_coeff_remain(cabac_data_t *const cabac, const uint32_t symbol,
                              const uint32_t r_param, const unsigned int cutoff);
uint32_t uvg_cabac_write_ep_ex_golomb(struct encoder_state_t * const state, cabac_data_t *const data,
                uint32_t symbol, uint32_t count);
void uvg_cabac_write_unary_max_symbol(cabac_data_t *const data, cabac_ctx_t *const ctx,
                                      uint32_t symbol, const int32_t offset,
                                      const uint32_t max_symbol, double* bits_out);
void uvg_cabac_write_unary_max_symbol_ep(cabac_data_t *const data, unsigned int symbol, const unsigned int max_symbol);

#define CTX_PROB_BITS 15
#define CTX_PROB_BITS_0 10
#define CTX_PROB_BITS_1 14
#define CTX_MASK_0 (~(~0u << CTX_PROB_BITS_0) << (CTX_PROB_BITS - CTX_PROB_BITS_0))
#define CTX_MASK_1 (~(~0u << CTX_PROB_BITS_1) << (CTX_PROB_BITS - CTX_PROB_BITS_1))

// Floating point fractional bits, derived from kvz_entropy_bits
extern const float uvg_f_entropy_bits[512];
#define CTX_ENTROPY_FBITS(ctx, val) uvg_f_entropy_bits[(CTX_STATE(ctx)<<1) ^ (val)]

#define CABAC_FBITS_UPDATE(cabac, ctx, val, bits, name) do { \
  if((cabac)->only_count) (bits) += uvg_f_entropy_bits[(CTX_STATE(ctx)<<1) ^ (val)]; \
  if((cabac)->update) {\
    (cabac)->cur_ctx = ctx;\
    CABAC_BIN((cabac), (val), (name));\
  } \
} while(0)

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

#ifdef UVG_DEBUG_PRINT_CABAC
extern uint32_t uvg_cabac_bins_count;
extern bool uvg_cabac_bins_verbose;
#define CABAC_BIN(data, value, name) { \
    uint32_t prev_state = CTX_STATE(data->cur_ctx); \
    if(uvg_cabac_bins_verbose && !(data)->only_count) {printf("%d %d  [%d:%d]  %s = %u, range = %u LPS = %u state = %u -> ", \
           uvg_cabac_bins_count++, (data)->range, (data)->range-CTX_LPS((data)->cur_ctx,(data)->range), CTX_LPS((data)->cur_ctx,(data)->range), (name), (uint32_t)(value), (data)->range, CTX_LPS((data)->cur_ctx,(data)->range), prev_state); }\
    uvg_cabac_encode_bin((data), (value)); \
    if(uvg_cabac_bins_verbose && !(data)->only_count) printf("%u\n", CTX_STATE((data)->cur_ctx)); }
    

  #define CABAC_BINS_EP(data, value, bins, name) { \
    uint32_t prev_state = (!(data)->only_count) ? CTX_STATE(data->cur_ctx) : 0; \
    uvg_cabac_encode_bins_ep((data), (value), (bins)); \
    if(uvg_cabac_bins_verbose && !data->only_count) { printf("%d %s = %u(%u bins), state = %u -> %u\n", \
           uvg_cabac_bins_count, (name), (uint32_t)(value), (bins), prev_state, CTX_STATE((data)->cur_ctx));  uvg_cabac_bins_count+=(bins);}}

  #define CABAC_BIN_EP(data, value, name) { \
    uint32_t prev_state = (!(data)->only_count) ? CTX_STATE((data)->cur_ctx) : 0;; \
    uvg_cabac_encode_bin_ep((data), (value)); \
    if(uvg_cabac_bins_verbose && !(data)->only_count) {printf("%d %s = %u, state = %u -> %u\n", \
           uvg_cabac_bins_count++, (name), (uint32_t)(value), prev_state, CTX_STATE((data)->cur_ctx)); }}
#else
  #define CABAC_BIN(data, value, name) \
    uvg_cabac_encode_bin((data), (value));
  #define CABAC_BINS_EP(data, value, bins, name) \
    uvg_cabac_encode_bins_ep((data), (value), (bins));
  #define CABAC_BIN_EP(data, value, name) \
    uvg_cabac_encode_bin_ep((data), (value));
#endif

#endif
