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

#include "context.h"

#include "tables.h"



static const uint8_t  INIT_SPLIT_FLAG[4][9] = {
  { 122, 124, 141, 108, 125, 156, 138, 126, 143, },
  { 93, 139, 171, 124, 125, 141, 139, 141, 158, },
  { 138, 154, 172, 124, 140, 142, 154, 127, 175, },
  { 9, 13, 8, 8, 13, 12, 5, 10, 12, },
};

static const uint8_t  INIT_QT_SPLIT_FLAG[4][6] = {
  { 138, 140, 142, 136, 138, 140, },
  { 139, 126, 142, 107, 138, 125, },
  { 139, 125, 127, 136, 153, 126, },
  { 0, 8, 8, 12, 12, 8, },
};

static const uint8_t INIT_SKIP_FLAG[4][3] = {
  { 197, 214, 216, },
  { 197, 198, 185, },
  { 40, 138, 154, },
  { 5, 8, 8, },
};

static const uint8_t INIT_MERGE_FLAG_EXT[4][1] = {
  { 111, },
  { 111, },
  { 153, },
  { 5, },
};


static const uint8_t INIT_MERGE_IDX_EXT[4][1] = {
  { 138, },
  { 154, },
  { 153, },
  { 8, },
};

static const uint8_t INIT_PART_SIZE[4][4] = {
  {  CNU, CNU, CNU, CNU,},
  {  CNU, CNU, CNU, CNU,},
  {  CNU, CNU, CNU, CNU,},
  { DWS, DWS, DWS, DWS, }
};

static const uint8_t INIT_PRED_MODE[4][2] = {
  { 192, 168, },
  { 165, 139, },
  { CNU, CNU, },
  { 5, 2, },
};

static const uint8_t INIT_INTRA_PRED_MODE[4] = {
  154, 154, 170, 6
};

static const uint8_t INIT_CHROMA_PRED_MODE[4][3] = {
  {  137, 139, 140,},
  {  138, 139, 169,},
  {  154, 139, 154,},
  {    5,   8,   9,},
};

static const uint8_t INIT_CU_QP_DELTA_ABS[4][3] = {
  {  154, 154, 154,},
  {  154, 154, 154,},
  {  154, 154, 154,},
  { DWS, DWS, DWS, }
};

static const uint8_t INIT_INTER_DIR[4][5] = {
  { 111, 125, 110, 94, 192, },
  { 126, 111, 110, 94, 208, },
  { CNU, CNU, CNU, CNU, CNU, },
  { 0, 0, 4, 5, 0, },
};

static const uint8_t INIT_REF_PIC[4][2] = {
  { 125, 139, },
  { 138, 168, },
  { CNU, CNU, },
  { 4, 5, },
};

static const uint8_t INIT_MVD[4][2] = {
  { 169, 183, },
  { 155, 154, },
  { 141, 156, },
  { 9, 5, },
};

static const uint8_t INIT_QT_ROOT_CBF[4][1] = {
  { 109, },
  { 95, },
  { 110, },
  { 4, },
};

static const uint8_t INIT_QT_CBF[4][15] = {
  { 141, 127, 139, 140, 163, 154, CNU, CNU, CNU, 161, 154,},
  { 142, 127, 139, 140, 164, 154, CNU, CNU, CNU, 192, 154,},
  { CNU, 111, 124, 111, 109, CNU, CNU, CNU, CNU, 151, 155,},
  { 1, 5, 9, 8, 5, 8, DWS, DWS, DWS, 5, 5, },
};

static const uint8_t INIT_SIG_COEFF_GROUP[4][8] = {
    { 105, 155, 91, 155, CNU, CNU, CNU, CNU, },
    { 106, 156, 90, 141, CNU, CNU, CNU, CNU, },
    { 107, 158, 76, 127, CNU, CNU, CNU, CNU, },
    { 8, 5, 5, 8, DWS, DWS,DWS, DWS,},
};

static const uint8_t INIT_SIG_FLAG[6][4][18] = {
  {
    { 88, 166, 152, 182, 168, 154, 0, 167, 182, 168, 183, 155, 193, 213, 183, 183, 169, 185, },
    { 132, 152, 167, 168, 183, 140, 177, 182, 168, 154, 169, 155, 180, 213, 183, 169, 184, 156, },
    { 89, 138, 153, 139, 154, 140, 134, 139, 139, 140, 140, 141, 137, 170, 169, 170, 141, 157, },
    { 12, 9, 9, 9, 9, 10, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 9, },
  },
  {
    { 72, 167, 153, 168, 154, 155, 180, 199, 183, 199, 199, 186, },
    { 133, 138, 153, 139, 154, 140, 181, 229, 169, 229, 170, 157, },
    { 43, 153, 168, 169, 154, 155, 152, 215, 155, 201, 171, 143, },
    { 9, 9, 12, 9, 13, 13, 5, 5, 8, 8, 8, 9, },
  },
  {
    { 152, 156, 201, 186, 186, 187, 182, 248, 188, 232, 188, 205, 182, 223, 223, 223, 223, 223, },
    { 123, 142, 157, 172, 172, 218, 138, 249, 248, 248, 219, 223, 139, 223, 223, 223, 223, 223, },
    { 93, 142, 157, 143, 188, 175, 138, 238, 205, 238, 253, 237, 139, 223, 223, 223, 223, 253, },
    { 9, 12, 8, 8, 8, 8, 8, 8, 8, 8, 8, 5, 8, 0, 0, 0, 0, 0, },
  },
  {
    { 182, 171, 143, 158, 172, 189, 183, 223, 223, 223, 223, 223, },
    { 168, 156, 173, 216, 172, 219, 169, 223, 223, 223, 223, 223, },
    { 152, 173, 157, 187, 204, 253, 170, 223, 223, 223, 223, 223, },
    { 8, 9, 12, 8, 8, 8, 4, 0, 2, 2, 2, 2, },
  },
  {
    { 123, 173, 223, 191, 232, 251, 212, 223, 223, 236, 206, 223, 192, 223, 223, 223, 223, 223, },
    { 123, 175, 223, 175, 218, 223, 138, 223, 223, 223, 222, 223, 196, 223, 223, 223, 223, 223, },
    { 107, 174, 223, 238, 251, 223, 63, 223, 223, 238, 223, 238, 12, 223, 223, 223, 223, 223, },
    { 8, 8, 4, 8, 8, 8, 8, 0, 0, 4, 8, 5, 4, 2, 2, 2, 2, 1, },
  },
  {
    { 167, 201, 223, 248, 219, 223, 181, 223, 223, 223, 223, 223, },
    { 167, 171, 223, 175, 248, 223, 152, 223, 223, 223, 223, 223, },
    { 166, 234, 223, 236, 248, 223, 108, 223, 223, 223, 223, 223, },
    { 8, 8, 5, 8, 8, 8, 5, 1, 2, 2, 2, 2, },
  }
};

static const uint8_t INIT_PARITY_FLAG[2][4][21] =
{
  {
    { 121, 105, 136, 152, 138, 183, 90, 122, 167, 153, 168, 135, 152, 153, 168, 139, 151, 153, 139, 168, 154, },
    { 121, 119, 136, 137, 138, 153, 104, 122, 138, 153, 139, 106, 138, 153, 168, 139, 137, 153, 168, 139, 139, },
    { 121, 135, 137, 152, 138, 153, 91, 137, 138, 153, 139, 151, 138, 153, 139, 139, 138, 168, 139, 154, 139, },
    { 8, 9, 12, 13, 13, 13, 10, 13, 13, 13, 13, 13, 13, 13, 13, 13, 10, 13, 13, 13, 13, },
  },
  {
    { 151, 120, 152, 138, 153, 153, 136, 168, 154, 168, 154, },
    { 135, 120, 137, 138, 138, 153, 136, 153, 168, 139, 154, },
    { 136, 135, 152, 153, 138, 153, 136, 168, 154, 139, 154, },
    { 8, 10, 12, 12, 13, 13, 10, 10, 13, 13, 13, },
  }
};

static const uint8_t INIT_GTX_FLAG[4][4][21] =
{
  {
    { 31, 73, 118, 75, 152, 109, 42, 44, 105, 107, 109, 0, 119, 136, 152, 124, 118, 136, 138, 153, 140, },
    { 14, 116, 86, 119, 106, 152, 0, 72, 120, 151, 138, 116, 90, 107, 152, 153, 104, 107, 123, 153, 154, },
    { 90, 72, 119, 135, 137, 138, 43, 60, 106, 137, 109, 58, 106, 108, 109, 124, 121, 138, 139, 154, 155, },
    { 4, 1, 8, 8, 4, 2, 5, 9, 9, 8, 9, 9, 9, 9, 8, 9, 9, 8, 9, 8, 8, },
  },
  {
    { 119, 101, 134, 151, 107, 123, 118, 122, 124, 140, 155, },
    { 117, 0, 90, 106, 92, 93, 147, 136, 138, 154, 140, },
    { 194, 40, 120, 122, 122, 138, 103, 121, 153, 154, 155, },
    { 2, 5, 8, 8, 8, 6, 6, 8, 8, 8, 7, },
  },
  {
    { 43, 177, 181, 168, 154, 170, 133, 167, 139, 154, 155, 164, 153, 154, 169, 155, 181, 183, 169, 185, 186, },
    { 101, 133, 137, 153, 139, 140, 134, 138, 139, 154, 155, 136, 153, 154, 140, 170, 138, 154, 155, 170, 186, },
    { 134, 120, 123, 153, 139, 140, 92, 124, 154, 125, 111, 138, 154, 140, 155, 141, 154, 140, 185, 171, 143, },
    { 8, 5, 9, 9, 12, 9, 9, 10, 13, 12, 10, 9, 10, 10, 10, 10, 8, 9, 8, 8, 10, },
  },
  {
    { 0, 178, 153, 154, 140, 140, 196, 170, 186, 157, 188, },
    { 0, 135, 153, 139, 125, 140, 182, 155, 156, 142, 159, },
    { 163, 136, 153, 154, 125, 140, 183, 170, 201, 187, 174, },
    { 6, 9, 10, 12, 12, 10, 5, 9, 8, 8, 9, },
  }
};

static const uint8_t INIT_LAST_X[4][29] = {
  { 111, 111, 110, 111, 111, 139, 111, 126, 111, 139, 126, 126, 111, 111, 169, 154, 111, 110, 110, 139, CNU, CNU, CNU, CNU, CNU, 122, 124, 63, CNU,},
  { 125, 110, 109, 125, 125, 123, 111, 111,  95, 123, 126, 111, 110,  95, 169, 154, 140, 139, 139, 138, CNU, CNU, CNU, CNU, CNU, 138, 123, 92, CNU,},
  { 125, 140, 124, 111, 111, 109, 111, 126, 125, 123, 111, 141, 111, 125,  79, 155, 142, 170, 140, 183, CNU, CNU, CNU, CNU, CNU, 138, 108, 47, CNU,},
  { 8, 5, 5, 5, 4, 4, 5, 4, 4, 0, 5, 1, 0, 0, 0, 1, 1, 0, 0, 0, DWS, DWS, DWS, DWS, DWS, 2, 1, 1, DWS, },
};

static const uint8_t INIT_LAST_Y[4][29] = {
  { 125, 125, 139, 125, 111, 139, 111, 111, 110, 110, 140, 126, 125, 125, 140, 139, 111, 110, 124, 181, CNU, CNU, CNU, CNU, CNU, 122, 124, 123, CNU,},
    { 95, 95, 109, 110, 110, 108, 125, 111, 124, 123, 140, 111, 110, 124, 139, 125, 126, 110, 124, 182, CNU, CNU, CNU, CNU, CNU, 108, 123, 121, CNU,},
    { 110, 110, 109, 125, 111, 123, 111, 126, 95, 108, 111, 127, 111, 95, 78, 169, 157, 141, 125, 138, CNU, CNU, CNU, CNU, CNU, 123, 123, 91, CNU,},
    { 8, 5, 8, 5, 5, 4, 5, 5, 4, 0, 5, 5, 1, 0, 0, 1, 4, 1, 0, 0, DWS, DWS, DWS, DWS, DWS, 2, 2, 2, DWS,},
};

static const uint8_t INIT_MVP_IDX[4][1] = {
  { 153, },
  { 168, },
  { 168, },
  { 10, },
};

static const uint8_t INIT_SAO_MERGE_FLAG[4] = { 47, 244, 199, 0 };

static const uint8_t INIT_SAO_TYPE_IDX[4] = { 47, 95, 95, 0 };

static const uint8_t INIT_CU_TRANSQUANT_BYPASS[4][1] = {
  { 154, },
  { 154, },
  { 154, },
  { DWS, }
};

static const uint8_t INIT_INTRA_SUBPART_MODE[4][2] = {
  { 152, 154, },
  { 166, 154, },
  { 152, 154, },
  { 8, 5, },
};

static const uint16_t g_inistateToCount[128] = {
  614,   647,   681,   718,   756,   797,   839,   884,   932,   982,   1034,  1089,  1148,  1209,  1274,  1342,
  1414,  1490,  1569,  1653,  1742,  1835,  1933,  2037,  2146,  2261,  2382,  2509,  2643,  2785,  2934,  3091,
  3256,  3430,  3614,  3807,  4011,  4225,  4452,  4690,  4941,  5205,  5483,  5777,  6086,  6412,  6755,  7116,
  7497,  7898,  8320,  8766,  9235,  9729,  10249, 10798, 11375, 11984, 12625, 13300, 14012, 14762, 15551, 16384,
  16384, 17216, 18005, 18755, 19467, 20142, 20783, 21392, 21969, 22518, 23038, 23532, 24001, 24447, 24869, 25270,
  25651, 26012, 26355, 26681, 26990, 27284, 27562, 27826, 28077, 28315, 28542, 28756, 28960, 29153, 29337, 29511,
  29676, 29833, 29982, 30124, 30258, 30385, 30506, 30621, 30730, 30834, 30932, 31025, 31114, 31198, 31277, 31353,
  31425, 31493, 31558, 31619, 31678, 31733, 31785, 31835, 31883, 31928, 31970, 32011, 32049, 32086, 32120, 32153
};


/**
 * \brief Initialize struct cabac_ctx.
 */
void kvz_ctx_init(cabac_ctx_t *ctx, uint32_t qp, uint32_t init_value, uint8_t rate)
{
  int slope = (init_value >> 4) * 5 - 45;
  int offset = ((init_value & 15) << 3) - 16;
  int init_state = MIN(MAX(1, ((slope * (int)qp) >> 4) + offset), 126);
    
  const int p1 = g_inistateToCount[init_state < 0 ? 0 : init_state > 127 ? 127 : init_state];

  ctx->state[0] = p1 & CTX_MASK_0;
  ctx->state[1] = p1 & CTX_MASK_1;

  CTX_SET_LOG2_WIN(ctx, rate); 

}

/**
 * \brief Initialize cabac context to be used for coding
 * \param encoder encoder control struct
 * \param slice type of slice we are coding (P/B/I)
 */

void kvz_init_contexts(encoder_state_t *state, int8_t QP, int8_t slice)
{
  cabac_data_t * const cabac = &state->cabac;
  uint16_t i, ii;

  // Initialize contexts
  kvz_ctx_init(&cabac->ctx.sao_merge_flag_model, QP, INIT_SAO_MERGE_FLAG[slice], INIT_SAO_MERGE_FLAG[3]);
  kvz_ctx_init(&cabac->ctx.sao_type_idx_model, QP, INIT_SAO_TYPE_IDX[slice], INIT_SAO_TYPE_IDX[3]);

  kvz_ctx_init(&cabac->ctx.cu_merge_flag_ext_model, QP, INIT_MERGE_FLAG_EXT[slice][0], INIT_MERGE_FLAG_EXT[3][0]);
  kvz_ctx_init(&cabac->ctx.cu_merge_idx_ext_model, QP, INIT_MERGE_IDX_EXT[slice][0], INIT_MERGE_IDX_EXT[3][0]);
  kvz_ctx_init(&cabac->ctx.cu_pred_mode_model, QP, INIT_PRED_MODE[slice][0], INIT_PRED_MODE[3][0]);
  kvz_ctx_init(&cabac->ctx.cu_transquant_bypass, QP, INIT_CU_TRANSQUANT_BYPASS[slice][0], INIT_CU_TRANSQUANT_BYPASS[3][0]);
  
  kvz_ctx_init(&cabac->ctx.intra_mode_model, QP, INIT_INTRA_PRED_MODE[slice], INIT_INTRA_PRED_MODE[3]);


  kvz_ctx_init(&cabac->ctx.intra_subpart_model[0], QP, INIT_INTRA_SUBPART_MODE[slice][0], INIT_INTRA_SUBPART_MODE[3][0]);
  kvz_ctx_init(&cabac->ctx.intra_subpart_model[1], QP, INIT_INTRA_SUBPART_MODE[slice][1], INIT_INTRA_SUBPART_MODE[3][1]);


  for (i = 0; i < 3; i++) {
    kvz_ctx_init(&cabac->ctx.cu_skip_flag_model[i], QP, INIT_SKIP_FLAG[slice][i], INIT_SKIP_FLAG[3][i]);
    kvz_ctx_init(&cabac->ctx.chroma_pred_model[i], QP, INIT_CHROMA_PRED_MODE[slice][i], INIT_CHROMA_PRED_MODE[3][i]);
  }
  
  for (i = 0; i < 6; i++) {
    kvz_ctx_init(&cabac->ctx.qt_split_flag_model[i], QP, INIT_QT_SPLIT_FLAG[slice][i], INIT_QT_SPLIT_FLAG[3][i]);
  }

  for (i = 0; i < 8; i++) {
    kvz_ctx_init(&cabac->ctx.sig_coeff_group_model[i], QP, INIT_SIG_COEFF_GROUP[slice][i], INIT_SIG_COEFF_GROUP[3][i]);
  }
  

  for (i = 0; i < 9; i++) {
    kvz_ctx_init(&cabac->ctx.split_flag_model[i], QP, INIT_SPLIT_FLAG[slice][i], INIT_SPLIT_FLAG[3][i]);
  }
  

  //TODO: ignore P/B contexts on intra frame
  kvz_ctx_init(&cabac->ctx.cu_qt_root_cbf_model, QP, INIT_QT_ROOT_CBF[slice][0], INIT_QT_ROOT_CBF[3][0]);

  for (i = 0; i < 2; i++) {
    kvz_ctx_init(&cabac->ctx.qt_cbf_model_cr[i], QP, INIT_QT_CBF[slice][i + 9], INIT_QT_CBF[3][i + 9]);
    kvz_ctx_init(&cabac->ctx.cu_mvd_model[i], QP, INIT_MVD[slice][i], INIT_MVD[3][i]);
    kvz_ctx_init(&cabac->ctx.cu_ref_pic_model[i], QP, INIT_REF_PIC[slice][i], INIT_REF_PIC[3][i]);
    kvz_ctx_init(&cabac->ctx.mvp_idx_model[i], QP, INIT_MVP_IDX[slice][i], INIT_MVP_IDX[3][i]);
    
  }

  for (i = 0; i < 3; i++) {
    kvz_ctx_init(&cabac->ctx.cu_qp_delta_abs[i], QP, INIT_CU_QP_DELTA_ABS[slice][i], INIT_CU_QP_DELTA_ABS[3][i]);
  }

  for (i = 0; i < 4; i++) {
    kvz_ctx_init(&cabac->ctx.qt_cbf_model_luma[i], QP, INIT_QT_CBF[slice][i], INIT_QT_CBF[3][i]);
    kvz_ctx_init(&cabac->ctx.part_size_model[i], QP, INIT_PART_SIZE[slice][i], INIT_PART_SIZE[3][i]);

    kvz_ctx_init(&cabac->ctx.cu_ctx_last_y_chroma[i], QP, INIT_LAST_X[slice][i + 25], INIT_LAST_X[3][i + 25]);
    kvz_ctx_init(&cabac->ctx.cu_ctx_last_x_chroma[i], QP, INIT_LAST_Y[slice][i + 25], INIT_LAST_Y[3][i + 25]);
  }

  for (i = 0; i < 5; i++) {
    kvz_ctx_init(&cabac->ctx.qt_cbf_model_cb[i], QP, INIT_QT_CBF[slice][i + 4], INIT_QT_CBF[3][i + 4]);

    kvz_ctx_init(&cabac->ctx.inter_dir[i], QP, INIT_INTER_DIR[slice][i], INIT_INTER_DIR[3][i]);
  }

  for (i = 0; i < 25; i++) {
    kvz_ctx_init(&cabac->ctx.cu_ctx_last_y_luma[i], QP, INIT_LAST_Y[slice][i], INIT_LAST_Y[3][i]);
    kvz_ctx_init(&cabac->ctx.cu_ctx_last_x_luma[i], QP, INIT_LAST_X[slice][i], INIT_LAST_X[3][i]);
  }

  for (ii = 0; ii < 3; ii++) {
    for (i = 0; i < 18; i++) {
      kvz_ctx_init(&cabac->ctx.cu_sig_model_luma[ii][i], QP, INIT_SIG_FLAG[ii*2][slice][i], INIT_SIG_FLAG[ii * 2][3][i]);
      if (i < 12) kvz_ctx_init(&cabac->ctx.cu_sig_model_chroma[ii][i], QP, INIT_SIG_FLAG[ii*2+1][slice][i], INIT_SIG_FLAG[ii * 2 + 1][3][i]);
    }
  }

  for (i = 0; i < 21; i++) {
    kvz_ctx_init(&cabac->ctx.cu_parity_flag_model_luma[i], QP, INIT_PARITY_FLAG[0][slice][i], INIT_PARITY_FLAG[0][3][i]);
    if (i < 11) kvz_ctx_init(&cabac->ctx.cu_parity_flag_model_chroma[i], QP, INIT_PARITY_FLAG[1][slice][i], INIT_PARITY_FLAG[1][3][i]);
  }
  for (ii = 0; ii < 2; ii++) {
    for (i = 0; i < 21; i++) {
      kvz_ctx_init(&cabac->ctx.cu_gtx_flag_model_luma[ii][i], QP, INIT_GTX_FLAG[ii * 2][slice][i], INIT_GTX_FLAG[ii * 2][3][i]);
      if (i < 11) kvz_ctx_init(&cabac->ctx.cu_gtx_flag_model_chroma[ii][i], QP, INIT_GTX_FLAG[ii * 2 + 1][slice][i], INIT_GTX_FLAG[ii * 2 + 1][3][i]);
    }
  }
}

void kvz_context_copy(encoder_state_t * const target_state, const encoder_state_t * const source_state) {
  cabac_data_t * const target_cabac = &target_state->cabac;
  const cabac_data_t * const source_cabac = &source_state->cabac;
  
  if (target_cabac == source_cabac) return;

  target_cabac->ctx = source_cabac->ctx;
}


uint32_t kvz_context_get_sig_coeff_group( uint32_t *sig_coeff_group_flag,
                                      uint32_t pos_x,
                                      uint32_t pos_y,
                                      int32_t width)
{
  uint32_t uiRight = 0;
  uint32_t uiLower = 0;
  width >>= 2;
  if (pos_x < (uint32_t)width - 1) uiRight = (sig_coeff_group_flag[pos_y * width + pos_x + 1] != 0);
  if (pos_y < (uint32_t)width - 1) uiLower = (sig_coeff_group_flag[(pos_y  + 1 ) * width + pos_x] != 0);

  return uiRight || uiLower;
}


/**
 * \brief Pattern decision for context derivation process of significant_coeff_flag
 * \param sig_coeff_group_flag pointer to prior coded significant coeff group
 * \param pos_x column of current coefficient group
 * \param pos_y row of current coefficient group
 * \param width width of the block
 * \returns pattern for current coefficient group
 */

int32_t kvz_context_calc_pattern_sig_ctx(const uint32_t *sig_coeff_group_flag, uint32_t pos_x, uint32_t pos_y, int32_t width)
{
  uint32_t sigRight = 0;
  uint32_t sigLower = 0;

  if (width == 4) return -1;

  width >>= 2;
  if (pos_x < (uint32_t)width - 1) sigRight = (sig_coeff_group_flag[pos_y * width + pos_x + 1] != 0);
  if (pos_y < (uint32_t)width - 1) sigLower = (sig_coeff_group_flag[(pos_y  + 1 ) * width + pos_x] != 0);

  return sigRight + (sigLower<<1);
}


/**
 * \brief Context derivation process of coeff_abs_significant_flag
 * \param pattern_sig_ctx pattern for current coefficient group
 * \param scan_idx pixel scan type in use
 * \param pos_x column of current scan position
 * \param pos_y row of current scan position
 * \param block_type log2 value of block size if square block, or 4 otherwise
 * \param width width of the block
 * \param texture_type texture type (TEXT_LUMA...)
 * \returns ctx_inc for current scan position
 */

int32_t kvz_context_get_sig_ctx_inc(int32_t pattern_sig_ctx, uint32_t scan_idx, int32_t pos_x,
                                int32_t pos_y, int32_t block_type, int8_t texture_type)
{
  const int32_t ctx_ind_map[16] =
  {
    0, 1, 4, 5,
    2, 3, 4, 5,
    6, 6, 8, 8,
    7, 7, 8, 8
  };

  int32_t cnt,offset,pos_x_in_subset,pos_y_in_subset;

  if (pos_x + pos_y == 0) return 0;

  if (block_type == 2) return ctx_ind_map[4 * pos_y + pos_x];

  offset = (block_type == 3) ? ((scan_idx == SCAN_DIAG) ? 9 : 15) : ((texture_type == 0) ? 21 : 12);
  pos_x_in_subset = pos_x - ((pos_x>>2)<<2);
  pos_y_in_subset = pos_y - ((pos_y>>2)<<2);

  if (pattern_sig_ctx == 0) {
    cnt = (pos_x_in_subset + pos_y_in_subset <= 2) ? ((pos_x_in_subset + pos_y_in_subset==0) ? 2 : 1) : 0;
  } else if (pattern_sig_ctx==1) {
    cnt = (pos_y_in_subset <= 1) ? ((pos_y_in_subset == 0) ? 2 : 1) : 0;
  } else if (pattern_sig_ctx==2) {
    cnt = (pos_x_in_subset <= 1) ? ((pos_x_in_subset == 0) ? 2 : 1) : 0;
  } else {
    cnt = 2;
  }
  return (( texture_type == 0 && ((pos_x>>2) + (pos_y>>2)) > 0 ) ? 3 : 0) + offset + cnt;
}


/**
* \brief Context derivation process of coeff_abs_significant_flag
* \param coeff     pointer to the current coefficient
* \param pos_x     column of current scan position
* \param pos_y     row of current scan position
* \param width     width of the block
* \param height    height of the block
* \param type      texture type (TEXT_LUMA...)
* \param temp_diag temporary output value used in the next steps
* \param temp_sum  temporary output value used in the next steps
* \returns context index for current scan position
*/
uint32_t kvz_context_get_sig_ctx_idx_abs(const coeff_t* coeff, int32_t pos_x, int32_t pos_y,
                                         uint32_t height, uint32_t width, int8_t type,
                                         int32_t* temp_diag, int32_t* temp_sum)
{
  const coeff_t* data = coeff + pos_x + pos_y * width;
  const int     diag = pos_x + pos_y;
  int           num_pos = 0;
  int           sum_abs = 0;
#define UPDATE(x) {int a=abs(x);sum_abs+=MIN(4-(a&1),a);num_pos+=(a?1:0);}
  if (pos_x < width - 1)
  {
    UPDATE(data[1]);
    if (pos_x < width - 2)
    {
      UPDATE(data[2]);
    }
    if (pos_y < height - 1)
    {
      UPDATE(data[width + 1]);
    }
  }
  if (pos_y < height - 1)
  {
    UPDATE(data[width]);
    if (pos_y < height - 2)
    {
      UPDATE(data[width << 1]);
    }
  }
#undef UPDATE
  int ctx_ofs = MIN(sum_abs, 5) + (diag < 2 ? 6 : 0);
  if (type == 0 /* Luma */)
  {
    ctx_ofs += diag < 5 ? 6 : 0;
  }
  
  *temp_diag = diag;
  *temp_sum = sum_abs - num_pos;
  return ctx_ofs;
}

/**
* \brief Calculate Go rice parameter for remainder coefficient value coding
* \param coeff     pointer to the current coefficient
* \param pos_x     column of current scan position
* \param pos_y     row of current scan position
* \param width     width of the block
* \param height    height of the block
* \returns context go rice parameter
*/
uint32_t kvz_go_rice_par_abs(const coeff_t* coeff, int32_t pos_x, int32_t pos_y,
                             uint32_t height, uint32_t width)
{
#define UPDATE(x) sum+=abs(x)-(x?1:0)

  const coeff_t* data = coeff + pos_x + pos_y * width;
  int           sum = 0;
  if (pos_x < width - 1)
  {
    UPDATE(data[1]);
    if (pos_x < width - 2)
    {
      UPDATE(data[2]);
    }
    if (pos_y < height - 1)
    {
      UPDATE(data[width + 1]);
    }
  }
  if (pos_y < height - 1)
  {
    UPDATE(data[width]);
    if (pos_y < height - 2)
    {
      UPDATE(data[width << 1]);
    }
  }
#undef UPDATE
  return  g_go_rice_pars[MIN(sum, 31)];
}