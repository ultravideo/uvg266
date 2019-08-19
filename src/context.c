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
  {  18,  27,  15,  11,  28,  30,  19,  22,  23, },
  {  18,  27,  53,  12,   6,  30,  13,  15,  31, },
  {  19,  28,  38,  12,  29,  38,  28,  38,  31, },
  {  12,  13,   8,   8,  13,  12,   5,  10,   9, },
};

static const uint8_t  INIT_QT_SPLIT_FLAG[4][6] = {
  {  26,  36,  38,  33,  34,  21, },
  {  20,   7,  23,  18,  19,   6, },
  {  12,   6,  15,  33,  27,  22, },
  {   0,   8,   8,  12,  12,  12, },
};

static const uint8_t INIT_SKIP_FLAG[4][3] = {
  {  57,  60,  53, },
  {  57,  59,  45, },
  {   0,  34,  36, },
  {   5,   4,   8, },
};

static const uint8_t INIT_MERGE_FLAG_EXT[4][1] = {
  {   6, },
  {   6, },
  {  19, },
  {   5, },
};


static const uint8_t INIT_MERGE_IDX_EXT[4][1] = {
  {  41, },
  {  43, },
  {  34, },
  {   4, },
};

static const uint8_t INIT_PART_SIZE[4][4] = {
  {  CNU, CNU, CNU, CNU,},
  {  CNU, CNU, CNU, CNU,},
  {  CNU, CNU, CNU, CNU,},
  { DWS, DWS, DWS, DWS, }
};

static const uint8_t INIT_PRED_MODE[4][2] = {
  {  40,  50, },
  {  40,  35, },
  { CNU, CNU, },
  {   5,   2, },
};

static const uint8_t MULTI_REF_LINE_MODE[4][2] = {
  {  25,  50, },
  {  25,  57, },
  {  25,  51, },
  {   6,   8, },
};

static const uint8_t INIT_INTRA_PRED_MODE[4] = {
  154, 154, 170, 6
};

static const uint8_t INIT_INTRA_LUMA_PLANAR_MODE[4][2] = {
  {  13,  21, },
  {  12,  13, },
  {  13,  28, },
  {   4,   5, },
 };

static const uint8_t INIT_CHROMA_PRED_MODE[4][1] = {
  {  25, },
  {  33, },
  {  19, },
  {   6, },
};

static const uint8_t INIT_CU_QP_DELTA_ABS[4][2] = {
  { CNU, CNU, },
  { CNU, CNU, },
  { CNU, CNU, },
  { DWS, DWS, },
};

static const uint8_t INIT_INTER_DIR[4][5] = {
  {   6,  13,   5,   4,  25, },
  {   7,   6,   5,   4,  33, },
  { CNU, CNU, CNU, CNU, CNU, },
  {   0,   0,   1,   4,   0, },
};

static const uint8_t INIT_REF_PIC[4][2] = {
  {  13,  20, },
  {  27,  35, },
  { CNU, CNU, },
  {   0,   4, },
};

static const uint8_t INIT_MVD[4][2] = {
  {  51,  58, },
  {  44,  43, },
  {  14,  45, },
  {   9,   5, },
};

static const uint8_t INIT_QT_ROOT_CBF[4][1] = {
  {  12, },
  {   5, },
  {   6, },
  {   4, },
};


static const uint8_t INIT_QT_CBF[4][7] = {
    {  15, CNU,   5,  14, 25,   9,  44},
    {  15, CNU,  20,   7, 25,  25,  29},
    {   7, CNU,   5,   7, 12,  33,  21},
    {   5, DWS,   8,   8,  5,   2,   1},
};

static const uint8_t BDPCM_MODE_INIT[4][2] = {
  { CNU, CNU, },
  { CNU, CNU, },
  { CNU, CNU, },
  { DWS, DWS, },
};

static const uint8_t INIT_SIG_COEFF_GROUP[4][4] = {
    {  25,  37,   25,  37},
    {  25,  30,   25,  52},
    {  18,  31,   25,   7},
    {   8,   5,    5,   8},
};

static const uint8_t INIT_SIG_FLAG[6][4][12] = {
  {
    {  17,  41,  49,  51,   1,  49,  50,  37,  48,  51,  58,  45, },
    {  17,  41,  42,  29,  25,  49,  43,  37,  33,  51,  51,  30, },
    {  25,  19,  28,  14,  25,  20,  29,  30,  19,  52,  30,  38, },
    {  12,   9,   9,  10,   9,   9,   9,  10,   8,   8,   8,  10, },
  },
  {
    {   9,  49,  42,  21,  48,  59,  59,  53, },
    {  17,  19,  20,  29,  41,  59,  60,  38, },
    {  25,  27,  28,  37,  49,  53,  53,  46, },
    {   9,   9,   9,  13,   5,   5,   8,   9, },
  },
  {
    {  26,  45,  53,  46,  49,  54,  61,  39,  42,  39,  39,  39, },
    {  19,  38,  38,  46,  34,  54,  54,  39,   6,  39,  39,  39, },
    {  11,  38,  46,  54,  27,  39,  39,  39,  28,  39,  39,  39, },
    {   9,  12,   8,   8,   8,   8,   8,   5,   8,   0,   0,   0, },
  },
  {
    {  34,  45,  38,  31,  58,  39,  39,  39, },
    {  35,  45,  53,  54,  51,  39,  39,  39, },
    {  19,  46,  38,  39,  52,  39,  39,  39, },
    {   8,  12,   8,   8,   4,   0,   0,   0, },
  },
  {
    {  19,  54,  39,  39,  50,  39,  39,  39,   0,  39,  39,  39, },
    {  19,  39,  54,  39,  19,  39,  39,  39,  56,  39,  39,  39, },
    {  18,  39,  39,  39,  11,  39,  39,  39,   0,  39,  39,  39, },
    {   8,   8,   8,   8,   8,   0,   4,   4,   0,   0,   0,   0, },
  },
  {
    {  34,  38,  54,  39,  41,  39,  39,  39, },
    {  34,  38,  62,  39,  26,  39,  39,  39, },
    {  26,  39,  39,  39,  19,  39,  39,  39, },
    {   8,   8,   8,   8,   4,   0,   0,   0, },
  }
};

static const uint8_t INIT_PARITY_FLAG[2][4][21] =
{
  {
    {  33,  40,  25,  41,  26,  42,  25,  33,  26,  34,  27,  25,  41,  42,  42,  35,  33,  27,  35,  42,  43, },
    {  18,  17,  33,  18,  34,  42,  25,  33,  26,  42,  27,  25,  34,  42,  42,  20,  26,  27,  42,  20,  20, },
    {  33,  25,  18,  26,  34,  27,  25,  26,  19,  42,  35,  33,  19,  27,  35,  20,  34,  42,  20,  43,  20, },
    {   8,   9,  12,  13,  13,  13,  10,  13,  13,  13,  13,  13,  13,  13,  13,  13,  10,  13,  13,  13,  13, },
  },
  {
    {  33,  25,  26,  19,  19,  27,  33,  42,  43,  27,  43, },
    {  25,  25,  26,  11,  19,  27,  33,  42,  50,  20,  43, },
    {  33,  25,  26,  42,  19,  27,  26,  50,  43,  20,  43, },
    {   9,  13,  12,  12,  13,  13,  13,  13,  13,  13,  13, },
  }
};

static const uint8_t INIT_GTX_FLAG[4][4][21] =
{
  {
    {  25,   0,   0,  17,  25,  18,   0,   9,  25,  33,  19,   0,  25,  33,  26,  20,  25,  33,  34,  42,  29, },
    {  17,   0,   1,  17,  25,  18,   0,   9,  25,  33,  34,   9,  25,  18,  26,  20,  25,  18,  19,  27,  21, },
    {  25,   1,  40,  25,  33,  11,  17,  25,  25,  18,   4,  17,  33,  11,   4,   5,  33,  19,  20,  28,  22, },
    {   1,   5,   9,   9,   9,   6,   5,   9,  10,  10,   9,   9,   9,   9,   9,   9,   6,   8,   9,   8,   9, },
  },
  {
    {  25,   1,  40,  33,  18,   4,  25,  33,  27,  36,  37, },
    {  17,   9,  25,  10,   3,   4,  17,  33,  19,  28,  29, },
    {  48,   9,  25,  18,  26,  27,  25,  26,  35,  28,  37, },
    {   1,   5,   8,   8,   8,   6,   6,   9,   8,   8,   9, },
  },
  {
    {   0,   0,  33,  34,  35,  36,  25,  49,  35,  28,  29,  40,  42,  43,  36,  37,  56,  58,  59,  45,  38, },
    {   0,  17,  26,  19,  20,  21,  25,  34,  20,  28,  29,  33,  27,  28,  29,  37,  34,  28,  44,  37,  38, },
    {  25,  25,  11,  27,  20,  21,  18,  12,  28,  21,  22,  34,  28,  29,  29,  30,  28,  29,  45,  30,  23, },
    {   9,   5,  10,  13,  13,  10,   9,  10,  13,  13,  13,   9,  10,  10,  10,  10,   8,   9,   8,  10,  13, },
  },
  {
    {   0,  40,  42,  20,  21,  29,  49,  52,  53,  38,  46, },
    {   0,  25,  27,  20,  13,   6,  57,  52,  30,  38,  31, },
    {  40,  33,  27,  28,  21,  37,  51,  37,  53,  38,  46, },
    {   9,   9,  10,  12,  12,  10,   5,   9,   9,   9,   9, },
  }
};

static const uint8_t INIT_LAST_X[4][23] = {
    {  14,   6,   5,   7,   7,  12,   7,   7,   6,  12,  22,   7,   6,  14,  20,  28,   7,  13,  13,  20,   11,   5,   3},
    {   6,  13,  12,   6,   6,   4,  14,  14,   5,  12,  29,  14,  13,   5,  36,  28,  14,  13,  20,  19,   19,   4,  18},
    {   6,  13,  12,   6,  14,  12,  14,  14,  29,   4,  14,   7,  14,  29,   4,  29,  30,  37,  29,  58,   12,  11,   3},
    {   8,   5,   4,   5,   4,   4,   5,   4,   1,   0,   5,   1,   0,   0,   0,   1,   1,   0,   0,   0,    2,   1,   1},
};

static const uint8_t INIT_LAST_Y[4][23] = {
    {  13,  13,  20,   6,   6,  12,  14,  14,   5,  13,  14,   7,   5,  12,  21,  13,   7,  13,  12,  41,   11,   5,  19},
    {   5,   5,  12,   6,   6,  19,   6,  14,   5,  19,  29,   7,  13,   5,  36,  21,   7,  13,   5,  27,   11,   4,  18},
    {  13,   5,   4,   6,   6,  11,  14,  14,   5,  11,  14,   7,  14,   5,   3,  21,  45,  45,  21,  34,   12,   4,   3},
    {   8,   5,   8,   5,   5,   4,   5,   5,   4,   0,   5,   5,   1,   0,   0,   1,   4,   0,   0,   0,    6,   2,   2},
};

static const uint8_t INIT_MVP_IDX[4][1] = {
  {  34, },
  {  49, },
  {  42, },
  {  12, },
};

static const uint8_t INIT_JOINT_CB_CR_FLAG[4][3] = { 
  {  51,  44,  45, },
  {  36,  44,  45, },
  {  43,  29,  51, },
  {   1,   1,   0, }, 
};

static const uint8_t INIT_SAO_MERGE_FLAG[4] = { 2, 60, 59, 0 };

static const uint8_t INIT_SAO_TYPE_IDX[4] = { 10, 5, 5, 0 };

static const uint8_t INIT_CU_TRANSQUANT_BYPASS[4][1] = {
  { CNU, },
  { CNU, },
  { CNU, },
  { DWS, },
};

static const uint8_t INIT_INTRA_SUBPART_MODE[4][2] = {
  {  48,  43, },
  {  33,  43, },
  {  33,  43, },
  {   9,   2, },
};

/*
static const uint16_t g_inistateToCount[128] = {
  614,   647,   681,   718,   756,   797,   839,   884,   932,   982,   1034,  1089,  1148,  1209,  1274,  1342,
  1414,  1490,  1569,  1653,  1742,  1835,  1933,  2037,  2146,  2261,  2382,  2509,  2643,  2785,  2934,  3091,
  3256,  3430,  3614,  3807,  4011,  4225,  4452,  4690,  4941,  5205,  5483,  5777,  6086,  6412,  6755,  7116,
  7497,  7898,  8320,  8766,  9235,  9729,  10249, 10798, 11375, 11984, 12625, 13300, 14012, 14762, 15551, 16384,
  16384, 17216, 18005, 18755, 19467, 20142, 20783, 21392, 21969, 22518, 23038, 23532, 24001, 24447, 24869, 25270,
  25651, 26012, 26355, 26681, 26990, 27284, 27562, 27826, 28077, 28315, 28542, 28756, 28960, 29153, 29337, 29511,
  29676, 29833, 29982, 30124, 30258, 30385, 30506, 30621, 30730, 30834, 30932, 31025, 31114, 31198, 31277, 31353,
  31425, 31493, 31558, 31619, 31678, 31733, 31785, 31835, 31883, 31928, 31970, 32011, 32049, 32086, 32120, 32153
};*/


/**
 * \brief Initialize struct cabac_ctx.
 */
void kvz_ctx_init(cabac_ctx_t *ctx, uint32_t qp, uint32_t init_value, uint8_t rate)
{

  int slope = (init_value >> 3) - 4;
  int offset = ((init_value & 7) * 18) + 1;
  int inistate = ((slope   * (qp - 16)) >> 1) + offset;
  int state_clip = inistate < 1 ? 1 : inistate > 127 ? 127 : inistate;
  const int p1 = (state_clip << 8);
  /*
  int slope = (init_value >> 4) * 5 - 45;
  int offset = ((init_value & 15) << 3) - 16;
  int init_state = ((slope * (int)qp) >> 4) + offset;
    
  const int p1 = g_inistateToCount[init_state < 0 ? 0 : init_state > 127 ? 127 : init_state];
  */

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

  

  kvz_ctx_init(&cabac->ctx.multi_ref_line[0], QP, MULTI_REF_LINE_MODE[slice][0], MULTI_REF_LINE_MODE[3][0]);
  kvz_ctx_init(&cabac->ctx.multi_ref_line[1], QP, MULTI_REF_LINE_MODE[slice][1], MULTI_REF_LINE_MODE[3][1]);

  kvz_ctx_init(&cabac->ctx.chroma_pred_model[0], QP, INIT_CHROMA_PRED_MODE[slice][0], INIT_CHROMA_PRED_MODE[3][0]);


  for (i = 0; i < 3; i++) {
    kvz_ctx_init(&cabac->ctx.cu_skip_flag_model[i], QP, INIT_SKIP_FLAG[slice][i], INIT_SKIP_FLAG[3][i]);
    kvz_ctx_init(&cabac->ctx.joint_bc_br[i], QP, INIT_JOINT_CB_CR_FLAG[slice][i], INIT_JOINT_CB_CR_FLAG[3][i]);   
    
  }

  for (i = 0; i < 4; i++) {
    kvz_ctx_init(&cabac->ctx.sig_coeff_group_model[i], QP, INIT_SIG_COEFF_GROUP[slice][i], INIT_SIG_COEFF_GROUP[3][i]);
  }


  for (i = 0; i < 6; i++) {
    kvz_ctx_init(&cabac->ctx.qt_split_flag_model[i], QP, INIT_QT_SPLIT_FLAG[slice][i], INIT_QT_SPLIT_FLAG[3][i]);
  }



  for (i = 0; i < 9; i++) {
    kvz_ctx_init(&cabac->ctx.split_flag_model[i], QP, INIT_SPLIT_FLAG[slice][i], INIT_SPLIT_FLAG[3][i]);
  }
  

  //TODO: ignore P/B contexts on intra frame
  kvz_ctx_init(&cabac->ctx.cu_qt_root_cbf_model, QP, INIT_QT_ROOT_CBF[slice][0], INIT_QT_ROOT_CBF[3][0]);
  kvz_ctx_init(&cabac->ctx.mvp_idx_model, QP, INIT_MVP_IDX[slice][0], INIT_MVP_IDX[3][0]);

  kvz_ctx_init(&cabac->ctx.qt_cbf_model_cb[0], QP, INIT_QT_CBF[slice][5], INIT_QT_CBF[3][5]);

  for (i = 0; i < 2; i++) {
    kvz_ctx_init(&cabac->ctx.qt_cbf_model_cr[i], QP, INIT_QT_CBF[slice][i + 6], INIT_QT_CBF[3][i + 6]);
    kvz_ctx_init(&cabac->ctx.cu_mvd_model[i], QP, INIT_MVD[slice][i], INIT_MVD[3][i]);
    kvz_ctx_init(&cabac->ctx.cu_ref_pic_model[i], QP, INIT_REF_PIC[slice][i], INIT_REF_PIC[3][i]);    
    kvz_ctx_init(&cabac->ctx.luma_planar_model[i], QP, INIT_INTRA_LUMA_PLANAR_MODE[slice][i], INIT_INTRA_LUMA_PLANAR_MODE[3][i]);
    kvz_ctx_init(&cabac->ctx.bdpcm_mode[i], QP, BDPCM_MODE_INIT[slice][i], BDPCM_MODE_INIT[3][i]);
    kvz_ctx_init(&cabac->ctx.cu_qp_delta_abs[i], QP, INIT_CU_QP_DELTA_ABS[slice][i], INIT_CU_QP_DELTA_ABS[3][i]);
  }

  for (i = 0; i < 3; i++) {
    kvz_ctx_init(&cabac->ctx.cu_ctx_last_y_chroma[i], QP, INIT_LAST_Y[slice][i + 20], INIT_LAST_Y[3][i + 20]);
    kvz_ctx_init(&cabac->ctx.cu_ctx_last_x_chroma[i], QP, INIT_LAST_X[slice][i + 20], INIT_LAST_X[3][i + 20]);
  }

  for (i = 0; i < 4; i++) {
    
    kvz_ctx_init(&cabac->ctx.part_size_model[i], QP, INIT_PART_SIZE[slice][i], INIT_PART_SIZE[3][i]);

    kvz_ctx_init(&cabac->ctx.qt_cbf_model_luma[i], QP, INIT_QT_CBF[slice][i], INIT_QT_CBF[3][i]);
  }

  for (i = 0; i < 5; i++) {  

    kvz_ctx_init(&cabac->ctx.inter_dir[i], QP, INIT_INTER_DIR[slice][i], INIT_INTER_DIR[3][i]);
  }

  for (i = 0; i < 20; i++) {
    kvz_ctx_init(&cabac->ctx.cu_ctx_last_y_luma[i], QP, INIT_LAST_Y[slice][i], INIT_LAST_Y[3][i]);
    kvz_ctx_init(&cabac->ctx.cu_ctx_last_x_luma[i], QP, INIT_LAST_X[slice][i], INIT_LAST_X[3][i]);
  }

  for (ii = 0; ii < 3; ii++) {
    for (i = 0; i < 12; i++) {
      kvz_ctx_init(&cabac->ctx.cu_sig_model_luma[ii][i], QP, INIT_SIG_FLAG[ii*2][slice][i], INIT_SIG_FLAG[ii * 2][3][i]);
      if (i < 8) kvz_ctx_init(&cabac->ctx.cu_sig_model_chroma[ii][i], QP, INIT_SIG_FLAG[ii*2+1][slice][i], INIT_SIG_FLAG[ii * 2 + 1][3][i]);
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
  uint32_t position = pos_y * width + pos_x;
  if (pos_x + 1 < (uint32_t)width) uiRight = sig_coeff_group_flag[position + 1];
  if (pos_y + 1 < (uint32_t)width) uiLower = sig_coeff_group_flag[position + width];

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
#define UPDATE(x) {int a=abs(x);sum_abs+=MIN(4+(a&1),a);num_pos+=(a?1:0);}
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
* \brief Calculate slot of Go rice parameter for remainder coefficient value coding
* \param coeff     pointer to the current coefficient
* \param pos_x     column of current scan position
* \param pos_y     row of current scan position
* \param width     width of the block
* \param height    height of the block
* \returns context go rice parameter
*/
uint32_t kvz_abs_sum(const coeff_t* coeff, int32_t pos_x, int32_t pos_y,
                             uint32_t height, uint32_t width, uint32_t baselevel)
{
#define UPDATE(x) sum+=abs(x)/*-(x?1:0)*/


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
  return  MAX(MIN(sum - 5 * (int32_t)baselevel, 31),0);
  /*return  MIN(sum, 31);*/
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
  uint32_t height, uint32_t width, uint32_t baselevel)
{
//#define UPDATE(x) sum+=abs(x)/*-(x?1:0)*/
//
//  const coeff_t* data = coeff + pos_x + pos_y * width;
//  int           sum = 0;
//  if (pos_x < width - 1)
//  {
//    UPDATE(data[1]);
//    if (pos_x < width - 2)
//    {
//      UPDATE(data[2]);
//    }
//    if (pos_y < height - 1)
//    {
//      UPDATE(data[width + 1]);
//    }
//  }
//  if (pos_y < height - 1)
//  {
//    UPDATE(data[width]);
//    if (pos_y < height - 2)
//    {
//      UPDATE(data[width << 1]);
//    }
//  }
//#undef UPDATE
  uint32_t check = kvz_abs_sum(coeff, pos_x, pos_y, height, width, baselevel);
  return  g_go_rice_pars[check];
  /*return  g_go_rice_pars[kvz_abs_sum(coeff, pos_x, pos_y, height, width, baselevel)];*/
}