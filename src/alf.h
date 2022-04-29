#ifndef ALF_H_
#define ALF_H_

#include "checkpoint.h"
#include "cu.h"
#include "encoder.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "uvg266.h"
#include "videoframe.h"
#include "image.h"
#include "nal.h"


//ALF applied to 4x4 luma samples
//Filtering in CTB level
//Signalled in slice header if used or not
//Signalled in CTB level if ALF used for luma (if applied, then also for both chroma samples)

//To reduce bits overhead, filter coefficients of different classification can be merged.
//In slice header, the indices of the APSs used for the current slice are signaled.

#define MAX_NUM_ALF_CLASSES             25
#define MAX_NUM_ALF_LUMA_COEFF          13
#define MAX_NUM_ALF_CHROMA_COEFF        7
#define MAX_NUM_CC_ALF_FILTERS          4
#define MAX_NUM_CC_ALF_CHROMA_COEFF     8
#define ALF_NUM_FIXED_FILTER_SETS       16
#define ALF_UNUSED_CLASS_IDX            255
#define ALF_UNUSED_TRANSPOSE_IDX        255
#define CLASSIFICATION_BLK_SIZE         32
#define ALF_VB_POS_ABOVE_CTUROW_LUMA    4
#define ALF_VB_POS_ABOVE_CTUROW_CHMA    2
#define ALF_CTB_MAX_NUM_APS             8
#define MAX_ALF_PADDING_SIZE            4
#define MAX_ALF_NUM_CLIPPING_VALUES     4
#define MAX_NUM_ALF_ALTERNATIVES_CHROMA 8
#define CCALF_CANDS_COEFF_NR            8
#define CCALF_DYNAMIC_RANGE             6
#define CCALF_BITS_PER_COEFF_LEVEL      3
#define NUM_APS_TYPE_LEN                0 //1 //3
#define CC_ALF_NUM_COEFF                8

static const int cc_alf_small_tab[CCALF_CANDS_COEFF_NR] = { 0, 1, 2, 4, 8, 16, 32, 64 };

static const int g_fixed_filter_set_coeff[64][MAX_NUM_ALF_LUMA_COEFF] =
{
  { 0,   0,   2,  -3,   1,  -4,   1,   7,  -1,   1,  -1,   5,  0 },
  { 0,   0,   0,   0,   0,  -1,   0,   1,   0,   0,  -1,   2,  0 },
  { 0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,  0 },
  { 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   1,  0 },
  { 2,   2,  -7,  -3,   0,  -5,  13,  22,  12,  -3,  -3,  17,  0 },
  {-1,   0,   6,  -8,   1,  -5,   1,  23,   0,   2,  -5,  10,  0 },
  { 0,   0,  -1,  -1,   0,  -1,   2,   1,   0,   0,  -1,   4,  0 },
  { 0,   0,   3, -11,   1,   0,  -1,  35,   5,   2,  -9,   9,  0 },
  { 0,   0,   8,  -8,  -2,  -7,   4,   4,   2,   1,  -1,  25,  0 },
  { 0,   0,   1,  -1,   0,  -3,   1,   3,  -1,   1,  -1,   3,  0 },
  { 0,   0,   3,  -3,   0,  -6,   5,  -1,   2,   1,  -4,  21,  0 },
  {-7,   1,   5,   4,  -3,   5,  11,  13,  12,  -8,  11,  12,  0 },
  {-5,  -3,   6,  -2,  -3,   8,  14,  15,   2,  -7,  11,  16,  0 },
  { 2,  -1,  -6,  -5,  -2,  -2,  20,  14,  -4,   0,  -3,  25,  0 },
  { 3,   1,  -8,  -4,   0,  -8,  22,   5,  -3,   2, -10,  29,  0 },
  { 2,   1,  -7,  -1,   2, -11,  23,  -5,   0,   2, -10,  29,  0 },
  {-6,  -3,   8,   9,  -4,   8,   9,   7,  14,  -2,   8,   9,  0 },
  { 2,   1,  -4,  -7,   0,  -8,  17,  22,   1,  -1,  -4,  23,  0 },
  { 3,   0,  -5,  -7,   0,  -7,  15,  18,  -5,   0,  -5,  27,  0 },
  { 2,   0,   0,  -7,   1, -10,  13,  13,  -4,   2,  -7,  24,  0 },
  { 3,   3, -13,   4,  -2,  -5,   9,  21,  25,  -2,  -3,  12,  0 },
  {-5,  -2,   7,  -3,  -7,   9,   8,   9,  16,  -2,  15,  12,  0 },
  { 0,  -1,   0,  -7,  -5,   4,  11,  11,   8,  -6,  12,  21,  0 },
  { 3,  -2,  -3,  -8,  -4,  -1,  16,  15,  -2,  -3,   3,  26,  0 },
  { 2,   1,  -5,  -4,  -1,  -8,  16,   4,  -2,   1,  -7,  33,  0 },
  { 2,   1,  -4,  -2,   1, -10,  17,  -2,   0,   2, -11,  33,  0 },
  { 1,  -2,   7, -15, -16,  10,   8,   8,  20,  11,  14,  11,  0 },
  { 2,   2,   3, -13, -13,   4,   8,  12,   2,  -3,  16,  24,  0 },
  { 1,   4,   0,  -7,  -8,  -4,   9,   9,  -2,  -2,   8,  29,  0 },
  { 1,   1,   2,  -4,  -1,  -6,   6,   3,  -1,  -1,  -3,  30,  0 },
  {-7,   3,   2,  10,  -2,   3,   7,  11,  19,  -7,   8,  10,  0 },
  { 0,  -2,  -5,  -3,  -2,   4,  20,  15,  -1,  -3,  -1,  22,  0 },
  { 3,  -1,  -8,  -4,  -1,  -4,  22,   8,  -4,   2,  -8,  28,  0 },
  { 0,   3, -14,   3,   0,   1,  19,  17,   8,  -3,  -7,  20,  0 },
  { 0,   2,  -1,  -8,   3,  -6,   5,  21,   1,   1,  -9,  13,  0 },
  {-4,  -2,   8,  20,  -2,   2,   3,   5,  21,   4,   6,   1,  0 },
  { 2,  -2,  -3,  -9,  -4,   2,  14,  16,   3,  -6,   8,  24,  0 },
  { 2,   1,   5, -16,  -7,   2,   3,  11,  15,  -3,  11,  22,  0 },
  { 1,   2,   3, -11,  -2,  -5,   4,   8,   9,  -3,  -2,  26,  0 },
  { 0,  -1,  10,  -9,  -1,  -8,   2,   3,   4,   0,   0,  29,  0 },
  { 1,   2,   0,  -5,   1,  -9,   9,   3,   0,   1,  -7,  20,  0 },
  { -2,  8,  -6,  -4,   3,  -9,  -8,  45,  14,   2, -13,   7,  0 },
  { 1,  -1,  16, -19,  -8,  -4,  -3,   2,  19,   0,   4,  30,  0 },
  { 1,   1,  -3,   0,   2, -11,  15,  -5,   1,   2,  -9,  24,  0 },
  { 0,   1,  -2,   0,   1,  -4,   4,   0,   0,   1,  -4,   7,  0 },
  { 0,   1,   2,  -5,   1,  -6,   4,  10,  -2,   1,  -4,  10,  0 },
  { 3,   0,  -3,  -6,  -2,  -6,  14,   8,  -1,  -1,  -3,  31,  0 },
  { 0,   1,   0,  -2,   1,  -6,   5,   1,   0,   1,  -5,  13,  0 },
  { 3,   1,   9, -19, -21,   9,   7,   6,  13,   5,  15,  21,  0 },
  { 2,   4,   3, -12, -13,   1,   7,   8,   3,   0,  12,  26,  0 },
  { 3,   1,  -8,  -2,   0,  -6,  18,   2,  -2,   3, -10,  23,  0 },
  { 1,   1,  -4,  -1,   1,  -5,   8,   1,  -1,   2,  -5,  10,  0 },
  { 0,   1,  -1,   0,   0,  -2,   2,   0,   0,   1,  -2,   3,  0 },
  { 1,   1,  -2,  -7,   1,  -7,  14,  18,   0,   0,  -7,  21,  0 },
  { 0,   1,   0,  -2,   0,  -7,   8,   1,  -2,   0,  -3,  24,  0 },
  { 0,   1,   1,  -2,   2, -10,  10,   0,  -2,   1,  -7,  23,  0 },
  { 0,   2,   2, -11,   2,  -4,  -3,  39,   7,   1, -10,   9,  0 },
  { 1,   0,  13, -16,  -5,  -6,  -1,   8,   6,   0,   6,  29,  0 },
  { 1,   3,   1,  -6,  -4,  -7,   9,   6,  -3,  -2,   3,  33,  0 },
  { 4,   0, -17,  -1,  -1,   5,  26,   8,  -2,   3, -15,  30,  0 },
  { 0,   1,  -2,   0,   2,  -8,  12,  -6,   1,   1,  -6,  16,  0 },
  { 0,   0,   0,  -1,   1,  -4,   4,   0,   0,   0,  -3,  11,  0 },
  { 0,   1,   2,  -8,   2,  -6,   5,  15,   0,   2,  -7,   9,  0 },
  { 1,  -1,  12, -15,  -7,  -2,   3,   6,   6,  -1,   7,  30,  0 },
};

static const int g_class_to_filter_mapping[ALF_NUM_FIXED_FILTER_SETS][MAX_NUM_ALF_CLASSES] =
{
  { 8,   2,   2,   2,   3,   4,  53,   9,   9,  52,   4,   4,   5,   9,   2,   8,  10,   9,   1,   3,  39,  39,  10,   9,  52 },
  { 11,  12,  13,  14,  15,  30,  11,  17,  18,  19,  16,  20,  20,   4,  53,  21,  22,  23,  14,  25,  26,  26,  27,  28,  10 },
  { 16,  12,  31,  32,  14,  16,  30,  33,  53,  34,  35,  16,  20,   4,   7,  16,  21,  36,  18,  19,  21,  26,  37,  38,  39 },
  { 35,  11,  13,  14,  43,  35,  16,   4,  34,  62,  35,  35,  30,  56,   7,  35,  21,  38,  24,  40,  16,  21,  48,  57,  39 },
  { 11,  31,  32,  43,  44,  16,   4,  17,  34,  45,  30,  20,  20,   7,   5,  21,  22,  46,  40,  47,  26,  48,  63,  58,  10 },
  { 12,  13,  50,  51,  52,  11,  17,  53,  45,   9,  30,   4,  53,  19,   0,  22,  23,  25,  43,  44,  37,  27,  28,  10,  55 },
  { 30,  33,  62,  51,  44,  20,  41,  56,  34,  45,  20,  41,  41,  56,   5,  30,  56,  38,  40,  47,  11,  37,  42,  57,   8 },
  { 35,  11,  23,  32,  14,  35,  20,   4,  17,  18,  21,  20,  20,  20,   4,  16,  21,  36,  46,  25,  41,  26,  48,  49,  58 },
  { 12,  31,  59,  59,   3,  33,  33,  59,  59,  52,   4,  33,  17,  59,  55,  22,  36,  59,  59,  60,  22,  36,  59,  25,  55 },
  { 31,  25,  15,  60,  60,  22,  17,  19,  55,  55,  20,  20,  53,  19,  55,  22,  46,  25,  43,  60,  37,  28,  10,  55,  52 },
  { 12,  31,  32,  50,  51,  11,  33,  53,  19,  45,  16,   4,   4,  53,   5,  22,  36,  18,  25,  43,  26,  27,  27,  28,  10 },
  { 5,   2,  44,  52,   3,   4,  53,  45,   9,   3,   4,  56,   5,   0,   2,   5,  10,  47,  52,   3,  63,  39,  10,   9,  52 },
  { 12,  34,  44,  44,   3,  56,  56,  62,  45,   9,  56,  56,   7,   5,   0,  22,  38,  40,  47,  52,  48,  57,  39,  10,   9 },
  { 35,  11,  23,  14,  51,  35,  20,  41,  56,  62,  16,  20,  41,  56,   7,  16,  21,  38,  24,  40,  26,  26,  42,  57,  39 },
  { 33,  34,  51,  51,  52,  41,  41,  34,  62,   0,  41,  41,  56,   7,   5,  56,  38,  38,  40,  44,  37,  42,  57,  39,  10 },
  { 16,  31,  32,  15,  60,  30,   4,  17,  19,  25,  22,  20,   4,  53,  19,  21,  22,  46,  25,  55,  26,  48,  63,  58,  55 },
};


//-------------------------typedef enums----------------------------
typedef enum { ALF_FILTER_5X5 = 0, ALF_FILTER_7X7 = 1, ALF_NUM_OF_FILTER_TYPES = 2 } alf_filter_type;
typedef enum { ALF_LUMA = 0, ALF_CHROMA = 1 } alf_type;

typedef enum {
  T_ALF_APS = 0,
  T_LMCS_APS = 1,
  T_SCALING_LIST_APS = 2,
} aps_type;

typedef enum {
  ALF_HOR = 0,
  ALF_VER = 1,
  ALF_DIAG0 = 2,
  ALF_DIAG1 = 3,
  NUM_DIRECTIONS = 4
} alf_directions;

typedef enum {
  CHANNEL_TYPE_LUMA = 0,
  CHANNEL_TYPE_CHROMA = 1,
  MAX_NUM_CHANNEL_TYPE = 2
} channel_type;

typedef enum {
  COMPONENT_Y = 0,
  COMPONENT_Cb = 1,
  COMPONENT_Cr = 2,
  MAX_NUM_COMPONENT = 3,
} alf_component_id;

//----------------------------------------------------------------

//-------------------------typedef structs----------------------------
#if _MSC_VER
#define PACK(__Type__,__Declaration__) __pragma(pack(push, 1)) __Type__ __Declaration__ __pragma(pack(pop));
#else
#define PACK(__Type__, __Declaration__) __Type__ __attribute__((__packed__)) __Declaration__ ;
#endif


PACK(typedef struct, alf_covariance {
  double pix_acc;
  int64_t ee[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF][MAX_ALF_NUM_CLIPPING_VALUES][MAX_ALF_NUM_CLIPPING_VALUES];
  int32_t y[MAX_NUM_ALF_LUMA_COEFF][MAX_ALF_NUM_CLIPPING_VALUES];
  int num_coeff;
  int num_bins;
} alf_covariance)

typedef struct clp_rng {
  int min;
  int max;
  int bd;
  int n;
} clp_rng;

typedef struct clp_rngs {
  clp_rng comp[MAX_NUM_COMPONENT]; ///< the bit depth as indicated in the SPS
  bool used;
  bool chroma;
} clp_rngs;

typedef struct alf_classifier {
  uint8_t class_idx;
  uint8_t transpose_idx;
} alf_classifier;

typedef struct cc_alf_filter_param {
  bool    cc_alf_filter_enabled[2];
  bool    cc_alf_filter_idx_enabled[2][MAX_NUM_CC_ALF_FILTERS];
  uint8_t cc_alf_filter_count[2];
  short   cc_alf_coeff[2][MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF];
  int     new_cc_alf_filter[2];
  int     number_valid_components;
} cc_alf_filter_param;

typedef struct alf_aps {
  int aps_id;
  int temporal_id;
  int layer_id;
  aps_type aps_type;                  // aps_params_type

  //sliceparams
  bool enabled_flag[MAX_NUM_COMPONENT];                           // alf_slice_enable_flag, alf_chroma_idc
  bool non_linear_flag[MAX_NUM_CHANNEL_TYPE];                     // alf_[luma/chroma]_clip_flag

  short luma_coeff[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF]; // alf_coeff_luma_delta[i][j]
  int16_t luma_clipp[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF]; // alf_clipp_luma_[i][j]

  int num_alternatives_chroma;                                                  // alf_chroma_num_alts_minus_one + 1
  short chroma_coeff[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF]; // alf_coeff_chroma[i]
  int16_t chroma_clipp[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF]; // alf_clipp_chroma[i]

  short filter_coeff_delta_idx[MAX_NUM_ALF_CLASSES];              // filter_coeff_delta[i]
  bool alf_luma_coeff_flag[MAX_NUM_ALF_CLASSES];                  // alf_luma_coeff_flag[i]
  int num_luma_filters;                                           // number_of_filters_minus1 + 1
  bool alf_luma_coeff_delta_flag;                                 // alf_luma_coeff_delta_flag
  bool new_filter_flag[MAX_NUM_CHANNEL_TYPE];

  cc_alf_filter_param cc_alf_aps_param;

} alf_aps;

typedef struct alf_info_t {
  cabac_data_t cabac_estimator;

  uvg_pixel *alf_fulldata_buf;
  uvg_pixel *alf_fulldata;
  uvg_pixel *alf_tmp_y;
  uvg_pixel *alf_tmp_u;
  uvg_pixel *alf_tmp_v;

  alf_covariance* alf_covariance; //Covariances of each CTU for luma and chroma components //[ctu_idx][class_idx]
  alf_covariance* alf_covariance_y; //Pointer to the first luma covaraince //[ctu_idx][class_idx]
  alf_covariance* alf_covariance_u; //Pointer to the first Cb covariance //[ctu_idx][class_idx]
  alf_covariance* alf_covariance_v; //Pointer tot he first Cr covariance //[ctu_idx][class_idx]
  alf_covariance alf_covariance_frame_luma[MAX_NUM_ALF_CLASSES]; //[class_idx]
  alf_covariance alf_covariance_frame_chroma[MAX_NUM_ALF_ALTERNATIVES_CHROMA]; //[class_idx]
  alf_covariance alf_covariance_merged[MAX_NUM_ALF_CLASSES + 2];
  alf_covariance* alf_covariance_cc_alf[MAX_NUM_COMPONENT]; // [compIdx-1][filterIdx][ctbAddr]
  alf_covariance alf_covariance_frame_cc_alf[MAX_NUM_COMPONENT - 1][MAX_NUM_CC_ALF_FILTERS];

  bool *ctu_enable_flag[MAX_NUM_COMPONENT + 1];
  bool *ctu_enable_flag_tmp[MAX_NUM_COMPONENT + 1];
  uint8_t* ctu_alternative[MAX_NUM_COMPONENT + 1];
  uint8_t* ctu_alternative_tmp[MAX_NUM_COMPONENT + 1];
  double *ctb_distortion_unfilter[MAX_NUM_COMPONENT + 1];

  int aps_id_start;

  int alf_clip_merged[ALF_NUM_OF_FILTER_TYPES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF];
  int** diff_filter_coeff; // [lumaClassIdx][coeffIdx]
  int** filter_coeff_set;  // [lumaClassIdx][coeffIdx]
  int** filter_clipp_set; // [lumaClassIdx][coeffIdx]
  short* alf_ctb_filter_set_index_tmp; //g_num_ctus_in_pic //voisi olla lokaali muuttuja?
  short* alf_ctb_filter_index;     //g_num_ctus_in_pic

  uint8_t* training_cov_control; //[ctuAddr] 
  uint64_t* training_distortion[MAX_NUM_CC_ALF_FILTERS + 1]; //[ctuAddr] 
  uint8_t* filter_control; //[ctuAddr] 
  uint8_t* best_filter_control; //[ctuAddr] 
  uint8_t* cc_alf_filter_control[3]; //[ctuAddr] 

  alf_classifier **classifier;
  alf_aps alf_param_temp;

} alf_info_t;

typedef struct param_set_map {
  bool b_changed;
  //uint8_t* p_nalu_data;
  struct alf_aps parameter_set;
} param_set_map;

typedef struct array_variables {
  short fixed_filter_set_coeff_dec[ALF_NUM_FIXED_FILTER_SETS][MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
  short chroma_coeff_final[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF];
  short coeff_final[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
  short coeff_aps_luma[ALF_CTB_MAX_NUM_APS][MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];

  int16_t chroma_clipp_final[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF];
  int16_t clip_default[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
  int16_t clipp_final[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
  int16_t clipp_aps_luma[ALF_CTB_MAX_NUM_APS][MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];

  short filter_indices[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES];

  unsigned bits_new_filter[MAX_NUM_CHANNEL_TYPE];
  short alf_clipping_values[MAX_NUM_CHANNEL_TYPE][MAX_ALF_NUM_CLIPPING_VALUES];
  int cc_reuse_aps_id[2];

  int filter_coeff_set[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF];
  int filter_clipp_set[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF];

  struct clp_rngs clp_rngs;

} array_variables;

//inits aps parameter set in videoframe
void uvg_set_aps_map(videoframe_t* frame, enum uvg_alf alf_type);

//resets cc alf parameter
void uvg_reset_cc_alf_aps_param(cc_alf_filter_param *cc_alf);

//starts alf encoding process
void uvg_alf_enc_process(encoder_state_t *const state);

//creates variables for alf_info_t structure in videoframe_t 
void uvg_alf_create(videoframe_t *frame, enum uvg_chroma_format chroma_format);
//frees allocated memory in alf_info_t structure
void uvg_alf_destroy(videoframe_t * const frame);

//writes alf bits to bitstream
void uvg_encode_alf_bits(encoder_state_t * const state, const int ctu_idx);

//writes apss to header
void uvg_encode_alf_adaptive_parameter_set(encoder_state_t * const state);

#endif //ALF_H_