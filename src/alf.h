#ifndef ALF_H_
#define ALF_H_

#include "checkpoint.h"
#include "cu.h"
#include "encoder.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "videoframe.h"
#include "image.h"
#include "nal.h"


//ALF applied to 4x4 luma samples
//Filtering in CTB level
//Signalled in slice header if used or not
//Signalled in CTB level if ALF used for luma (if applied, then also for both chroma samples)

//To reduce bits overhead, filter coefficients of different classification can be merged.
//In slice header, the indices of the APSs used for the current slice are signaled.

#define FULL_FRAME                      1
#define RUN_ALF_AFTER_FULL_FRAME        0

#define ALF_FIXED_FILTER_NUM            64
#define MAX_NUM_ALF_CLASSES             25
#define MAX_NUM_ALF_LUMA_COEFF          13
#define MAX_NUM_ALF_CHROMA_COEFF        7
#define MAX_ALF_FILTER_LENGTH           7
#define MAX_NUM_ALF_COEFF               (MAX_ALF_FILTER_LENGTH * MAX_ALF_FILTER_LENGTH / 2 + 1)
#define MAX_NUM_CC_ALF_FILTERS          4
#define MAX_NUM_CC_ALF_CHROMA_COEFF     8
#define ALF_NUM_FIXED_FILTER_SETS       16
#define CC_ALF_BITS_PER_COEFF_LEVEL     3
#define ALF_UNUSED_CLASS_IDX            255
#define ALF_UNUSED_TRANSPOSE_IDX        255
#define REG                             0.0001
#define REG_SQR                         0.0000001
#define MAX_NUM_APS                     32
#define CLASSIFICATION_BLK_SIZE         32
#define ALF_VB_POS_ABOVE_CTUROW_LUMA    4
#define ALF_VB_POS_ABOVE_CTUROW_CHMA    2
#define ALF_CTB_MAX_NUM_APS             8
#define MAX_ALF_PADDING_SIZE            4
#define MAX_ALF_NUM_CLIPPING_VALUES     4
#define SCALE_BITS                      15
#define MAX_NUM_ALF_ALTERNATIVES_CHROMA 8
#define NUM_APS_TYPE_LEN                0 //1 //3

static const int g_fixed_filter_set_coeff[ALF_FIXED_FILTER_NUM][MAX_NUM_ALF_LUMA_COEFF] =
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

static const int alf_pattern_5[13] = {
              0,
          1,  2,  3,
      4,  5,  6,  5,  4,
          3,  2,  1,
              0
};
static const int alf_weights_5[8] = {
              2,
          2,  2,  2,
      2,  2,  1,  1
};

/*#if !JVET_O0216_ALF_COEFF_EG3 || !JVET_O0064_SIMP_ALF_CLIP_CODING
static const int alf_golomb_idx_5[8] = {
              0,
          0,  1,  0,
      0,  1,  2,  2
};
#endif*/

static const int alf_pattern_7[25] = {
              0,
          1,  2,  3,
      4,  5,  6,  7,  8,
  9, 10, 11, 12, 11, 10, 9,
      8,  7,  6,  5,  4,
          3,  2,  1,
              0
};
static const int alf_weights_7[14] = {
              2,
          2,  2,  2,
      2,  2,  2,  2,  2,
  2,  2,  2,  1,  1
};

/*#if !JVET_O0216_ALF_COEFF_EG3 || !JVET_O0064_SIMP_ALF_CLIP_CODING
static const int alf_golomb_idx_7[14] = {
              0,
          0,  1,  0,
      0,  1,  2,  1,  0,
  0,  1,  2,  3,  3
};
#endif*/

//-------------------------typedef enums----------------------------
typedef enum { ALF_FILTER_5X5 = 0, ALF_FILTER_7X7 = 1, ALF_NUM_OF_FILTER_TYPES = 2 } alf_filter_type;
typedef enum { ALF_LUMA = 0, ALF_CHROMA = 1 } alf_type;

typedef enum {
  T_ALF_APS = 0,
  T_LMCS_APS = 1,
} aps_type_values;

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
typedef struct alf_covariance {
  int num_coeff;
  int num_bins;
  double y[MAX_ALF_NUM_CLIPPING_VALUES][MAX_NUM_ALF_LUMA_COEFF];
  double ee[MAX_ALF_NUM_CLIPPING_VALUES][MAX_ALF_NUM_CLIPPING_VALUES][MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF];
  double pix_acc;
} alf_covariance;

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

clp_rngs g_clp_rngs;

typedef struct alf_classifier {
  int class_idx;
  int transpose_idx;
} alf_classifier;

typedef struct alf_info_t {
  alf_covariance*** g_alf_covariance[MAX_NUM_COMPONENT]; //[component_id][filter_type][ctu_idx][class_idx]
  alf_covariance** g_alf_covariance_frame[MAX_NUM_CHANNEL_TYPE]; //[channel][filter_type][class_idx]
  alf_covariance g_alf_covariance_merged[ALF_NUM_OF_FILTER_TYPES][MAX_NUM_ALF_CLASSES + 2];
  int g_alf_clip_merged[ALF_NUM_OF_FILTER_TYPES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF];
  uint8_t* g_ctu_enable_flag[MAX_NUM_COMPONENT];
  uint8_t* g_ctu_alternative[MAX_NUM_COMPONENT]; //#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  double *g_ctb_distortion_unfilter[MAX_NUM_COMPONENT];
  int g_aps_id_start;
  int** g_diff_filter_coeff; // [lumaClassIdx][coeffIdx]
  int** g_filter_coeff_set;  // [lumaClassIdx][coeffIdx]
  int** g_filter_clipp_set; // [lumaClassIdx][coeffIdx]
  short* g_alf_ctb_filter_set_index_tmp; //g_num_ctus_in_pic //voisi olla lokaali muuttuja?
  short* g_alf_ctb_filter_index;     //g_num_ctus_in_pic
} alf_info_t;

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
  int aps_type;

  //sliceparams
  bool enabled_flag[MAX_NUM_COMPONENT];                           // alf_slice_enable_flag, alf_chroma_idc
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  bool non_linear_flag[MAX_NUM_CHANNEL_TYPE][MAX_NUM_ALF_ALTERNATIVES_CHROMA]; // alf_[luma/chroma]_clip_flag
/*#else
  bool non_linear_flag[MAX_NUM_CHANNEL_TYPE];                       // alf_nonlinear_enable_flag[Luma/Chroma]
#endif*/
  short luma_coeff[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF]; // alf_coeff_luma_delta[i][j]
  short luma_clipp[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF]; // alf_clipp_luma_[i][j]
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  int num_alternatives_chroma;                                                  // alf_chroma_num_alts_minus_one + 1
  short chroma_coeff[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF]; // alf_coeff_chroma[i]
  short chroma_clipp[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF]; // alf_clipp_chroma[i]
/*#else
  short chroma_coeff[MAX_NUM_ALF_CHROMA_COEFF];                   // alf_coeff_chroma[i]
  short chroma_clipp[MAX_NUM_ALF_CHROMA_COEFF];                   // alf_clipp_chroma[i]
#endif*/
  short filter_coeff_delta_idx[MAX_NUM_ALF_CLASSES];              // filter_coeff_delta[i]
  bool alf_luma_coeff_flag[MAX_NUM_ALF_CLASSES];                  // alf_luma_coeff_flag[i]
  int num_luma_filters;                                           // number_of_filters_minus1 + 1
  bool alf_luma_coeff_delta_flag;                                 // alf_luma_coeff_delta_flag
/*#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
  bool alf_luma_coeff_delta_prediction_flag;                      // alf_luma_coeff_delta_prediction_flag
  int fixed_filter_pattern;
  int fixed_filter_idx[MAX_NUM_ALF_CLASSES];
  int fixed_filter_set_index;
#endif*/
  //std::vector<AlfFilterShape>* filterShapes;
  int t_layer;
  bool new_filter_flag[MAX_NUM_CHANNEL_TYPE];

  struct cc_alf_filter_param cc_alf_aps_param;

} alf_aps;

typedef struct param_set_map {
  bool b_changed;
  //uint8_t* p_nalu_data;
  struct alf_aps parameter_set;
} param_set_map;
//---------------------------------------------------------------

//dunno
uint8_t *g_alf_ctu_enable_flag[MAX_NUM_COMPONENT];
uint8_t* g_alf_ctu_alternative[MAX_NUM_COMPONENT];
alf_covariance*** g_alf_covariance_cc_alf[2]; // [compIdx-1][shapeIdx][ctbAddr][filterIdx]
alf_covariance** g_alf_covariance_frame_cc_alf[2]; // [compIdx-1][shapeIdx][filterIdx]
uint8_t*               g_training_cov_control;
uint64_t**             g_unfiltered_distortion;  // for different block size
uint64_t*              g_training_distortion[MAX_NUM_CC_ALF_FILTERS];    // for current block size
uint8_t*               g_filter_control;         // current iterations filter control
uint8_t*               g_best_filter_control;     // best saved filter control
uint64_t*              g_luma_swing_greater_than_threshold_count;
uint64_t*              g_chroma_sample_count_near_mid_point;
//tarpeeton jos WSSD=0
double* g_luma_level_to_weight_plut; //Ei anneta arvoja miss‰‰n

//defaults / consts
static const int g_input_bit_depth[2] = { 8, 8 };
static const int g_alf_num_clipping_values[MAX_NUM_CHANNEL_TYPE] = { 4,4 };
static const int g_max_alf_num_clipping_values = 4;
static const double frac_bits_scale = 1.0 / (double)(1 << SCALE_BITS);
static unsigned g_bits_new_filter[MAX_NUM_CHANNEL_TYPE];
static int g_clip_default_enc[MAX_NUM_ALF_LUMA_COEFF];

//define during run
//static int **g_laplacian[NUM_DIRECTIONS]; // korvattu alemmalla versiolla
static int g_laplacian[NUM_DIRECTIONS][CLASSIFICATION_BLK_SIZE + 5][CLASSIFICATION_BLK_SIZE + 5]; //puhdista vasta koko encodauksen lopuksi
static double g_lambda[MAX_NUM_COMPONENT];
short g_clip_default[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
static short g_filter_indices[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES];
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
short g_chroma_coeff_final[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF];
//alf_aps* g_alf_aps_chroma; //pelkk‰ turha v‰lik‰si (?)
/*#else
static short g_chroma_coeff_final[MAX_NUM_ALF_LUMA_COEFF]; #endif*/
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
short g_chroma_clipp_final[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF];
/*#else
static short g_chroma_clipp_final[MAX_NUM_ALF_LUMA_COEFF];
#endif*/
static short g_coeff_final[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
static int16_t g_clipp_final[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
static short g_coeff_aps_luma[ALF_CTB_MAX_NUM_APS][MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
static int16_t g_clipp_aps_luma[ALF_CTB_MAX_NUM_APS][MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
static short g_fixed_filter_set_coeff_dec[ALF_NUM_FIXED_FILTER_SETS][MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];

//once ever
static int kvz_bit_depth;
static int chroma_scale_x;
static int chroma_scale_y;
static int g_num_ctus_in_pic;
static int g_alf_vb_luma_ctu_height;
static int g_alf_vb_chma_ctu_height;
static int g_alf_vb_luma_pos;
static int g_alf_vb_chma_pos;
static short g_alf_clipping_values[MAX_NUM_CHANNEL_TYPE][MAX_ALF_NUM_CLIPPING_VALUES];
static alf_classifier **g_classifier;
static bool g_created = false;
static uint32_t g_frame_count = MAX_INT;
uint8_t* g_cc_alf_filter_control[2];
int g_aps_id_cc_alf_start[2];

//once per frame
alf_covariance*** g_alf_covariance[MAX_NUM_COMPONENT]; //[component_id][filter_type][ctu_idx][class_idx]
alf_covariance** g_alf_covariance_frame[MAX_NUM_CHANNEL_TYPE]; //[channel][filter_type][class_idx]
alf_covariance g_alf_covariance_merged[ALF_NUM_OF_FILTER_TYPES][MAX_NUM_ALF_CLASSES + 2];
int g_alf_clip_merged[ALF_NUM_OF_FILTER_TYPES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF];
uint8_t* g_ctu_enable_flag[MAX_NUM_COMPONENT];
uint8_t* g_ctu_alternative[MAX_NUM_COMPONENT]; //#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
double *g_ctb_distortion_unfilter[MAX_NUM_COMPONENT];
static int g_aps_id_start = ALF_CTB_MAX_NUM_APS;
int** g_diff_filter_coeff; // [lumaClassIdx][coeffIdx]
int** g_filter_coeff_set;  // [lumaClassIdx][coeffIdx]
int** g_filter_clipp_set; // [lumaClassIdx][coeffIdx]
short* g_alf_ctb_filter_set_index_tmp; //g_num_ctus_in_pic //voisi olla lokaali muuttuja?
short* g_alf_ctb_filter_index;     //g_num_ctus_in_pic
static uint32_t g_curr_frame = MAX_INT;
static uint32_t g_old_frame = 0;
#if !FULL_FRAME
alf_aps alf_param;
#endif
struct cc_alf_filter_param g_cc_alf_filter_param;

//temps
static alf_aps g_alf_aps_temp;
static alf_aps g_alf_aps_temp_nl;
uint8_t* g_ctu_enable_flag_tmp[MAX_NUM_COMPONENT]; //kerran
//uint8_t* g_ctu_enable_flag_tmp2[MAX_NUM_COMPONENT];
uint8_t* g_ctu_alternative_tmp[MAX_NUM_COMPONENT]; //kerran
static int g_filter_tmp[MAX_NUM_ALF_LUMA_COEFF];
static int g_clip_tmp[MAX_NUM_ALF_LUMA_COEFF];
//kvz_picture *tmp_rec_pic;

kvz_pixel *alf_fulldata;
kvz_pixel *alf_tmp_y;
kvz_pixel *alf_tmp_u;
kvz_pixel *alf_tmp_v;

//cabac temps
cabac_data_t ctx_start;
cabac_data_t ctx_start_cc_alf;
cabac_data_t cabac_estimator;

// encoder

//kvz_alf_encoder_ctb
//int best_aps_ids[ALF_CTB_MAX_NUM_APS]; //frame
//int size_of_best_aps_ids;
#if !FULL_FRAME
int new_aps_id;
int aps_ids[ALF_CTB_MAX_NUM_APS];
double d_dist_org_new_filter;
int blocks_using_new_filter;
#endif // !FULL_FRAME

int size_of_aps_ids;

//-------------------------init function----------------------------
#if !FULL_FRAME
void kvz_alf_init(encoder_state_t *const state);
  //alf_info_t *alf);
#endif

//------------------------------------------------------------------

//-------------------------help functions---------------------------

void set_config(kvz_config *const cfg);
bool is_crossed_by_virtual_boundaries(const int x_pos, const int y_pos, const int width, const int height, bool* clip_top, bool* clip_bottom, bool* clip_left, 
                                      bool* clip_right, int* num_hor_vir_bndry, int* num_ver_vir_bndry, int hor_vir_bndry_pos[], int ver_vir_bndry_pos[], encoder_state_t *const state);
void init_ctu_alternative_chroma(uint8_t* ctu_alts[MAX_NUM_COMPONENT]);
int16_t clip_alf(const int16_t clip, const int16_t ref, const int16_t val0, const int16_t val1);
int alf_clip_pixel(const int a, const clp_rng clp_rng);
int16_t alf_clip3(const int16_t min_val, const int16_t max_val, const int16_t a);
void get_clip_max(alf_covariance *cov, int *clip_max);
void reduce_clip_cost(alf_covariance *cov, int *clip);
void set_ey_from_clip(alf_covariance *cov, const int* clip, double ee[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double y[MAX_NUM_ALF_LUMA_COEFF], int siz);
double optimize_filter(alf_covariance *cov, int* clip, double *f, bool optimize_clip);
double optimize_filter_clip(alf_covariance *cov, int* clip);
double optimize_filter_gns_calc(alf_covariance *cov, const int* clip, double *f, int size);
void gns_backsubstitution(double r[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double* z, int size, double* A);
void gns_transpose_backsubstitution(double u[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double* rhs, double* x, int order);
int gns_cholesky_dec(double inp_matr[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double out_matr[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], int num_eq);
int gns_solve_by_chol(double lhs[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double rhs[MAX_NUM_ALF_LUMA_COEFF], double *x, int num_eq);
int gns_solve_by_chol_clip_gns(alf_covariance *cov, const int *clip, double *x, int num_eq);
double calc_error_for_coeffs(alf_covariance *cov, const int *clip, const int *coeff, const int num_coeff, const int bit_depth);
//int get_golomb_k_min(channel_type channel, const int numFilters, int kMinTab[MAX_NUM_ALF_LUMA_COEFF], int bitsCoeffScan[11/*m_MAX_SCAN_VAL*/][16/*m_MAX_EXP_GOLOMB*/]);
int length_golomb(int coeff_val, int k, bool signed_coeff);
double get_dist_coeff_force_0(bool* coded_var_bins, double error_force_0_coeff_tab[MAX_NUM_ALF_CLASSES][2], int* bits_var_bin, const int num_filters);
double get_dist_force_0(channel_type channel, const int num_filters, double error_tab_force_0_coeff[MAX_NUM_ALF_CLASSES][2], bool* coded_var_bins);
int get_cost_filter_coeff_force_0(channel_type channel, int **p_diff_q_filter_coeff_int_pp, const int num_filters, bool* coded_var_bins);
int get_cost_filter_coeff(channel_type channel, int **p_diff_q_filter_coeff_int_pp, const int num_filters);
int get_cost_filter_clipp(channel_type channel, int **p_diff_q_filter_coeff_int_pp, const int num_filters);
//int get_tb_length(int ui_symbol, const int ui_max_symbol);//#if !JVET_O0491_HLS_CLEANUP
int get_non_filter_coeff_rate(alf_aps *aps);
int length_filter_coeffs(channel_type channel, const int num_filters, int **filter_coeff);
double calculate_error(alf_covariance *cov, const int *clip, const double *coeff);
double calculate_error_opt_filt(alf_covariance *cov, const int *clip);
//int get_coeff_rate(alf_aps *aps, bool is_chroma);//#if !JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
int get_chroma_coeff_rate(alf_aps* aps, int alt_idx); 
//int length_truncated_unary(int symbol, int max_symbol);//#if !JVET_O0491_HLS_CLEANUP
double get_filtered_distortion(alf_covariance* cov, const int num_classes, const int num_filters_minus1, const int num_coeff);
double get_unfiltered_distortion_cov_channel(alf_covariance* cov, channel_type channel);
double get_unfiltered_distortion_cov_classes(alf_covariance* cov, const int num_classes);
void get_frame_stats(channel_type channel, int i_shape_idx
#if !FULL_FRAME
  , int ctu_idx
#endif
  );
void get_frame_stat(alf_covariance* frame_cov, alf_covariance** ctb_cov, uint8_t* ctb_enable_flags, uint8_t* ctb_alt_idx, const int num_classes, int alt_idx
#if !FULL_FRAME
  , int ctu_idx
#endif
  );

void copy_cov(alf_covariance *dst, alf_covariance *src);
void copy_alf_param(alf_aps *dst, alf_aps *src);
void copy_alf_param_w_channel(alf_aps* dst, alf_aps* src, channel_type channel);
void copy_aps(alf_aps *dst, alf_aps *src);
void copy_aps_to_map(param_set_map *dst, alf_aps *src, int8_t aps_id);
//bool compare_alf_param(const alf_aps* aps_1, const alf_aps* aps_2);
void reset_alf_param(alf_aps *src);
void add_alf_cov(alf_covariance *dst, alf_covariance *src);
void add_alf_cov_lhs_rhs(alf_covariance *dst, alf_covariance *lhs, alf_covariance *rhs);
void reset_alf_covariance(alf_covariance *alf, int num_bins);
void reset_cc_alf_aps_param(cc_alf_filter_param *cc_alf);
void copy_pixels(kvz_pixel *src, int x_src_start, int y_src_start, int src_stride,
  kvz_pixel *dst, int x_dst_start, int y_dst_start, int dst_stride,
  int width, int height);
void adjust_pixels(kvz_pixel *src, int x_start, int x_end, int y_start, int y_end,
  int stride, int pic_width, int pic_height);
void adjust_pixels_CTU_plus_4_pix(kvz_pixel *src, int x_start, int x_end, int y_start, int y_end,
                   int stride, int pic_width, int pic_height);
//Need to adjust
void adjust_pixels_chroma(kvz_pixel *src, int x_start, int x_end, int y_start, int y_end, 
                  int stride, int pic_width, int pic_height);

void set_ctu_enable_flag(uint8_t **flags, channel_type channel, 
#if !FULL_FRAME
  int ctu_idx,
#endif
  uint8_t value);
void copy_ctu_enable_flag(uint8_t **flags_dst, uint8_t **flags_src, channel_type channel
#if !FULL_FRAME
  , int ctu_idx
#endif
  );

//-------------------------------------------------------------------

//-------------------------encoding functions------------------------

void apply_cc_alf_filter(encoder_state_t *const state, alf_component_id comp_id, const kvz_pixel *dst_pixels,
  const kvz_pixel *recYuvExt, uint8_t *filterControl,
  const short filterSet[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF],
  const int   selectedFilterIdx
  );

//is_crossed_by_virtual_boundaries -osuus ep‰t‰ydellinen
void kvz_alf_enc_process(encoder_state_t *const state
#if !FULL_FRAME
  , const lcu_order_element_t *const lcu
#endif
  );

double kvz_alf_derive_ctb_alf_enable_flags(encoder_state_t * const state,
  channel_type channel,
  const int i_shape_idx,
  double *dist_unfilter,
  const int num_classes,
#if !FULL_FRAME
  int ctu_idx,
#endif
  const double chroma_weight
  );

void kvz_alf_enc_create(encoder_state_t const *state);

void kvz_alf_reconstruct(encoder_state_t const *state
#if !FULL_FRAME
  , const lcu_order_element_t *const lcu
#endif
  );

void kvz_alf_enc_destroy(videoframe_t * const frame);

void kvz_alf_encoder(encoder_state_t *const state,
#if !FULL_FRAME
  const lcu_order_element_t *lcu,
#endif
  alf_aps *aps,
  channel_type channel,
  const double lambda_chroma_weight
  );

//isIntra, PendingRasInit, IDRorBLA <--- ? selvit‰ n‰m‰
void kvz_alf_get_avai_aps_ids_luma(encoder_state_t *const state,
  int *newApsId,
  int *aps_ids,
  int *size_of_aps_ids);

void kvz_alf_derive_stats_for_filtering(encoder_state_t *const state
#if !FULL_FRAME
  ,  const lcu_order_element_t *const lcu
#endif
  );

//mik‰ on alf_WSSD?
void kvz_alf_get_blk_stats(encoder_state_t *const state,
  channel_type channel,
  alf_covariance **alfCovariace,
  alf_classifier **g_classifier,
  kvz_pixel *org,
  int32_t org_stride,
  kvz_pixel *rec,
  int32_t rec_stride,
  const int x_pos,
  const int y_pos,
  const int x_dst,
  const int y_dst,
  const int width,
  const int height,
  int vb_ctu_height,
  int vb_pos);

void kvz_alf_calc_covariance(int16_t e_local[MAX_NUM_ALF_LUMA_COEFF][MAX_ALF_NUM_CLIPPING_VALUES],
  const kvz_pixel *rec,
  const int stride,
  const channel_type channel,
  const int transpose_idx,
  int vb_distance);

double kvz_alf_get_filter_coeff_and_cost(encoder_state_t *const state,
  channel_type channel,
  double dist_unfilter,
  int *ui_coeff_bits,
  int i_shape_idx,
  bool b_re_collect_stat,
  bool only_filter_cost
#if !FULL_FRAME
  , int ctu_idx
#endif
  );

int kvz_alf_derive_filter_coefficients_prediction_mode(channel_type channel,
  int **filter_set,
  int **filter_coeff_diff,
  const int num_filters);

void kvz_alf_merge_classes(channel_type channel,
  alf_covariance* cov,
  alf_covariance* cov_merged,
  int clip_merged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF],
  const int num_classes,
  short filter_indices[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES]);

double kvz_alf_merge_filters_and_cost(encoder_state_t *const state,
  alf_aps *alf_aps,
  channel_type channel,
  int *ui_coeff_bits,
  alf_covariance *cov_frame,
  alf_covariance *cov_merged,
  int clip_merged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF]);

double kvz_alf_derive_filter_coeffs(alf_aps *aps,
  channel_type channel,
  alf_covariance *cov,
  alf_covariance *covMerged,
  short* filter_indices,
  int numFilters,
  double errorTabForce0Coeff[MAX_NUM_ALF_CLASSES][2],
  int clip_merged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF]);

double kvz_alf_derive_coeff_quant(channel_type channel,
  int *filter_clipp,
  int *filter_coeff_quant,
  const alf_covariance* cov,
  const int bit_depth,
  const bool oprimize_clip);

//Muutama asia pit‰‰ viel‰ p‰ivitt‰‰
//bookmarks
void kvz_alf_encoder_ctb(encoder_state_t *const state,
  alf_aps *aps,
#if !FULL_FRAME
  int ctu_idx,
#endif
  const double lambda_chroma_weight
  );

void kvz_alf_reconstructor(encoder_state_t const *state);

/*
void apply_cc_alf_filter(encoder_state_t *const state, alf_component_id comp_id, const kvz_pixel *dst_pixels,
  const kvz_pixel *recYuvExt, uint8_t *filterControl,
  const short filterSet[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF],
  const int   selectedFilterIdx)
*/
//----------------------------------------------------------------------

//-------------------------cabac writer functions------------------------

void kvz_cabac_reset_bits(cabac_data_t * const data);

void code_alf_ctu_enable_flags_channel(encoder_state_t * const state,
  cabac_data_t * const cabac,
  channel_type channel,
  alf_aps *aps);

void code_alf_ctu_enable_flags_component(encoder_state_t * const state,
  cabac_data_t * const cabac,
  alf_component_id component_id,
  alf_aps *aps);

void code_alf_ctu_enable_flag(encoder_state_t * const state,
  cabac_data_t * const cabac,
  uint32_t ctu_rs_addr,
  alf_component_id component_id,
  alf_aps *aps);

void code_alf_ctu_filter_index(encoder_state_t * const state,
  cabac_data_t * const cabac,
  uint32_t ctu_rs_addr,
  bool alf_enable_luma);

void code_alf_ctu_alternatives_channel(encoder_state_t * const state,
  cabac_data_t * const cabac,
  channel_type channel,
  alf_aps* aps
#if !FULL_FRAME
  , int ctu_idx
#endif
  );

void code_alf_ctu_alternatives_component(encoder_state_t * const state,
  cabac_data_t * const cabac,
  alf_component_id comp_id,
  alf_aps* aps
#if !FULL_FRAME
  , int ctu_idx
#endif  
  );

void code_alf_ctu_alternative_ctu(encoder_state_t * const state,
  cabac_data_t * const cabac,
  uint32_t ctu_rs_addr,
  const alf_component_id comp_idx,
  const alf_aps* aps);

void kvz_encode_alf_bits(encoder_state_t * const state, const int ctu_idx);

void encoder_state_write_adaptation_parameter_set(encoder_state_t * const state,
  alf_aps *aps);

void encode_alf_aps_flags(encoder_state_t * const state,
  alf_aps* aps);

void encode_alf_aps_filter(encoder_state_t * const state,
  alf_aps* aps,
  const bool is_chroma,
  const int alt_idx);

/*
void alf_golomb_encode(encoder_state_t * const state,
  int coeff,
  int k,
  const bool signed_coeff);
  */
void encode_alf_adaptive_parameter_set(encoder_state_t * const state);

void encode_alf_aps_lmcs(encoder_state_t * const state);

void encode_alf_aps_scaling_list(encoder_state_t * const state);

void encode_alf_aps(encoder_state_t * const state);

//---------------------------------------------------------------------

//-------------------------CTU functions--------------------------------

//ei varmuutta miten alf_param_tmp pit‰isi toimia t‰ss‰ tilanteessa
void kvz_alf_reconstruct_coeff_aps(encoder_state_t *const state,
  bool luma,
  bool chroma,
  bool is_rdo);

void kvz_alf_reconstruct_coeff(encoder_state_t *const state,
  alf_aps *aps,
  channel_type channel,
  const bool is_rdo,
  const bool is_redo);

void kvz_alf_create(encoder_state_t const *state);

void kvz_alf_destroy(videoframe_t * const frame);

void kvz_alf_derive_classification(encoder_state_t *const state,
  const int width,
  const int height,
  int x_pos,
  int y_pos,
  const int blk_dst_x,
  const int blk_dst_y);

/*#if !JVET_O0525_REMOVE_PCM
//VTM6.0 (ei muuttunut)
//OK
//Turha jos PCM on pois p‰‰lt‰.
//cu->ipcm?
void kvz_alf_reset_pcm_blk_class_info(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  const int width,
  const int height,
  int x_pos,
  int y_pos);*/

void kvz_alf_derive_classification_blk(encoder_state_t *const state,
  const int shift,
  const int n_height,
  const int n_width,
  const int blk_pos_x,
  const int blk_pos_y,
  const int blk_dst_x,
  const int blk_dst_y,
  const int vb_ctu_height,
  int vb_pos);

//Pit‰isi testata; x,y,pixels,toimivuus
//kvz_alf_reconstructor-funktiossa VTM:ss‰ return, pois jotta
//p‰‰see t‰h‰n funktioon
//tarkista toimivuus
void kvz_alf_filter_block(encoder_state_t *const state,
  const kvz_pixel *src_pixels,
  kvz_pixel *dst_pixels,
  const int src_stride,
  const int dst_stride,
  const short* filter_set,
  const int16_t *fClipSet,
  clp_rng clp_rng,
  alf_component_id component_id,
  const int width,
  const int height,
  int x_pos,
  int y_pos,
  int blk_dst_x,
  int blk_dst_y,
  int vb_pos,
  const int vb_ctu_height);
//---------------------------------------------------------------------

//kvz_alf_filter_block funktion dst?
#endif //ALF_H_