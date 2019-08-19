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


//ALF applied to 4x4 luma samples
//Filtering in CTB level
//Signalled in slice header if used or not
//Signalled in CTB level if ALF used for luma (if applied, then also for both chroma samples)

//To reduce bits overhead, filter coefficients of different classification can be merged.
//In slice header, the indices of the APSs used for the current slice are signaled.

#define MAX_NUM_ALF_CLASSES             25
#define MAX_NUM_ALF_LUMA_COEFF          13
#define MAX_NUM_ALF_CHROMA_COEFF        7
#define MAX_ALF_FILTER_LENGTH           7
#define MAX_NUM_ALF_COEFF               (MAX_ALF_FILTER_LENGTH * MAX_ALF_FILTER_LENGTH / 2 + 1)
#define ALF_NUM_FIXED_FILTER_SETS       16
#define ALF_NUM_BITS                    8
#define ALF_UNUSED_CLASS_IDX            255
#define ALF_UNUSED_TRANSPOSE_IDX        255
#define REG                             0.0001
#define REG_SQR                         0.0000001
#define MAX_NUM_APS                     32
#define CLASSIFICATION_BLK_SIZE         32

/* #if JVET_N0242_NON_LINEAR_ALF
static const int g_fixed_filter_set_coeff[64][13] =
{
  { 0,   0,   2,  -3,   1,  -4,   1,   7,  -1,   1,  -1,   5, 0 },
  { 0,   0,   0,   0,   0,  -1,   0,   1,   0,   0,  -1,   2, 0 },
  { 0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0, 0 },
  { 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   1, 0 },
  { 2,   2,  -7,  -3,   0,  -5,  13,  22,  12,  -3,  -3,  17,  0 },
  { -1,   0,   6,  -8,   1,  -5,   1,  23,   0,   2,  -5,  10,  0 },
  { 0,   0,  -1,  -1,   0,  -1,   2,   1,   0,   0,  -1,   4, 0 },
  { 0,   0,   3, -11,   1,   0,  -1,  35,   5,   2,  -9,   9,  0 },
  { 0,   0,   8,  -8,  -2,  -7,   4,   4,   2,   1,  -1,  25,  0 },
  { 0,   0,   1,  -1,   0,  -3,   1,   3,  -1,   1,  -1,   3, 0 },
  { 0,   0,   3,  -3,   0,  -6,   5,  -1,   2,   1,  -4,  21,  0 },
  { -7,   1,   5,   4,  -3,   5,  11,  13,  12,  -8,  11,  12,  0 },
  { -5,  -3,   6,  -2,  -3,   8,  14,  15,   2,  -7,  11,  16,  0 },
  { 2,  -1,  -6,  -5,  -2,  -2,  20,  14,  -4,   0,  -3,  25,  0 },
  { 3,   1,  -8,  -4,   0,  -8,  22,   5,  -3,   2, -10,  29,  0 },
  { 2,   1,  -7,  -1,   2, -11,  23,  -5,   0,   2, -10,  29,  0 },
  { -6,  -3,   8,   9,  -4,   8,   9,   7,  14,  -2,   8,   9,  0 },
  { 2,   1,  -4,  -7,   0,  -8,  17,  22,   1,  -1,  -4,  23,  0 },
  { 3,   0,  -5,  -7,   0,  -7,  15,  18,  -5,   0,  -5,  27,  0 },
  { 2,   0,   0,  -7,   1, -10,  13,  13,  -4,   2,  -7,  24,  0 },
  { 3,   3, -13,   4,  -2,  -5,   9,  21,  25,  -2,  -3,  12,  0 },
  { -5,  -2,   7,  -3,  -7,   9,   8,   9,  16,  -2,  15,  12,  0 },
  { 0,  -1,   0,  -7,  -5,   4,  11,  11,   8,  -6,  12,  21,  0 },
  { 3,  -2,  -3,  -8,  -4,  -1,  16,  15,  -2,  -3,   3,  26,  0 },
  { 2,   1,  -5,  -4,  -1,  -8,  16,   4,  -2,   1,  -7,  33,  0 },
  { 2,   1,  -4,  -2,   1, -10,  17,  -2,   0,   2, -11,  33,  0 },
  { 1,  -2,   7, -15, -16,  10,   8,   8,  20,  11,  14,  11,  0 },
  { 2,   2,   3, -13, -13,   4,   8,  12,   2,  -3,  16,  24,  0 },
  { 1,   4,   0,  -7,  -8,  -4,   9,   9,  -2,  -2,   8,  29,  0 },
  { 1,   1,   2,  -4,  -1,  -6,   6,   3,  -1,  -1,  -3,  30,  0 },
  { -7,   3,   2,  10,  -2,   3,   7,  11,  19,  -7,   8,  10, 0 },
  { 0,  -2,  -5,  -3,  -2,   4,  20,  15,  -1,  -3,  -1,  22,  0 },
  { 3,  -1,  -8,  -4,  -1,  -4,  22,   8,  -4,   2,  -8,  28,  0 },
  { 0,   3, -14,   3,   0,   1,  19,  17,   8,  -3,  -7,  20,  0 },
  { 0,   2,  -1,  -8,   3,  -6,   5,  21,   1,   1,  -9,  13,  0 },
  { -4,  -2,   8,  20,  -2,   2,   3,   5,  21,   4,   6,   1, 0 },
  { 2,  -2,  -3,  -9,  -4,   2,  14,  16,   3,  -6,   8,  24,  0 },
  { 2,   1,   5, -16,  -7,   2,   3,  11,  15,  -3,  11,  22,  0 },
  { 1,   2,   3, -11,  -2,  -5,   4,   8,   9,  -3,  -2,  26,  0 },
  { 0,  -1,  10,  -9,  -1,  -8,   2,   3,   4,   0,   0,  29,  0 },
  { 1,   2,   0,  -5,   1,  -9,   9,   3,   0,   1,  -7,  20,  0 },
  { -2,   8,  -6,  -4,   3,  -9,  -8,  45,  14,   2, -13,   7, 0 },
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
*/
static const int g_fixed_filter_set_coeff[64][13] =
{
  { 0,   0,   2,  -3,   1,  -4,   1,   7,  -1,   1,  -1,   5, 112 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   0,   0,   0,   0,  -1,   0,   1,   0,   0,  -1,   2, 126 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0, 126 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   1, 128 - (1 << (ALF_NUM_BITS - 1)) },
  { 2,   2,  -7,  -3,   0,  -5,  13,  22,  12,  -3,  -3,  17,  34 - (1 << (ALF_NUM_BITS - 1)) },
  { -1,   0,   6,  -8,   1,  -5,   1,  23,   0,   2,  -5,  10,  80 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   0,  -1,  -1,   0,  -1,   2,   1,   0,   0,  -1,   4, 122 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   0,   3, -11,   1,   0,  -1,  35,   5,   2,  -9,   9,  60 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   0,   8,  -8,  -2,  -7,   4,   4,   2,   1,  -1,  25,  76 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   0,   1,  -1,   0,  -3,   1,   3,  -1,   1,  -1,   3, 122 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   0,   3,  -3,   0,  -6,   5,  -1,   2,   1,  -4,  21,  92 - (1 << (ALF_NUM_BITS - 1)) },
  { -7,   1,   5,   4,  -3,   5,  11,  13,  12,  -8,  11,  12,  16 - (1 << (ALF_NUM_BITS - 1)) },
  { -5,  -3,   6,  -2,  -3,   8,  14,  15,   2,  -7,  11,  16,  24 - (1 << (ALF_NUM_BITS - 1)) },
  { 2,  -1,  -6,  -5,  -2,  -2,  20,  14,  -4,   0,  -3,  25,  52 - (1 << (ALF_NUM_BITS - 1)) },
  { 3,   1,  -8,  -4,   0,  -8,  22,   5,  -3,   2, -10,  29,  70 - (1 << (ALF_NUM_BITS - 1)) },
  { 2,   1,  -7,  -1,   2, -11,  23,  -5,   0,   2, -10,  29,  78 - (1 << (ALF_NUM_BITS - 1)) },
  { -6,  -3,   8,   9,  -4,   8,   9,   7,  14,  -2,   8,   9,  14 - (1 << (ALF_NUM_BITS - 1)) },
  { 2,   1,  -4,  -7,   0,  -8,  17,  22,   1,  -1,  -4,  23,  44 - (1 << (ALF_NUM_BITS - 1)) },
  { 3,   0,  -5,  -7,   0,  -7,  15,  18,  -5,   0,  -5,  27,  60 - (1 << (ALF_NUM_BITS - 1)) },
  { 2,   0,   0,  -7,   1, -10,  13,  13,  -4,   2,  -7,  24,  74 - (1 << (ALF_NUM_BITS - 1)) },
  { 3,   3, -13,   4,  -2,  -5,   9,  21,  25,  -2,  -3,  12,  24 - (1 << (ALF_NUM_BITS - 1)) },
  { -5,  -2,   7,  -3,  -7,   9,   8,   9,  16,  -2,  15,  12,  14 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,  -1,   0,  -7,  -5,   4,  11,  11,   8,  -6,  12,  21,  32 - (1 << (ALF_NUM_BITS - 1)) },
  { 3,  -2,  -3,  -8,  -4,  -1,  16,  15,  -2,  -3,   3,  26,  48 - (1 << (ALF_NUM_BITS - 1)) },
  { 2,   1,  -5,  -4,  -1,  -8,  16,   4,  -2,   1,  -7,  33,  68 - (1 << (ALF_NUM_BITS - 1)) },
  { 2,   1,  -4,  -2,   1, -10,  17,  -2,   0,   2, -11,  33,  74 - (1 << (ALF_NUM_BITS - 1)) },
  { 1,  -2,   7, -15, -16,  10,   8,   8,  20,  11,  14,  11,  14 - (1 << (ALF_NUM_BITS - 1)) },
  { 2,   2,   3, -13, -13,   4,   8,  12,   2,  -3,  16,  24,  40 - (1 << (ALF_NUM_BITS - 1)) },
  { 1,   4,   0,  -7,  -8,  -4,   9,   9,  -2,  -2,   8,  29,  54 - (1 << (ALF_NUM_BITS - 1)) },
  { 1,   1,   2,  -4,  -1,  -6,   6,   3,  -1,  -1,  -3,  30,  74 - (1 << (ALF_NUM_BITS - 1)) },
  { -7,   3,   2,  10,  -2,   3,   7,  11,  19,  -7,   8,  10,  14 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,  -2,  -5,  -3,  -2,   4,  20,  15,  -1,  -3,  -1,  22,  40 - (1 << (ALF_NUM_BITS - 1)) },
  { 3,  -1,  -8,  -4,  -1,  -4,  22,   8,  -4,   2,  -8,  28,  62 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   3, -14,   3,   0,   1,  19,  17,   8,  -3,  -7,  20,  34 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   2,  -1,  -8,   3,  -6,   5,  21,   1,   1,  -9,  13,  84 - (1 << (ALF_NUM_BITS - 1)) },
  { -4,  -2,   8,  20,  -2,   2,   3,   5,  21,   4,   6,   1,   4 - (1 << (ALF_NUM_BITS - 1)) },
  { 2,  -2,  -3,  -9,  -4,   2,  14,  16,   3,  -6,   8,  24,  38 - (1 << (ALF_NUM_BITS - 1)) },
  { 2,   1,   5, -16,  -7,   2,   3,  11,  15,  -3,  11,  22,  36 - (1 << (ALF_NUM_BITS - 1)) },
  { 1,   2,   3, -11,  -2,  -5,   4,   8,   9,  -3,  -2,  26,  68 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,  -1,  10,  -9,  -1,  -8,   2,   3,   4,   0,   0,  29,  70 - (1 << (ALF_NUM_BITS - 1)) },
  { 1,   2,   0,  -5,   1,  -9,   9,   3,   0,   1,  -7,  20,  96 - (1 << (ALF_NUM_BITS - 1)) },
  { -2,   8,  -6,  -4,   3,  -9,  -8,  45,  14,   2, -13,   7,  54 - (1 << (ALF_NUM_BITS - 1)) },
  { 1,  -1,  16, -19,  -8,  -4,  -3,   2,  19,   0,   4,  30,  54 - (1 << (ALF_NUM_BITS - 1)) },
  { 1,   1,  -3,   0,   2, -11,  15,  -5,   1,   2,  -9,  24,  92 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   1,  -2,   0,   1,  -4,   4,   0,   0,   1,  -4,   7, 120 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   1,   2,  -5,   1,  -6,   4,  10,  -2,   1,  -4,  10, 104 - (1 << (ALF_NUM_BITS - 1)) },
  { 3,   0,  -3,  -6,  -2,  -6,  14,   8,  -1,  -1,  -3,  31,  60 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   1,   0,  -2,   1,  -6,   5,   1,   0,   1,  -5,  13, 110 - (1 << (ALF_NUM_BITS - 1)) },
  { 3,   1,   9, -19, -21,   9,   7,   6,  13,   5,  15,  21,  30 - (1 << (ALF_NUM_BITS - 1)) },
  { 2,   4,   3, -12, -13,   1,   7,   8,   3,   0,  12,  26,  46 - (1 << (ALF_NUM_BITS - 1)) },
  { 3,   1,  -8,  -2,   0,  -6,  18,   2,  -2,   3, -10,  23,  84 - (1 << (ALF_NUM_BITS - 1)) },
  { 1,   1,  -4,  -1,   1,  -5,   8,   1,  -1,   2,  -5,  10, 112 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   1,  -1,   0,   0,  -2,   2,   0,   0,   1,  -2,   3, 124 - (1 << (ALF_NUM_BITS - 1)) },
  { 1,   1,  -2,  -7,   1,  -7,  14,  18,   0,   0,  -7,  21,  62 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   1,   0,  -2,   0,  -7,   8,   1,  -2,   0,  -3,  24,  88 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   1,   1,  -2,   2, -10,  10,   0,  -2,   1,  -7,  23,  94 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   2,   2, -11,   2,  -4,  -3,  39,   7,   1, -10,   9,  60 - (1 << (ALF_NUM_BITS - 1)) },
  { 1,   0,  13, -16,  -5,  -6,  -1,   8,   6,   0,   6,  29,  58 - (1 << (ALF_NUM_BITS - 1)) },
  { 1,   3,   1,  -6,  -4,  -7,   9,   6,  -3,  -2,   3,  33,  60 - (1 << (ALF_NUM_BITS - 1)) },
  { 4,   0, -17,  -1,  -1,   5,  26,   8,  -2,   3, -15,  30,  48 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   1,  -2,   0,   2,  -8,  12,  -6,   1,   1,  -6,  16, 106 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   0,   0,  -1,   1,  -4,   4,   0,   0,   0,  -3,  11, 112 - (1 << (ALF_NUM_BITS - 1)) },
  { 0,   1,   2,  -8,   2,  -6,   5,  15,   0,   2,  -7,   9,  98 - (1 << (ALF_NUM_BITS - 1)) },
  { 1,  -1,  12, -15,  -7,  -2,   3,   6,   6,  -1,   7,  30,  50 - (1 << (ALF_NUM_BITS - 1)) },
};

static const int g_class_to_filter_mapping[16][25] =
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
static const int alf_golomb_idx[8] = {
              0,
          0,  1,  0,
      0,  1,  2,  2
};

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
static const int alf_golomb_idx_7[14] = {
              0,
          0,  1,  0,
      0,  1,  2,  1,  0,
  0,  1,  2,  3,  3
};

typedef enum { ALF_FILTER_5X5 = 0, ALF_FILTER_7X7 = 1, ALF_NUM_OF_FILTER_TYPES = 2 } alf_filter_type;
typedef enum { ALF_LUMA = 0, ALF_CHROMA = 1 } alf_type;

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

static short g_fixed_filter_set_coeff_dec[ALF_NUM_FIXED_FILTER_SETS][MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];

static short g_coeff_aps_luma[6][MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
static short g_coeff_final[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];

//short g_luma_coeff[max_num_alf_classes * max_num_alf_luma_coeff];
//short g_chroma_coeff[max_num_alf_chroma_coeff];
static short g_chroma_coeff_final[MAX_NUM_ALF_LUMA_COEFF];
//short g_filter_coeff_delta_idx[max_num_alf_classes];
static short g_filter_indices[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES];
static int g_filter_tmp[MAX_NUM_ALF_LUMA_COEFF];
//int g_fixed_filter_idx[max_num_alf_classes];
static int g_bits_coeff_scan[11/*m_MAX_SCAN_VAL*/][16/*m_MAX_EXP_GOLOMB*/];

//static short* g_alf_ctb_filter_set_index;
short* g_alf_ctb_filter_set_index_tmp;
short* g_alf_ctb_filter_index;
//int g_fixed_filter_set_index;
//int g_fixed_filter_pattern;
//int g_num_luma_filters;
//bool g_alf_luma_coeff_delta_flag;
//bool g_alf_luma_coeff_delta_prediction_flag;
//bool g_alf_luma_coeff_flag[max_num_alf_classes];

int** g_diff_filter_coeff;
int** g_filter_coeff_set;
int* g_filter_coeff_quant;

static int g_aps_id_start;
//static double* g_ctb_distortion_fixed_filter; //Ei k‰ytetty miss‰‰n

static unsigned g_bits_new_filter[MAX_NUM_CHANNEL_TYPE];
double *g_ctb_distortion_unfilter[MAX_NUM_COMPONENT];
static double g_lambda[MAX_NUM_COMPONENT];

static int g_k_min_tab[MAX_NUM_ALF_LUMA_COEFF];

static const double frac_bits_scale = 1.0 / (double)(1 << 15/*SCALE_BITS*/);

static int alf_input_bit_depth[2] = { 8, 8 };

static int num_ctus_in_pic;

typedef struct alf_covariance {
  int num_coeff;
  double *y;
  double **ee;
  double pix_acc;
} alf_covariance;

alf_covariance*** g_alf_covariance[MAX_NUM_COMPONENT]; //[component_id][filter_type][ctu_idx][class_idx]
alf_covariance** g_alf_covariance_frame[MAX_NUM_CHANNEL_TYPE]; //[channel][filter_type][class_idx]
alf_covariance g_alf_covariance_merged[ALF_NUM_OF_FILTER_TYPES][MAX_NUM_ALF_CLASSES + 2];
uint8_t* g_ctu_enable_flag_tmp[MAX_NUM_COMPONENT];
uint8_t* g_ctu_enable_flag[MAX_NUM_COMPONENT];
//uint8_t g_alf_slice_param_temp[MAX_NUM_COMPONENT];

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
  //alf_type blk_type;
  //alf_filter_type filter_type;
  alf_component_id component_id;
  //alf_directions direction;
  const kvz_pixel *src_buf;
  const kvz_pixel *dst_buf;
  alf_classifier **classifier;
  int **laplacian[NUM_DIRECTIONS];
  bool created;
  //clp_rngs clp_rngs;
} alf_info_t;

typedef struct alf_aps {
  int aps_id;

  //sliceparam
  bool enabled_flag[MAX_NUM_COMPONENT];                          // alf_slice_enable_flag, alf_chroma_idc
  short luma_coeff[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF]; // alf_coeff_luma_delta[i][j]
  short chroma_coeff[MAX_NUM_ALF_CHROMA_COEFF];                   // alf_coeff_chroma[i]
  short filter_coeff_delta_idx[MAX_NUM_ALF_CLASSES];                // filter_coeff_delta[i]
  bool alf_luma_coeff_flag[MAX_NUM_ALF_CLASSES];                   // alf_luma_coeff_flag[i]
  int num_luma_filters;                                          // number_of_filters_minus1 + 1
  bool alf_luma_coeff_delta_flag;                                   // alf_luma_coeff_delta_flag
  bool alf_luma_coeff_delta_prediction_flag;                         // alf_luma_coeff_delta_prediction_flag
  //std::vector<AlfFilterShape>* filterShapes;
  int t_layer;
  bool new_filter_flag[MAX_NUM_CHANNEL_TYPE];
  int fixed_filter_pattern;
  int fixed_filter_idx[MAX_NUM_ALF_CLASSES];
  int fixed_filter_set_index;
} alf_aps;

static alf_aps g_alf_slice_aps_temp;

typedef struct param_set_map {
  bool b_changed;
  //uint8_t* p_nalu_data;
  struct alf_aps parameter_set;
} param_set_map;

void kvz_alf_init(encoder_state_t *const state,
  encoder_state_config_slice_t *slice,
  alf_info_t *alf);

//-------------------------help functions---------------------------
int alf_clip_pixel(const int a, const clp_rng clp_rng);
int alf_clip3(const int minVal, const int maxVal, const int a);
void gns_backsubstitution(double R[MAX_NUM_ALF_COEFF][MAX_NUM_ALF_COEFF], double* z, int size, double* A);
void gns_transpose_backsubstitution(double U[MAX_NUM_ALF_COEFF][MAX_NUM_ALF_COEFF], double* rhs, double* x, int order);
int gns_cholesky_dec(double **inpMatr, double outMatr[MAX_NUM_ALF_COEFF][MAX_NUM_ALF_COEFF], int numEq);
int gns_solve_by_chol(double **LHS, double *rhs, double *x, int numEq);
double calc_error_for_coeffs(double **ee, double *y, const int *coeff, const int numCoeff, const int bitDepth);
int get_golomb_k_min(channel_type channel, const int numFilters, int kMinTab[MAX_NUM_ALF_LUMA_COEFF], int bitsCoeffScan[11/*m_MAX_SCAN_VAL*/][16/*m_MAX_EXP_GOLOMB*/]);
int length_golomb(int coeffVal, int k);
double get_dist_force_0(channel_type channel, const int numFilters, double errorTabForce0Coeff[MAX_NUM_ALF_CLASSES][2], bool* codedVarBins);
int get_cost_filter_coeff_force_0(channel_type channel, int **pDiffQFilterCoeffIntPP, const int numFilters, bool* codedVarBins);
int get_cost_filter_coeff(channel_type channel, int **pDiffQFilterCoeffIntPP, const int numFilters);
int get_tb_length(int uiSymbol, const int uiMaxSymbol);
int get_non_filter_coeff_rate(alf_aps *aps);
double calculate_error(alf_covariance cov);
int get_coeff_rate(alf_aps *aps, bool isChroma);
int length_truncated_unary(int symbol, int max_symbol);
double get_unfiltered_distortion_cov_channel(alf_covariance* cov, channel_type channel);
double get_unfiltered_distortion_cov_classes(alf_covariance* cov, const int numClasses);
void get_frame_stats(channel_type channel, int iShapeIdx, int m_ctuEnableFlag);
void get_frame_stat(alf_covariance* frameCov, alf_covariance** ctbCov, int ctbEnableFlags, const int numClasses);
//-------------------------------------------------------------------

//-------------------------encoding functions------------------------
//
void kvz_alf_enc_process(encoder_state_t *const state,
  const lcu_order_element_t *const lcu);

//OK
void kvz_alf_enc_create(encoder_state_t const *state,
  const lcu_order_element_t *lcu);

void kvz_alf_enc_destroy(encoder_state_t const *state);

//Pit‰isi olla OK
void kvz_alf_encoder(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  alf_aps *aps,
  channel_type channel);

void kvz_alf_encoder_ctb(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  alf_aps *aps);

//p‰ivit‰
void kvz_alf_reconstructor(encoder_state_t const *state,
  const lcu_order_element_t *const lcu);

//OK
//g_alf_covariance_frame arvot riippuvat g_alf_covariance: arvoista
//eli funktiosta kvz_alf_get_blk_stats
//p‰ivit‰
void kvz_alf_derive_stats_for_filtering(encoder_state_t *const state,
  const lcu_order_element_t *const lcu);

//Onko lumalle ja chromalle sama stride?
//alf_covariace arvot riippuvat funktiosta kvz_alf_calc_covariance ja org- sek‰ reg-arvoista
//Muuten pit‰isi olla OK
void kvz_alf_get_blk_stats(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  channel_type channel,
  alf_covariance *alfCovariace,
  alf_classifier **classifier,
  kvz_pixel *org,
  int32_t org_stride,
  kvz_pixel *rec,
  int32_t rec_stride,
  const int x_pos,
  const int y_pos,
  const int width,
  const int height);

//OK
void kvz_alf_calc_covariance(int *ELocal,
  const kvz_pixel *rec,
  const int stride,
  const int *filterPattern,
  const int halfFilterLength,
  const int transposeIdx);

double kvz_alf_get_filter_coeff_and_cost(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  channel_type channel,
  double dist_unfilter,
  int *ui_coeff_bits,
  int i_shape_idx,
  bool b_re_collect_stat,
  bool only_filter_cost);

int kvz_alf_derive_filter_coefficients_prediction_mode(channel_type channel,
  int **filter_set,
  int** filter_coeff_diff,
  const int num_filters,
  int *pred_mode);

//Ongelmia alf_classifierin kanssa aka ongelma funktiossa kvz_alf_derive_classification_blk
void kvz_alf_merge_classes(channel_type channel,
  alf_covariance* cov,
  alf_covariance* cov_merged,
  const int num_classes,
  short filter_indices[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES]);

double kvz_alf_merge_filters_and_cost(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  channel_type channel,
  int *ui_coeff_bits,
  alf_covariance *cov_frame,
  alf_covariance *cov_merged);

//tarkista
double kvz_alf_derive_filter_coeffs(alf_aps *aps,
  channel_type channel,
  alf_covariance *cov,
  alf_covariance *covMerged,
  short* filter_indices,
  int numFilters,
  double errorTabForce0Coeff[MAX_NUM_ALF_CLASSES][2]);

double kvz_alf_derive_coeff_quant(channel_type channel,
  int *filterCoeffQuant,
  double **ee,
  double *y);

double kvz_alf_derive_ctb_alf_enable_flags(encoder_state_t * const state,
  const lcu_order_element_t *const lcu,
  channel_type channel,
  const int i_shape_idx,
  double *distUnfilter,
  const int numClasses);

void kvz_alf_get_avai_aps_ids_luma(encoder_state_t *const state,
  int *newApsId,
  int aps_ids[6],
  int *size_of_aps_ids);

//----------------------------------------------------------------------

//-------------------------cabac writer functions------------------------
void code_alf_ctu_enable_flag(encoder_state_t * const state,
  const lcu_order_element_t *const lcu,
  uint32_t ctu_rs_addr,
  int8_t component_id,
  alf_aps *aps);

void code_alf_ctu_filter_index(encoder_state_t * const state,
  const lcu_order_element_t *const lcu,
  uint32_t ctu_rs_addr,
  bool alf_enable_luma);
//---------------------------------------------------------------------

//-------------------------CTU functions--------------------------------
//Ei kutsuta koskaan?
//En p‰‰se VTM:ll‰ kyseiseen funktioon
void kvz_alf_process(encoder_state_t const *state,
  const lcu_order_element_t *lcu);

//Pit‰isi luoda APS ( adaptive parameter set )
//kvz_alf_process
//kvz_alf_encoder_ctb
//kvz_alf_reconstructor
void kvz_alf_reconstruct_coeff_aps(encoder_state_t *const state,
  bool luma,
  bool chroma,
  bool is_rdo);

//OK, kun joitakin hardcoded arvoja k‰ytetty
//g_num_luma_filters??
//g_luma_coeff??
//g_filter_coeff_delta_idx??
//g_fixed_filter_set_index??
//g_fixed_filter_idx??
//Ongelmat kantautuu funktiosta kvz_alf_derive_classification_blk
void kvz_alf_reconstruct_coeff(encoder_state_t *const state,
  alf_aps *aps,
  channel_type channel,
  const bool is_rdo,
  const bool is_redo);

//OK
void kvz_alf_create(encoder_state_t const *state,
  const lcu_order_element_t *lcu);

//OK
void kvz_alf_destroy(encoder_state_t const *state);

//OK
//alf_input_bit_depth = 10 VTM:ss‰
void kvz_alf_derive_classification(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  const int width,
  const int height,
  int x_pos,
  int y_pos);

//OK
//Turha jos PCM on pois p‰‰lt‰.
//cu->ipcm?
void kvz_alf_reset_pcm_blk_class_info(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  const int width,
  const int height,
  int x_pos,
  int y_pos);

//OK
void kvz_alf_derive_classification_blk(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  alf_classifier** classifier,
  const int shift,
  const int n_height,
  const int n_width,
  const int x,
  const int y);

//OK
//resetPCMBlkClassInfo-funktiossa VTM:ss‰ return, pois jotta
//p‰‰see t‰h‰n funktioon
void kvz_alf_filter_block(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  short* filter_set,
  clp_rng clp_rng,
  alf_component_id component_id,
  const int width,
  const int height,
  int x_pos,
  int y_pos);
//---------------------------------------------------------------------

//kvz_alf_filter_block funktion dst?
//chroman stride puolet luman stridest‰?
#endif //ALF_H_