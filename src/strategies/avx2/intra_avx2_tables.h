#ifndef INTRA_AVX2_TABLES_H
#define INTRA_AVX2_TABLES_H

#include "global.h"

// Test table
ALIGNED(32) const int8_t  intra_chroma_linear_interpolation_w4_m40[] = {
  16,  16,  16,  16,  16,  16,  16,  16,  32,   0,  32,   0,  32,   0,  32,   0,
};

#endif INTRA_AVX2_TABLES_H
