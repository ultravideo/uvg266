/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2015 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License version 2.1 as
 * published by the Free Software Foundation.
 *
 * Kvazaar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#include "greatest/greatest.h"

#include "test_strategies.h"

#include "src/image.h"

#include <math.h>
#include <stdlib.h>


//////////////////////////////////////////////////////////////////////////
// MACROS
#define NUM_SIZES 4
#define LCU_MAX_LOG_W 5
#define LCU_MIN_LOG_W 2
#define NUM_TRANSFORM 4
#define NUM_TESTS NUM_TRANSFORM*NUM_SIZES

//////////////////////////////////////////////////////////////////////////
// GLOBALS
static int16_t * dct_bufs[NUM_TESTS] = { 0 }; // SIMD aligned pointers.
static int16_t * dct_actual_bufs[NUM_TESTS] = { 0 }; // pointers returned by malloc.

static int16_t dct_result[NUM_TRANSFORM][NUM_SIZES][LCU_WIDTH*LCU_WIDTH] = { { 0 } };
static int16_t idct_result[NUM_TRANSFORM][NUM_SIZES][LCU_WIDTH*LCU_WIDTH] = { { 0 } };

static struct test_env_t {
  int log_width; // for selecting dim from bufs
  mts_dct_func* tested_func;
  const strategy_t * strategy;
  char msg[1024];
} test_env;


//////////////////////////////////////////////////////////////////////////
// SETUP, TEARDOWN AND HELPER FUNCTIONS
static void init_gradient(int x_px, int y_px, int width, int slope, int16_t *buf)
{
  for (int y = 0; y < width; ++y) {
    for (int x = 0; x < width; ++x) {
      int diff_x = x_px - x;
      int diff_y = y_px - y;
      int val = slope * sqrt(diff_x * diff_x + diff_y * diff_y) + 0.5;
      buf[y * width + x] = CLIP(0, 255, val);
    }
  }
}


static void setup_tests()
{
  for (int test = 0; test < NUM_TESTS; ++test) {

    dct_actual_bufs[test] = malloc(LCU_WIDTH*LCU_WIDTH*sizeof(int16_t) + SIMD_ALIGNMENT);
    dct_bufs[test] = ALIGNED_POINTER(dct_actual_bufs[test], SIMD_ALIGNMENT);
  }

  for (int test = 0; test < NUM_TESTS; ++test) {
      const int width = LCU_WIDTH;
      init_gradient(width, width, width, 255 / width, dct_bufs[test]);
  }

   

  // Select buffer width according to function name for dct function.
  int block = 0;
  for (int s = 0; s < strategies.count; ++s)
  {
    strategy_t *strat = &strategies.strategies[s];
    mts_dct_func* mts_generic = 0;
    if (strcmp(strat->type, "mts_dct") == 0  &&
      strcmp(strat->strategy_name, "generic") == 0)
    {
      mts_generic = strat->fptr;
      for (block = 0; block < NUM_SIZES; block++) {
        for (int trafo = 0; trafo < NUM_TRANSFORM; trafo++) {
          cu_info_t tu;
          tu.type = CU_INTRA;
          tu.tr_idx = MTS_DST7_DST7 + trafo;
          mts_generic(KVZ_BIT_DEPTH, COLOR_Y, &tu, 1 << (LCU_MIN_LOG_W + block), dct_bufs[trafo*NUM_SIZES+block], dct_result[trafo][block], KVZ_MTS_BOTH);
        }
      }      
    }
  }

  block = 0;
  for (int s = 0; s < strategies.count; ++s)
  {
    strategy_t *strat = &strategies.strategies[s];
    mts_idct_func* idct_generic = 0;
    if (strcmp(strat->type, "mts_idct") == 0  &&
      strcmp(strat->strategy_name, "generic") == 0)
    {
      
      idct_generic = strat->fptr;
      for (block = 0; block < NUM_SIZES; block++) {
        for (int trafo = 0; trafo < NUM_TRANSFORM; trafo++) {
          cu_info_t tu;
          tu.type = CU_INTRA;
          tu.tr_idx = MTS_DST7_DST7 + trafo;
          idct_generic(KVZ_BIT_DEPTH, COLOR_Y, &tu, 1 << (LCU_MIN_LOG_W + block), dct_bufs[trafo * NUM_SIZES + block], idct_result[trafo][block], KVZ_MTS_BOTH);
        }
      }
      
    }
  }
}

static void tear_down_tests()
{
  for (int test = 0; test < NUM_TESTS; ++test) {
    free(dct_actual_bufs[test]);
  }
}


//////////////////////////////////////////////////////////////////////////
// TESTS
TEST dct(void)
{
  char testname[100];
  for (int blocksize = 0; blocksize < NUM_SIZES; blocksize++) {
    for (int trafo = 0; trafo < NUM_TRANSFORM; trafo++) {      
      sprintf(testname, "Block: %d x %d, trafo: %d", 1 << (LCU_MIN_LOG_W + blocksize), 1 << (LCU_MIN_LOG_W + blocksize), trafo);
      cu_info_t tu;
      tu.type = CU_INTRA;
      tu.tr_idx = MTS_DST7_DST7 + trafo;

      int16_t* buf = dct_bufs[trafo * NUM_SIZES + blocksize];
      ALIGNED(32) int16_t test_result[LCU_WIDTH * LCU_WIDTH] = { 0 };

      test_env.tested_func(KVZ_BIT_DEPTH, COLOR_Y, &tu, 1 << (LCU_MIN_LOG_W + blocksize), buf, test_result, KVZ_MTS_BOTH);

      for (int i = 0; i < LCU_WIDTH * LCU_WIDTH; ++i) {
        ASSERT_EQm(testname, test_result[i], dct_result[trafo][blocksize][i]);
      }
      //fprintf(stderr, "PASS: %s\r\n", testname);
    }
  }
  
  PASS();
}

TEST idct(void)
{
  char testname[100];
  for (int blocksize = 0; blocksize < NUM_SIZES; blocksize++) {
    for (int trafo = 0; trafo < NUM_TRANSFORM; trafo++) {
      sprintf(testname, "Block: %d x %d, trafo: %d", 1 << (LCU_MIN_LOG_W + blocksize), 1 << (LCU_MIN_LOG_W + blocksize), trafo);
      cu_info_t tu;
      tu.type = CU_INTRA;
      tu.tr_idx = MTS_DST7_DST7 + trafo;

      int16_t* buf = dct_bufs[trafo * NUM_SIZES + blocksize];
      ALIGNED(32) int16_t test_result[LCU_WIDTH * LCU_WIDTH] = { 0 };

      test_env.tested_func(KVZ_BIT_DEPTH, COLOR_Y, &tu, 1 << (LCU_MIN_LOG_W + blocksize), buf, test_result, KVZ_MTS_BOTH);

      for (int i = 0; i < LCU_WIDTH * LCU_WIDTH; ++i) {
        ASSERT_EQm(testname, test_result[i], idct_result[trafo][blocksize][i]);
      }
      //fprintf(stderr, "PASS: %s\r\n", testname);
    }
  }

  PASS();
}


//////////////////////////////////////////////////////////////////////////
// TEST FIXTURES
SUITE(mts_tests)
{
  //SET_SETUP(sad_setup);
  //SET_TEARDOWN(sad_teardown);

  setup_tests();

  // Loop through all strategies picking out the intra sad ones and run
  // select strategies though all tests
  for (volatile unsigned i = 0; i < strategies.count; ++i) {
    const strategy_t * strategy = &strategies.strategies[i];

    test_env.tested_func = strategies.strategies[i].fptr;
    test_env.strategy = strategy;

    // Call different tests depending on type of function.
    // This allows for selecting a subset of tests with -t parameter.
   if (strcmp(strategy->type, "mts_dct") == 0)
    {
      //fprintf(stderr, "Test: %s\r\n", strategy->strategy_name);
      RUN_TEST(dct);
    }
    else if (strcmp(strategy->type, "mts_idct") == 0)
    {
      //fprintf(stderr, "Test: %s\r\n", strategy->strategy_name);
      RUN_TEST(idct);
    }
  }

  tear_down_tests();
}
