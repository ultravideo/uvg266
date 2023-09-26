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

static int16_t dct_result[NUM_TRANSFORM][NUM_SIZES][LCU_WIDTH*LCU_WIDTH] = { { { 0 } } };
static int16_t idct_result[NUM_TRANSFORM][NUM_SIZES][LCU_WIDTH*LCU_WIDTH] = { { { 0 } } };

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
          tu.lfnst_idx = 0;
          tu.cr_lfnst_idx = 0;
          tu.intra.isp_mode = 0;
          mts_generic(UVG_BIT_DEPTH, COLOR_Y, &tu, 1 << (LCU_MIN_LOG_W + block), 1 << (LCU_MIN_LOG_W + block), dct_bufs[trafo*NUM_SIZES+block], dct_result[trafo][block], UVG_MTS_BOTH);
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
          tu.lfnst_idx = 0;
          tu.cr_lfnst_idx = 0;
          tu.intra.isp_mode = 0;
          idct_generic(UVG_BIT_DEPTH, COLOR_Y, &tu, 1 << (LCU_MIN_LOG_W + block), 1 << (LCU_MIN_LOG_W + block), dct_bufs[trafo * NUM_SIZES + block], idct_result[trafo][block], UVG_MTS_BOTH);
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
    size_t size = 1 << (LCU_MIN_LOG_W + blocksize);
    for (int trafo = 0; trafo < NUM_TRANSFORM; trafo++) {      
      sprintf(testname, "Block: %d x %d, trafo: %d", 1 << (LCU_MIN_LOG_W + blocksize), 1 << (LCU_MIN_LOG_W + blocksize), trafo);
      cu_info_t tu;
      tu.type = CU_INTRA;
      tu.tr_idx = MTS_DST7_DST7 + trafo;
      tu.lfnst_idx = 0;
      tu.cr_lfnst_idx = 0;
      tu.intra.isp_mode = 0;

      int16_t* buf = dct_bufs[trafo * NUM_SIZES + blocksize];
      ALIGNED(32) int16_t test_result[LCU_WIDTH * LCU_WIDTH] = { 0 };

      test_env.tested_func(UVG_BIT_DEPTH, COLOR_Y, &tu, 1 << (LCU_MIN_LOG_W + blocksize), 1 << (LCU_MIN_LOG_W + blocksize), buf, test_result, UVG_MTS_BOTH);

      for (int y = 0; y < size; ++y) {
        if (y>= 16) break;
        for (int x = 0; x < size; ++x) {
          if (x >= 16) break;
          int i = y * size + x;
          ASSERT_EQm(testname, test_result[i], dct_result[trafo][blocksize][i]);
        }
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
      tu.lfnst_idx = 0;
      tu.cr_lfnst_idx = 0;
      tu.intra.isp_mode = 0;

      int16_t* buf = dct_bufs[trafo * NUM_SIZES + blocksize];
      ALIGNED(32) int16_t test_result[LCU_WIDTH * LCU_WIDTH] = { 0 };

      test_env.tested_func(UVG_BIT_DEPTH, COLOR_Y, &tu, 1 << (LCU_MIN_LOG_W + blocksize), 1 << (LCU_MIN_LOG_W + blocksize), buf, test_result, UVG_MTS_BOTH);

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
