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
* \ingroup Reconstruction
* \file
* MIP weight matrix data.
*/

/** \file     MipData.h
\brief    weight and bias data for matrix-based intra prediction (MIP)
*/
#pragma once

#define MIP_SHIFT_MATRIX 6
#define MIP_OFFSET_MATRIX 32
// MIP weight tables for AVX2.

// This is the same table as used in generic version, but 16-bit.
static ALIGNED(32) const uint16_t uvg_mip_sid0_weights[16][16][4] =
{
  {
    {   32,   30,   90,   28},
    {   32,   32,   72,   28},
    {   34,   77,   53,   30},
    {   51,  124,   36,   37},
    {   31,   31,   95,   37},
    {   33,   31,   70,   50},
    {   52,   80,   25,   60},
    {   78,  107,    1,   65},
    {   31,   29,   37,   95},
    {   38,   34,   19,  101},
    {   73,   85,    0,   81},
    {   92,   99,    0,   65},
    {   34,   29,   14,  111},
    {   48,   48,    7,  100},
    {   80,   91,    0,   74},
    {   89,   97,    0,   64}
  },
  {
    {   31,   23,   34,   29},
    {   31,   43,   34,   31},
    {   30,   95,   34,   32},
    {   29,  100,   35,   33},
    {   31,   23,   34,   29},
    {   31,   43,   34,   31},
    {   30,   95,   34,   32},
    {   29,   99,   35,   33},
    {   31,   24,   35,   29},
    {   31,   44,   34,   31},
    {   30,   95,   35,   32},
    {   29,   99,   35,   33},
    {   31,   24,   35,   30},
    {   31,   44,   35,   31},
    {   30,   95,   35,   32},
    {   29,   99,   35,   33}
  },
  {
    {   32,   32,   36,   58},
    {   32,   29,   26,   66},
    {   36,   37,   23,   61},
    {   79,   84,    3,   37},
    {   32,   32,   30,   69},
    {   33,   29,   24,   71},
    {   44,   16,   21,   70},
    {   96,   18,    0,   57},
    {   32,   31,   24,   74},
    {   33,   30,   23,   71},
    {   36,   24,   24,   71},
    {   59,    9,   16,   68},
    {   32,   32,   23,   75},
    {   33,   30,   24,   70},
    {   32,   30,   25,   71},
    {   36,   26,   25,   70}
  },
  {
    {   32,   33,   34,   32},
    {   32,   30,   22,   38},
    {   29,   46,   25,   38},
    {   53,  123,   28,   22},
    {   32,   33,   30,   37},
    {   32,   30,   21,   38},
    {   32,   40,   24,   38},
    {   64,  116,   26,   17},
    {   32,   32,   23,   49},
    {   32,   30,   21,   39},
    {   34,   39,   24,   37},
    {   72,  109,   23,   16},
    {   33,   31,   17,   60},
    {   32,   31,   21,   39},
    {   35,   41,   24,   37},
    {   72,  106,   22,   18}
  },
  {
    {   34,   25,   89,   20},
    {   38,   32,   47,   24},
    {   40,   86,   29,   27},
    {   38,   98,   32,   29},
    {   34,   31,   94,   40},
    {   44,   25,   83,   27},
    {   54,   72,   43,   16},
    {   47,   94,   33,   22},
    {   33,   31,   36,   94},
    {   43,   23,   51,   76},
    {   62,   55,   64,   25},
    {   57,   89,   38,   15},
    {   32,   32,   28,  101},
    {   38,   26,   33,   94},
    {   55,   38,   68,   47},
    {   59,   80,   52,   16}
  },
  {
    {   28,   30,   68,   29},
    {   23,   48,   23,   48},
    {   39,   98,   16,   42},
    {   84,   86,   20,   17},
    {   25,   31,   52,   74},
    {   38,   68,    5,   70},
    {   95,   78,    7,   21},
    {  127,   54,   12,    0},
    {   30,   47,   14,  107},
    {   79,   76,    0,   53},
    {  127,   59,    7,    1},
    {  127,   51,    9,    0},
    {   50,   71,    1,   96},
    {  109,   69,    7,   25},
    {  127,   56,    9,    0},
    {  123,   53,   13,    0}
  },
  {
    {   40,   20,   72,   18},
    {   48,   29,   44,   18},
    {   53,   81,   35,   18},
    {   48,   96,   33,   22},
    {   45,   23,   79,   49},
    {   61,   21,   56,   49},
    {   72,   52,   32,   48},
    {   65,   69,   20,   50},
    {   41,   27,   29,   96},
    {   49,   22,   28,   94},
    {   52,   22,   28,   93},
    {   49,   27,   27,   92},
    {   37,   29,   26,   98},
    {   39,   28,   28,   97},
    {   38,   28,   30,   97},
    {   38,   29,   30,   95}
  },
  {
    {   33,   27,   43,   27},
    {   32,   29,   31,   31},
    {   31,   73,   33,   31},
    {   35,  104,   34,   28},
    {   32,   30,   63,   22},
    {   33,   26,   33,   29},
    {   33,   57,   33,   30},
    {   37,  100,   35,   27},
    {   32,   31,   85,   25},
    {   34,   25,   39,   25},
    {   35,   39,   32,   28},
    {   40,   91,   35,   25},
    {   32,   30,   77,   50},
    {   34,   26,   54,   22},
    {   37,   31,   34,   27},
    {   45,   75,   34,   23}
  },
  {
    {   34,   25,   77,   19},
    {   36,   34,   56,   24},
    {   41,   83,   39,   30},
    {   47,   96,   28,   35},
    {   34,   31,   70,   65},
    {   38,   29,   53,   77},
    {   43,   36,   37,   83},
    {   48,   39,   28,   83},
    {   33,   31,   31,   98},
    {   33,   31,   30,   99},
    {   34,   30,   31,   98},
    {   36,   29,   31,   96},
    {   32,   32,   30,   97},
    {   32,   32,   31,   96},
    {   31,   33,   33,   96},
    {   32,   33,   34,   94}
  },
  {
    {   30,   30,   93,   19},
    {   31,   59,   67,   34},
    {   31,   79,   36,   59},
    {   30,   67,   17,   79},
    {   30,   38,   68,   69},
    {   29,   40,   43,   91},
    {   26,   35,   32,  101},
    {   23,   32,   30,  101},
    {   26,   34,   30,  101},
    {   23,   33,   30,  102},
    {   20,   32,   31,  102},
    {   18,   33,   32,  102},
    {   23,   33,   31,  100},
    {   20,   34,   32,  100},
    {   18,   35,   33,  100},
    {   18,   35,   33,  100}
  },
  {
    {   31,   54,   90,   26},
    {   32,   60,   53,   61},
    {   34,   49,   37,   84},
    {   34,   39,   35,   89},
    {   35,   38,   41,   88},
    {   35,   35,   32,   96},
    {   35,   31,   33,   96},
    {   35,   32,   35,   94},
    {   34,   34,   30,   97},
    {   35,   32,   33,   95},
    {   35,   32,   34,   94},
    {   35,   34,   34,   93},
    {   34,   34,   34,   93},
    {   35,   34,   34,   93},
    {   35,   34,   34,   92},
    {   36,   34,   35,   91}
  },
  {
    {   32,   29,   54,   24},
    {   31,   32,   34,   29},
    {   31,   43,   34,   29},
    {   32,   67,   36,   28},
    {   31,   34,   69,   37},
    {   31,   35,   46,   33},
    {   30,   35,   39,   33},
    {   30,   42,   39,   36},
    {   31,   35,   39,   88},
    {   30,   38,   41,   84},
    {   30,   39,   40,   81},
    {   39,   46,   38,   78},
    {   31,   36,   34,   96},
    {   34,   38,   37,   93},
    {   55,   42,   38,   82},
    {   89,   53,   38,   65}
  },
  {
    {   32,   33,   43,   29},
    {   32,   30,   29,   33},
    {   31,   47,   31,   33},
    {   33,  100,   31,   31},
    {   32,   33,   74,   25},
    {   32,   32,   34,   31},
    {   32,   33,   30,   33},
    {   32,   68,   30,   32},
    {   32,   31,   91,   40},
    {   32,   32,   58,   26},
    {   31,   31,   30,   32},
    {   31,   42,   30,   33},
    {   32,   31,   49,   85},
    {   32,   31,   83,   35},
    {   31,   33,   48,   29},
    {   31,   36,   32,   33}
  },
  {
    {   31,   29,   81,   35},
    {   32,   28,   34,   50},
    {   31,   75,   16,   43},
    {   34,  103,   29,   32},
    {   32,   32,   53,   78},
    {   31,   28,   36,   88},
    {   30,   52,   18,   73},
    {   52,   88,   17,   35},
    {   32,   32,   35,   94},
    {   30,   31,   35,   95},
    {   36,   29,   31,   92},
    {  100,   43,   16,   40},
    {   32,   32,   35,   93},
    {   30,   32,   38,   93},
    {   55,   18,   37,   83},
    {  127,    0,   30,   40}
  },
  {
    {   31,   22,   47,   30},
    {   31,   48,   25,   34},
    {   30,   95,   31,   32},
    {   32,  103,   33,   32},
    {   30,   24,   57,   31},
    {   30,   47,   26,   34},
    {   31,   95,   31,   32},
    {   43,   97,   35,   25},
    {   29,   26,   44,   63},
    {   37,   38,   24,   47},
    {   74,   63,   28,   20},
    {  110,   58,   34,    3},
    {   46,   22,    5,  108},
    {   93,    5,    9,   77},
    {  127,    0,   17,   52},
    {  127,    0,   15,   50}
  },
  {
    {   32,   27,   68,   24},
    {   35,   23,   35,   28},
    {   35,   64,   29,   29},
    {   37,  104,   33,   28},
    {   32,   32,   91,   40},
    {   36,   23,   67,   36},
    {   49,   23,   39,   28},
    {   60,   67,   30,   20},
    {   32,   32,   36,   95},
    {   35,   29,   38,   93},
    {   50,   16,   30,   84},
    {   72,   16,   15,   65},
    {   32,   32,   27,  100},
    {   33,   32,   29,  100},
    {   37,   29,   30,   98},
    {   48,   21,   29,   90}
  }
};

// Weight vectors for MIP size_id 1.
static ALIGNED(32) const uint16_t uvg_mip_sid1_weights[] = {
 30,  63,  30,  60,  29,  45,  30,  39,  46,  37,  66,  38,  74,  42,  62,  58,  // mode 0, offset 0
 25,  33,  32,  31,  32,  32,  32,  33,  33,  34,  32,  33,  32,  33,  32,  33,
 30,  66,  29,  54,  28,  48,  28,  41,  55,  39,  69,  40,  71,  43,  72,  46,
 32,  30,  33,  31,  32,  33,  32,  34,  30,  36,  31,  33,  32,  33,  32,  33,
 30,  66,  29,  55,  27,  46,  27,  42,  56,  40,  69,  39,  72,  43,  69,  48,
 32,  33,  33,  33,  33,  33,  32,  34,  28,  33,  30,  32,  32,  33,  32,  33,
 30,  63,  29,  56,  27,  47,  27,  42,  55,  40,  66,  40,  69,  44,  65,  50,
 32,  33,  33,  33,  33,  33,  32,  34,  35,  30,  33,  30,  33,  32,  32,  33,
 32,  33,  33,  56,  33,  77,  33,  37,  30,  31,  28,  30,  52,  26,  80,  41,  // mode 1, offset 128
 74,  30,  41,  29,  29,  34,  31,  34,  31,  32,  32,  32,  30,  32,  30,  32,
 32,  32,  33,  31,  33,  47,  33,  61,  33,  31,  31,  30,  28,  29,  44,  28,
 59,  76,  78,  40,  53,  27,  34,  32,  28,  31,  28,  32,  31,  31,  31,  31,
 32,  31,  32,  31,  33,  27,  33,  33,  34,  30,  34,  29,  34,  29,  34,  30,
 26,  64,  45,  86,  73,  55,  62,  33,  76,  27,  36,  29,  25,  32,  30,  31,
 32,  31,  32,  31,  32,  30,  33,  28,  34,  30,  35,  29,  36,  29,  37,  30,
 30,  29,  27,  53,  40,  80,  58,  60,  58,  74,  77,  35,  44,  31,  31,  33,
 32,  51,  32,  95,  32,  27,  32,  34,  27,  32,  42,  29,  99,  34,  21, 104,  // mode 2, offset 256
 27,  50,  29,  42,  31,  41,  31,  42,  29,  32,  30,  32,  29,  32,  30,  32,
 32,  45,  32,  77,  32,  38,  32,  30,  30,  32,  38,  30,  78,  33,  30,  87,
  9,  88,   9,  76,  14,  67,  20,  59,  40,  30,  38,  30,  37,  30,  38,  31,
 33,  37,  34,  44,  36,  39,  37,  31,  32,  32,  34,  31,  45,  31,  31,  54,
 27,  18,  25,  17,  24,  15,  25,  14, 106,  34, 108,  31, 108,  30, 101,  32,
 36,  33,  39,  32,  44,  33,  47,  30,  32,  30,  32,  29,  31,  27,  31,  32,
 29,  37,  27,  37,  25,  37,  25,  34,  13, 110,  15, 108,  16, 106,  19, 102,
 32,  48,  32,  33,  32,  29,  33,  33,  35,  35,  59,  40,  47,  65,  31,  81,  // mode 3, offset 384
 47,  68,  27,  71,  24,  62,  26,  50,  31,  31,  33,  30,  37,  30,  42,  32,
 32,  30,  32,  20,  33,  30,  36,  34,  40,  38,  46,  50,  29,  66,  27,  69,
 30,  70,  26,  55,  25,  41,  26,  31,  55,  31,  64,  31,  72,  33,  67,  39,
 33,  28,  36,  27,  43,  30,  51,  27,  36,  40,  33,  50,  26,  57,  28,  55,
 30,  26,  31,  20,  28,  17,  22,  23,  85,  47,  79,  53,  67,  62,  49,  70,
 38,  29,  51,  31,  69,  23,  77,  13,  32,  39,  28,  43,  30,  40,  35,  38,
 28,  30,  24,  31,  15,  38,   8,  43,  22, 104,  17, 102,  10,  95,   8,  90,
 32,  38,  32,  40,  32,  37,  33,  34,  32,  33,  37,  32,  46,  35,  30,  62,  // mode 4, offset 512
101,  40, 100,  36,  94,  33,  81,  35,  29,  32,  30,  32,  30,  31,  30,  31,
 32,  32,  32,  31,  33,  33,  33,  32,  33,  32,  33,  33,  33,  33,  34,  36,
 22, 102,  26, 104,  31, 103,  37,  94,  39,  29,  34,  28,  32,  28,  33,  28,
 32,  33,  32,  34,  33,  33,  33,  33,  32,  32,  33,  33,  34,  33,  33,  36,
 34,  24,  33,  30,  31,  37,  30,  46,  99,  36,  98,  32,  95,  29,  85,  31,
 32,  33,  32,  34,  32,  33,  33,  33,  32,  33,  33,  33,  34,  34,  32,  37,
 30,  34,  31,  32,  31,  29,  32,  30,  23, 104,  30,  98,  39,  91,  47,  82,
 32,  52,  33,  19,  33,  30,  34,  35,  48,  31,  62,  50,  20,  74,  23,  56,  // mode 5, offset 640
 38,  76,  25,  50,  29,  29,  31,  25,  26,  32,  51,  31,  54,  51,  41,  76,
 33,  25,  35,  28,  37,  35,  38,  32,  38,  39,  25,  47,  22,  38,  33,  29,
 28,  39,  31,  23,  31,  27,  30,  31,  83,  35,  57,  74,  30, 101,  27, 103,
 34,  32,  38,  33,  40,  32,  40,  32,  27,  37,  28,  32,  33,  27,  34,  27,
 32,  25,  30,  31,  29,  33,  28,  33,  41,  92,  18, 111,  18, 111,  23, 105,
 35,  32,  38,  31,  40,  32,  40,  32,  30,  33,  33,  30,  33,  29,  33,  30,
 31,  33,  29,  33,  29,  34,  29,  34,  20, 107,  21, 106,  22, 105,  24, 101,
 32,  28,  33,  30,  33,  60,  33,  63,  31,  33,  28,  33,  26,  33,  44,  36,  // mode 6, offset 768
 92,  33,  71,  26,  47,  28,  37,  31,  30,  31,  32,  30,  33,  30,  33,  30,
 33,  30,  33,  28,  33,  30,  33,  38,  31,  33,  29,  34,  26,  33,  29,  32,
 43,  90,  71,  71,  86,  45,  74,  32,  33,  29,  26,  30,  28,  30,  33,  29,
 33,  32,  34,  31,  34,  31,  33,  32,  30,  32,  29,  33,  29,  33,  28,  34,
 29,  41,  26,  71,  37,  88,  55,  75,  95,  27,  73,  22,  46,  25,  36,  28,
 34,  31,  35,  32,  34,  33,  34,  34,  30,  32,  28,  33,  28,  33,  28,  34,
 33,  27,  33,  23,  30,  35,  33,  53,  43,  89,  77,  59,  91,  37,  74,  31,
 33,  49,  33,  71,  32,  23,  31,  33,  26,  32,  72,  24,  70,  68,  21, 106,  // mode 7, offset 896
 26,  52,  30,  32,  32,  32,  33,  32,  28,  31,  34,  31,  32,  32,  32,  33,
 34,  47,  34,  44,  32,  27,  30,  33,  32,  29,  89,  28,  46,  89,  20, 107,
  5,  86,  28,  37,  33,  31,  33,  33,  44,  26,  33,  30,  31,  32,  32,  33,
 35,  39,  34,  27,  31,  31,  29,  32,  42,  27,  87,  43,  32, 100,  22, 106,
 26,  24,  30,  34,  32,  33,  33,  33,  92,  35,  38,  31,  30,  32,  32,  33,
 35,  29,  34,  24,  31,  33,  29,  33,  47,  32,  69,  60,  31,  99,  25, 103,
 32,  32,  34,  33,  32,  33,  33,  33,  17, 100,  28,  44,  32,  31,  32,  35,
};


// Weight vectors for MIP size_id 2.
static ALIGNED(32) const uint16_t uvg_mip_sid2_weights[] = {
  0,  42,   0,  71,   0,  77,   0,  64,  37,  33,  39,  34,  46,  35,  60,  35,  // mode 0, offset 0
 27,  44,  24,  36,  33,  30,  33,  31,  33,  35,  35,  36,  34,  36,  32,  36,
  0,  49,   0,  42,   0,  40,   0,  38,  71,  38,  66,  50,  52,  67,  43,  75,
 32,  32,  33,  31,  33,  31,  33,  32,  31,  36,  32,  36,  32,  35,  32,  35,
  0,  56,   0,  70,   0,  65,   0,  59,  40,  33,  49,  34,  57,  36,  60,  39,
 26,  43,  30,  28,  34,  28,  33,  30,  38,  36,  38,  38,  33,  39,  31,  38,
  0,  55,   0,  51,   0,  46,   0,  42,  60,  43,  61,  47,  62,  51,  60,  55,
 33,  30,  33,  30,  34,  30,  33,  31,  31,  38,  32,  37,  32,  37,  32,  37,
  0,  60,   0,  68,   0,  62,   0,  58,  42,  34,  52,  35,  58,  37,  59,  41,
 30,  37,  35,  22,  34,  28,  33,  30,  43,  38,  37,  40,  31,  40,  30,  39,
  0,  56,   0,  53,   0,  49,   0,  45,  59,  44,  60,  45,  65,  45,  64,  47,
 34,  30,  33,  30,  33,  30,  33,  31,  31,  38,  31,  38,  31,  38,  32,  38,
  0,  59,   0,  66,   0,  61,   0,  59,  44,  35,  53,  36,  58,  38,  57,  41,
 31,  34,  35,  25,  34,  29,  33,  30,  43,  41,  31,  43,  30,  40,  31,  39,
  0,  57,   0,  54,   0,  51,   0,  48,  58,  43,  61,  43,  64,  43,  64,  45,
 33,  30,  33,  31,  33,  31,  33,  32,  31,  39,  31,  39,  31,  39,  31,  39,
  0,  57,   0,  65,   0,  63,   0,  61,  45,  35,  54,  37,  56,  38,  56,  41,
 30,  35,  33,  33,  34,  30,  34,  30,  40,  44,  24,  44,  29,  39,  32,  39,
  0,  58,   0,  54,   0,  51,   0,  48,  58,  42,  62,  41,  65,  42,  63,  43,
 33,  31,  33,  31,  33,  31,  33,  32,  31,  39,  31,  39,  31,  39,  31,  39,
  0,  55,   0,  65,   0,  65,   0,  63,  46,  35,  53,  37,  54,  38,  55,  39,
 30,  36,  32,  36,  33,  31,  33,  30,  38,  47,  26,  40,  30,  38,  32,  38,
  0,  59,   0,  54,   0,  49,   0,  48,  58,  40,  64,  40,  66,  40,  64,  42,
 33,  31,  33,  31,  32,  32,  32,  32,  31,  39,  30,  40,  30,  41,  30,  41,
  0,  54,   0,  64,   0,  65,   0,  63,  46,  35,  52,  36,  53,  37,  55,  38,
 30,  34,  32,  34,  33,  32,  33,  31,  39,  49,  34,  35,  32,  37,  31,  39,
  0,  59,   0,  54,   0,  49,   0,  47,  60,  38,  64,  38,  66,  39,  64,  42,
 33,  31,  33,  32,  33,  32,  32,  33,  31,  40,  30,  40,  29,  41,  29,  42,
  0,  51,   0,  61,   0,  63,   0,  62,  46,  35,  51,  36,  53,  37,  55,  37,
 31,  33,  32,  33,  32,  32,  33,  32,  37,  54,  38,  36,  34,  37,  32,  39,
  0,  58,   0,  53,   0,  49,   0,  46,  59,  37,  63,  38,  64,  40,  62,  42,
 33,  32,  33,  32,  33,  33,  33,  33,  31,  40,  31,  40,  30,  41,  30,  42,
  0,  39,   0,  60,   0,  73,   0,  60,  34,  33,  38,  32,  49,  31,  73,  30,  // mode 1, offset 512
 58,  44,  40,  51,  39,  48,  39,  46,  31,  32,  30,  31,  32,  31,  33,  32,
  0,  43,   0,  35,   0,  33,   0,  31,  87,  35,  78,  54,  47,  86,  17, 114,
 38,  45,  36,  45,  35,  44,  34,  44,  33,  32,  33,  32,  33,  32,  34,  33,
  0,  43,   0,  53,   0,  52,   0,  46,  37,  32,  50,  30,  66,  30,  78,  35,
 53,  70,  42,  72,  39,  70,  37,  68,  30,  31,  31,  30,  32,  30,  34,  30,
  0,  43,   0,  40,   0,  33,   0,  26,  75,  48,  62,  68,  37,  97,  14, 122,
 37,  66,  35,  65,  33,  62,  32,  59,  34,  30,  35,  30,  37,  31,  38,  33,
  0,  40,   0,  45,   0,  41,   0,  37,  39,  33,  54,  32,  70,  35,  73,  44,
 34,  87,  34,  84,  33,  83,  32,  82,  37,  30,  41,  29,  40,  29,  40,  30,
  0,  37,   0,  35,   0,  28,   0,  19,  65,  60,  48,  82,  27, 108,  11, 127,
 31,  81,  30,  79,  28,  76,  27,  70,  41,  29,  43,  29,  45,  30,  46,  32,
  0,  38,   0,  39,   0,  33,   0,  30,  40,  34,  54,  35,  65,  41,  65,  53,
 27,  73,  30,  73,  29,  75,  27,  76,  62,  28,  62,  28,  59,  28,  58,  29,
  0,  29,   0,  27,   0,  19,   0,   9,  53,  72,  35,  95,  19, 117,  16, 127,
 26,  77,  24,  77,  23,  74,  23,  68,  58,  29,  60,  28,  61,  30,  60,  34,
  0,  35,   0,  33,   0,  28,   0,  24,  40,  35,  51,  39,  57,  49,  52,  65,
 29,  44,  29,  49,  28,  53,  26,  56,  89,  30,  86,  30,  83,  30,  82,  30,
  0,  22,   0,  18,   0,  10,   0,   0,  39,  86,  22, 108,  13, 125,  19, 127,
 24,  58,  23,  59,  22,  58,  22,  56,  82,  30,  82,  31,  80,  33,  74,  40,
  0,  33,   0,  29,   0,  24,   0,  19,  40,  36,  46,  44,  45,  58,  37,  78,
 31,  28,  29,  31,  28,  34,  26,  37,  90,  45,  92,  43,  91,  43,  91,  43,
  0,  15,   0,  11,   0,   2,   0,   0,  22,  99,  11, 118,  11, 127,  17, 127,
 25,  38,  24,  39,  23,  41,  23,  43,  91,  42,  90,  44,  85,  48,  75,  55,
  0,  31,   0,  27,   0,  22,   0,  15,  37,  39,  37,  52,  30,  70,  19,  91,
 30,  28,  28,  30,  27,  32,  26,  33,  54,  82,  58,  79,  58,  79,  58,  79,
  0,  10,   0,   5,   0,   0,   0,   0,   8, 111,   2, 125,   9, 127,  13, 127,
 25,  34,  25,  35,  25,  36,  25,  39,  58,  79,  57,  80,  53,  84,  47,  88,
  0,  28,   0,  24,   0,  19,   0,  13,  29,  46,  24,  62,  14,  81,   4, 101,
 28,  39,  27,  41,  25,  43,  24,  44,   2, 123,   1, 125,   0, 126,   0, 127,
  0,   6,   0,   0,   0,   0,   0,   0,   0, 116,   0, 126,   4, 127,   9, 127,
 23,  45,  23,  45,  25,  44,  25,  44,   0, 127,   1, 127,   2, 127,   3, 127,
  0,  30,   0,  63,   0,  98,   0,  75,  32,  32,  26,  34,  26,  34,  61,  30,  // mode 2, offset 1024
 42,  34,  16,  38,  25,  34,  31,  32,  32,  32,  32,  32,  33,  32,  33,  32,
  0,  36,   0,  26,   0,  30,   0,  32,  94,  32,  76,  58,  39,  91,  23, 105,
 30,  33,  30,  33,  31,  32,  32,  32,  32,  32,  32,  32,  33,  31,  32,  32,
  0,  34,   0,  66,   0,  97,   0,  71,  30,  33,  24,  34,  28,  34,  65,  30,
 31,  52,  11,  41,  24,  34,  30,  32,  29,  32,  33,  32,  33,  32,  33,  32,
  0,  34,   0,  26,   0,  30,   0,  32,  92,  35,  70,  64,  37,  94,  23, 105,
 30,  33,  29,  34,  30,  33,  31,  33,  32,  32,  32,  32,  32,  31,  33,  31,
  0,  37,   0,  71,   0,  98,   0,  66,  29,  33,  22,  35,  29,  34,  70,  30,
  8,  79,   5,  50,  23,  34,  31,  31,  27,  32,  32,  32,  34,  32,  33,  32,
  0,  31,   0,  26,   0,  30,   0,  31,  92,  38,  66,  68,  34,  97,  22, 106,
 30,  33,  29,  34,  30,  34,  30,  34,  32,  32,  32,  31,  33,  31,  33,  31,
  0,  40,   0,  76,   0,  97,   0,  61,  28,  34,  21,  35,  32,  34,  75,  29,
  0,  76,   0,  55,  21,  37,  30,  32,  46,  28,  35,  32,  33,  33,  32,  32,
  0,  29,   0,  26,   0,  29,   0,  31,  92,  40,  62,  73,  32,  99,  22, 107,
 29,  33,  29,  34,  30,  34,  30,  34,  32,  32,  32,  31,  33,  30,  33,  31,
  0,  42,   0,  80,   0,  94,   0,  55,  27,  34,  20,  35,  36,  32,  80,  29,
  1,  48,   0,  48,  17,  40,  27,  35,  79,  25,  47,  31,  33,  33,  31,  32,
  0,  27,   0,  26,   0,  29,   0,  31,  90,  43,  58,  76,  30, 101,  21, 108,
 28,  34,  29,  33,  29,  34,  29,  35,  32,  31,  33,  30,  34,  30,  34,  30,
  0,  44,   0,  81,   0,  90,   0,  51,  26,  34,  21,  35,  41,  31,  82,  29,
  6,  30,   0,  41,  14,  41,  24,  37,  80,  40,  52,  35,  35,  33,  32,  32,
  0,  27,   0,  26,   0,  29,   0,  31,  87,  47,  54,  79,  29, 102,  21, 108,
 27,  35,  29,  34,  28,  34,  28,  35,  32,  31,  33,  30,  33,  30,  33,  31,
  0,  47,   0,  80,   0,  84,   0,  49,  26,  34,  24,  34,  45,  31,  81,  31,
  7,  34,   0,  41,  12,  40,  22,  37,  44,  75,  41,  50,  36,  36,  33,  32,
  0,  28,   0,  28,   0,  29,   0,  31,  81,  51,  51,  81,  30, 101,  22, 107,
 26,  35,  28,  34,  28,  35,  28,  35,  33,  31,  33,  30,  33,  31,  33,  32,
  0,  48,   0,  75,   0,  77,   0,  49,  27,  34,  27,  34,  47,  33,  75,  36,
 10,  40,   3,  42,  12,  40,  21,  37,  16,  97,  26,  66,  32,  43,  33,  35,
  0,  32,   0,  30,   0,  30,   0,  31,  72,  55,  49,  81,  32,  98,  24, 104,
 25,  36,  27,  35,  28,  35,  28,  35,  33,  32,  33,  31,  32,  32,  32,  33,
  0,  36,   0,  74,   0,  92,   0,  53,  29,  33,  20,  35,  35,  32,  80,  26,  // mode 3, offset 1536
 43,  47,  19,  47,  29,  31,  33,  28,  29,  31,  34,  32,  40,  34,  36,  37,
  0,  24,   0,  25,   0,  32,   0,  34,  91,  41,  57,  74,  28,  99,  20, 105,
 31,  31,  31,  32,  32,  32,  33,  32,  31,  38,  30,  37,  29,  36,  30,  35,
  0,  50,   0,  75,   0,  64,   0,  31,  26,  34,  28,  33,  58,  29,  85,  37,
 33,  74,  23,  46,  30,  26,  31,  27,  30,  31,  47,  33,  46,  40,  33,  44,
  0,  22,   0,  29,   0,  33,   0,  34,  67,  64,  35,  93,  20, 105,  19, 106,
 30,  31,  31,  32,  32,  33,  33,  32,  28,  42,  27,  40,  27,  37,  29,  36,
  0,  51,   0,  61,   0,  40,   0,  22,  29,  33,  42,  31,  70,  34,  72,  54,
 25,  72,  30,  31,  32,  24,  30,  31,  51,  30,  60,  39,  41,  50,  27,  50,
  0,  25,   0,  32,   0,  34,   0,  34,  44,  83,  23, 102,  18, 107,  19, 105,
 30,  33,  32,  33,  32,  33,  33,  32,  25,  44,  26,  40,  28,  37,  30,  35,
  0,  45,   0,  43,   0,  27,   0,  22,  35,  32,  53,  33,  67,  45,  53,  72,
 30,  39,  35,  24,  32,  29,  30,  33,  79,  33,  53,  55,  27,  61,  22,  52,
  0,  28,   0,  32,   0,  34,   0,  34,  31,  95,  20, 105,  18, 107,  20, 105,
 31,  33,  32,  33,  32,  32,  33,  31,  25,  43,  27,  38,  29,  36,  31,  35,
  0,  38,   0,  31,   0,  22,   0,  25,  40,  32,  55,  39,  57,  60,  39,  86,
 35,  23,  34,  29,  31,  35,  31,  35,  72,  54,  32,  73,  18,  64,  22,  49,
  0,  30,   0,  33,   0,  34,   0,  34,  24, 101,  19, 106,  18, 107,  20, 104,
 32,  33,  32,  32,  33,  31,  33,  31,  27,  40,  30,  36,  31,  35,  32,  34,
  0,  33,   0,  26,   0,  23,   0,  27,  42,  35,  51,  50,  46,  74,  32,  93,
 34,  28,  33,  34,  31,  35,  32,  34,  39,  82,  18,  80,  20,  59,  26,  44,
  0,  31,   0,  33,   0,  34,   0,  35,  22, 103,  19, 106,  19, 106,  21, 103,
 32,  32,  33,  31,  33,  31,  34,  31,  30,  37,  31,  35,  32,  34,  32,  34,
  0,  29,   0,  24,   0,  24,   0,  28,  41,  41,  44,  62,  37,  83,  28,  97,
 33,  34,  34,  35,  34,  33,  33,  32,  20,  92,  18,  73,  25,  52,  30,  40,
  0,  32,   0,  34,   0,  35,   0,  35,  23, 103,  20, 105,  20, 104,  22, 102,
 33,  31,  34,  30,  34,  30,  34,  30,  32,  36,  33,  34,  33,  33,  33,  34,
  0,  27,   0,  26,   0,  27,   0,  30,  38,  51,  37,  71,  33,  87,  28,  96,
 34,  34,  35,  34,  35,  32,  34,  31,  20,  86,  24,  64,  30,  47,  32,  39,
  0,  32,   0,  34,   0,  35,   0,  34,  24, 100,  23, 101,  23, 101,  24,  99,
 35,  30,  34,  30,  34,  30,  35,  30,  32,  36,  33,  34,  32,  34,  33,  34,
  0,  39,   0,  72,   0, 100,   0,  75,  30,  31,  21,  32,  23,  32,  63,  24,  // mode 4, offset 2048
 67,  33,  43,  39,  35,  39,  32,  38,  34,  31,  33,  31,  34,  31,  34,  32,
  0,  32,   0,  22,   0,  31,   0,  35,  98,  26,  77,  55,  37,  90,  22, 100,
 29,  37,  29,  36,  31,  35,  33,  33,  35,  32,  35,  31,  35,  32,  36,  33,
  0,  47,   0,  71,   0,  86,   0,  65,  29,  32,  24,  32,  31,  30,  63,  25,
 74,  54,  60,  50,  46,  48,  34,  46,  32,  31,  36,  30,  37,  30,  39,  30,
  0,  33,   0,  26,   0,  33,   0,  37,  85,  32,  64,  60,  33,  87,  23,  93,
 28,  43,  27,  39,  29,  35,  32,  33,  40,  30,  41,  30,  41,  31,  41,  32,
  0,  41,   0,  55,   0,  62,   0,  53,  32,  32,  31,  32,  37,  31,  55,  31,
 45,  84,  50,  70,  45,  61,  36,  55,  32,  32,  40,  30,  45,  29,  48,  29,
  0,  38,   0,  34,   0,  38,   0,  40,  63,  40,  49,  60,  30,  78,  24,  83,
 29,  48,  27,  43,  28,  38,  30,  36,  50,  28,  51,  29,  50,  31,  48,  33,
  0,  35,   0,  39,   0,  41,   0,  41,  33,  33,  35,  33,  39,  34,  43,  37,
 29,  75,  34,  68,  36,  61,  33,  54,  58,  29,  59,  29,  62,  29,  64,  28,
  0,  41,   0,  42,   0,  42,   0,  42,  43,  45,  36,  56,  30,  65,  28,  68,
 30,  48,  27,  44,  27,  41,  28,  37,  65,  29,  63,  30,  60,  33,  56,  36,
  0,  33,   0,  31,   0,  31,   0,  35,  34,  33,  36,  34,  37,  35,  35,  39,
 31,  42,  31,  44,  32,  43,  32,  40,  88,  30,  84,  31,  83,  31,  82,  31,
  0,  40,   0,  44,   0,  44,   0,  43,  32,  44,  30,  48,  30,  52,  30,  55,
 31,  38,  30,  37,  28,  37,  29,  35,  81,  31,  78,  33,  72,  36,  66,  40,
  0,  32,   0,  30,   0,  30,   0,  33,  33,  33,  34,  34,  34,  36,  32,  38,
 34,  25,  33,  25,  34,  25,  34,  25,  85,  48,  88,  44,  90,  41,  90,  40,
  0,  38,   0,  42,   0,  43,   0,  42,  29,  41,  29,  41,  30,  42,  31,  45,
 34,  26,  33,  27,  31,  28,  31,  30,  88,  40,  85,  41,  80,  43,  72,  47,
  0,  32,   0,  31,   0,  32,   0,  34,  33,  33,  32,  34,  32,  35,  31,  36,
 33,  26,  35,  20,  36,  17,  36,  17,  54,  79,  68,  68,  76,  62,  79,  59,
  0,  37,   0,  39,   0,  41,   0,  40,  29,  37,  29,  37,  30,  37,  31,  40,
 36,  18,  35,  20,  34,  22,  32,  26,  78,  58,  77,  58,  74,  58,  68,  59,
  0,  33,   0,  34,   0,  34,   0,  35,  31,  34,  30,  34,  31,  34,  31,  34,
 33,  29,  35,  23,  36,  20,  36,  18,  31,  98,  45,  88,  54,  82,  59,  78,
  0,  36,   0,  38,   0,  39,   0,  39,  31,  34,  30,  34,  31,  35,  31,  37,
 37,  19,  36,  20,  35,  22,  34,  24,  60,  76,  61,  74,  60,  73,  59,  71,
  0,  30,   0,  47,   0,  81,   0,  85,  33,  32,  30,  31,  28,  32,  46,  29,  // mode 5, offset 2560
 55,  32,  29,  36,  28,  34,  32,  32,  32,  32,  32,  32,  32,  32,  33,  32,
  0,  54,   0,  30,   0,  30,   0,  37,  82,  26,  90,  38,  56,  73,  21, 102,
 32,  32,  31,  32,  31,  33,  32,  32,  33,  32,  33,  32,  32,  32,  32,  32,
  0,  33,   0,  38,   0,  63,   0,  82,  32,  31,  32,  31,  30,  31,  37,  30,
 68,  39,  43,  34,  29,  34,  29,  33,  31,  31,  33,  31,  32,  32,  32,  32,
  0,  71,   0,  44,   0,  33,   0,  37,  63,  27,  86,  30,  72,  55,  37,  86,
 31,  32,  30,  33,  30,  32,  31,  32,  33,  32,  33,  32,  32,  31,  33,  31,
  0,  34,   0,  36,   0,  51,   0,  75,  33,  32,  33,  31,  30,  31,  31,  31,
 60,  61,  56,  38,  38,  33,  30,  33,  29,  32,  32,  31,  33,  32,  33,  32,
  0,  80,   0,  60,   0,  41,   0,  38,  47,  29,  73,  27,  78,  41,  53,  68,
 30,  32,  30,  33,  30,  33,  30,  32,  33,  31,  33,  31,  32,  31,  33,  31,
  0,  33,   0,  35,   0,  43,   0,  64,  33,  32,  33,  31,  32,  31,  30,  31,
 43,  77,  55,  54,  46,  39,  35,  34,  35,  30,  29,  32,  31,  32,  33,  32,
  0,  79,   0,  73,   0,  54,   0,  43,  37,  30,  57,  28,  73,  33,  64,  52,
 31,  32,  30,  32,  30,  32,  30,  32,  33,  31,  33,  31,  33,  31,  33,  31,
  0,  33,   0,  34,   0,  38,   0,  54,  33,  32,  33,  31,  33,  31,  31,  31,
 34,  68,  45,  70,  48,  52,  40,  39,  58,  28,  33,  31,  29,  32,  31,  32,
  0,  73,   0,  77,   0,  65,   0,  51,  32,  31,  45,  29,  63,  30,  66,  42,
 34,  34,  31,  32,  31,  31,  30,  32,  33,  31,  32,  32,  33,  31,  33,  31,
  0,  33,   0,  34,   0,  36,   0,  47,  32,  32,  33,  31,  33,  30,  31,  31,
 34,  44,  38,  66,  44,  62,  43,  48,  81,  31,  52,  28,  34,  31,  30,  32,
  0,  64,   0,  75,   0,  71,   0,  59,  31,  31,  38,  30,  53,  30,  61,  37,
 38,  38,  33,  34,  31,  32,  30,  32,  32,  32,  32,  32,  33,  32,  33,  32,
  0,  33,   0,  34,   0,  36,   0,  43,  32,  31,  33,  31,  33,  31,  32,  31,
 35,  31,  37,  49,  41,  60,  43,  54,  71,  54,  70,  33,  48,  30,  35,  31,
  0,  56,   0,  68,   0,  70,   0,  63,  31,  31,  35,  30,  45,  30,  55,  35,
 40,  44,  36,  37,  33,  34,  31,  33,  32,  32,  32,  32,  33,  32,  33,  32,
  0,  33,   0,  34,   0,  36,   0,  41,  32,  31,  32,  31,  33,  31,  33,  31,
 33,  34,  36,  38,  39,  50,  41,  53,  36,  87,  62,  52,  57,  36,  43,  33,
  0,  50,   0,  59,   0,  65,   0,  62,  33,  31,  35,  31,  42,  31,  49,  35,
 41,  48,  37,  41,  35,  36,  33,  34,  36,  32,  34,  32,  33,  32,  34,  33,
};