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

#ifndef FAST_COEFF_COST_H_
#define FAST_COEFF_COST_H_

#include <stdio.h>
#include "uvg266.h"
// #include "encoderstate.h"

#define MAX_FAST_COEFF_COST_QP 50

typedef struct {
  uint64_t wts_by_qp[MAX_FAST_COEFF_COST_QP];
} fast_coeff_table_t;

// Weights for 4 buckets (coeff 0, coeff 1, coeff 2, coeff >= 3), for QPs from
// 0 to MAX_FAST_COEFF_COST_QP
static const double default_fast_coeff_cost_wts[][4] = {
  // Just extend it by stretching the first actual values..
  {0.164240f, 4.161530f, 3.509033f, 6.928047f},
  {0.164240f, 4.161530f, 3.509033f, 6.928047f},
  {0.164240f, 4.161530f, 3.509033f, 6.928047f},
  {0.164240f, 4.161530f, 3.509033f, 6.928047f},
  {0.164240f, 4.161530f, 3.509033f, 6.928047f},
  {0.164240f, 4.161530f, 3.509033f, 6.928047f},
  {0.164240f, 4.161530f, 3.509033f, 6.928047f},
  {0.164240f, 4.161530f, 3.509033f, 6.928047f},
  {0.164240f, 4.161530f, 3.509033f, 6.928047f},
  {0.164240f, 4.161530f, 3.509033f, 6.928047f},
  // up to here                          
  {0.164240f, 4.161530f, 3.509033f, 6.928047f},
  {0.162844f, 4.055940f, 3.564467f, 6.861493f},
  {0.128729f, 4.311973f, 3.942837f, 6.935403f},
  {0.110956f, 4.433190f, 3.945753f, 6.877697f},
  {0.095026f, 4.483547f, 4.194173f, 6.781540f},
  {0.075046f, 4.633703f, 4.084193f, 6.698600f},
  {0.052426f, 4.967223f, 4.027210f, 6.549197f},
  {0.040219f, 5.141820f, 3.982650f, 6.461557f},
  {0.035090f, 5.192493f, 3.830950f, 6.418477f},
  {0.029845f, 5.211647f, 3.815457f, 6.345440f},
  {0.023522f, 5.322213f, 3.816537f, 6.360677f},
  {0.021305f, 5.225923f, 3.842700f, 6.325787f},
  {0.015878f, 5.183090f, 3.956003f, 6.329680f},
  {0.010430f, 5.099230f, 4.176803f, 6.305400f},
  {0.008433f, 5.030257f, 4.237587f, 6.270133f},
  {0.006500f, 4.969247f, 4.339397f, 6.217827f},
  {0.004929f, 4.923500f, 4.442413f, 6.183523f},
  {0.003715f, 4.915583f, 4.429090f, 6.125320f},
  {0.003089f, 4.883907f, 4.562790f, 6.156447f},
  {0.002466f, 4.881063f, 4.629883f, 6.142643f},
  {0.002169f, 4.882493f, 4.646313f, 6.127663f},
  {0.002546f, 4.793337f, 4.837413f, 6.199270f},
  {0.001314f, 4.808853f, 4.828337f, 6.243437f},
  {0.001154f, 4.862603f, 4.846883f, 6.205523f},
  {0.000984f, 4.866403f, 4.859330f, 6.240893f},
  {0.000813f, 4.856633f, 4.924527f, 6.293413f},
  {0.001112f, 4.789260f, 5.009880f, 6.433540f},
  {0.000552f, 4.760747f, 5.090447f, 6.599380f},
  {0.000391f, 4.961447f, 5.111033f, 6.756370f},
  {0.000332f, 4.980953f, 5.138127f, 6.867420f},
  {0.000201f, 5.181957f, 4.740160f, 6.460997f},
  {0.000240f, 5.185390f, 4.874840f, 6.819093f},
  {0.000130f, 5.270350f, 4.734213f, 6.826240f},
  {0.000104f, 5.371937f, 4.595087f, 6.659253f},
  {0.000083f, 5.362000f, 4.617470f, 6.837770f},
  {0.000069f, 5.285997f, 4.754993f, 7.159043f},
  {0.000049f, 5.488470f, 4.396107f, 6.727357f},
  {0.000058f, 4.958940f, 4.580460f, 6.477740f},
  {0.000028f, 5.521253f, 4.440493f, 7.205017f},
  {0.000000f, 0.000000f, 0.000000f, 0.000000f},
  {0.000019f, 5.811260f, 4.399110f, 7.336310f},
};

typedef struct encoder_state_t encoder_state_t;

int uvg_fast_coeff_table_parse(fast_coeff_table_t *fast_coeff_table, FILE *fast_coeff_table_f);
void uvg_fast_coeff_use_default_table(fast_coeff_table_t *fast_coeff_table);
uint64_t uvg_fast_coeff_get_weights(const encoder_state_t *state);

#endif // FAST_COEFF_COST_H_
