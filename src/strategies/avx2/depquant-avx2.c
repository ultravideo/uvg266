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

/*
* \file
*/

#include "strategies/avx2/depquant-avx2.h"
#include "strategyselector.h"


#if COMPILE_INTEL_AVX2 && defined X86_64
#include "dep_quant.h"

#include <immintrin.h>
#include "cu.h"
#include "encoderstate.h"
#include "intra.h"
#include "rdo.h"
#include "transform.h"
#include "generic/quant-generic.h"
#include "uvg_math.h"
static const int32_t g_goRiceBits[4][RICEMAX] = {
    { 32768,  65536,  98304, 131072, 163840, 196608, 262144, 262144, 327680, 327680, 327680, 327680, 393216, 393216, 393216, 393216, 393216, 393216, 393216, 393216, 458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752},
    { 65536,  65536,  98304,  98304, 131072, 131072, 163840, 163840, 196608, 196608, 229376, 229376, 294912, 294912, 294912, 294912, 360448, 360448, 360448, 360448, 360448, 360448, 360448, 360448, 425984, 425984, 425984, 425984, 425984, 425984, 425984, 425984},
    { 98304,  98304,  98304,  98304, 131072, 131072, 131072, 131072, 163840, 163840, 163840, 163840, 196608, 196608, 196608, 196608, 229376, 229376, 229376, 229376, 262144, 262144, 262144, 262144, 327680, 327680, 327680, 327680, 327680, 327680, 327680, 327680},
    {131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072, 163840, 163840, 163840, 163840, 163840, 163840, 163840, 163840, 196608, 196608, 196608, 196608, 196608, 196608, 196608, 196608, 229376, 229376, 229376, 229376, 229376, 229376, 229376, 229376},
};

static const int g_riceT[4] = { 32,128, 512, 2048 };
static const int g_riceShift[5] = { 0, 2, 4, 6, 8 };

static const uint32_t g_goRiceParsCoeff[32] = { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2,
                                         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3 };

static void check_rd_costs_avx2(const all_depquant_states* const state, const enum ScanPosType spt, const PQData* pqDataA, Decision* decisions, int start)
{
  int64_t temp_rd_cost_a[4] = {0, 0, 0, 0};
  int64_t temp_rd_cost_b[4] = {0, 0, 0, 0};
  int64_t temp_rd_cost_z[4] = {0, 0, 0, 0};

  __m256i pq_a_delta_dist = _mm256_setr_epi64x(pqDataA->deltaDist[0], pqDataA->deltaDist[0], pqDataA->deltaDist[3], pqDataA->deltaDist[3]);
  __m256i pq_b_delta_dist = _mm256_setr_epi64x(pqDataA->deltaDist[2], pqDataA->deltaDist[2], pqDataA->deltaDist[1], pqDataA->deltaDist[1]);

  __m256i rd_cost_a = _mm256_load_si256((__m256i const*)&state->m_rdCost[start]);
  __m256i rd_cost_b = rd_cost_a;
  __m256i rd_cost_z = rd_cost_a;

  rd_cost_a = _mm256_add_epi64(rd_cost_a, pq_a_delta_dist);
  rd_cost_b = _mm256_add_epi64(rd_cost_b, pq_b_delta_dist);


  if (state->all_gte_four) {
    // pqDataA
    // In case the both levels are smaller than 4 or gte 4 avx 2 can be used
    if (pqDataA->absLevel[0] < 4 && pqDataA->absLevel[3] < 4) {
      // The coeffFracBits arrays are 6 elements long, so we need to offset the indices and gather is only eficient way to load the data
      __m128i offsets = _mm_set_epi32(18 + pqDataA->absLevel[3], 12 + pqDataA->absLevel[3], 6 + pqDataA->absLevel[0], 0 + pqDataA->absLevel[0]);
      __m128i coeff_frac_bits = _mm_i32gather_epi32(&state->m_coeffFracBits[start][0], offsets, 4);
      // RD costs are 64 bit, so we need to extend the 32 bit values
      __m256i ext_frac_bits = _mm256_cvtepi32_epi64(coeff_frac_bits);
      rd_cost_a = _mm256_add_epi64(rd_cost_a, ext_frac_bits);
    }

    else if (pqDataA->absLevel[0] >= 4 && pqDataA->absLevel[3] >= 4) {
      __m128i value = _mm_set_epi32((pqDataA->absLevel[3] - 4) >> 1, (pqDataA->absLevel[3] - 4) >> 1, (pqDataA->absLevel[0] - 4) >> 1, (pqDataA->absLevel[0] - 4) >> 1);

      __m128i offsets = _mm_set_epi32(18 + pqDataA->absLevel[3], 12 + pqDataA->absLevel[3], 6 + pqDataA->absLevel[0], 0 + pqDataA->absLevel[0]);
      __m128i t = _mm_slli_epi32(value, 1);
      offsets = _mm_sub_epi32(offsets, t);
      __m128i coeff_frac_bits = _mm_i32gather_epi32(state->m_coeffFracBits[start], offsets, 4);

      __m128i max_rice = _mm_set1_epi32(31);
      value = _mm_min_epi32(value, max_rice);
      // In the original implementation the goRiceTab is selected beforehand, but since we need to load from
      // potentially four different locations, we need to calculate the offsets and use gather
      __m128i go_rice_tab = _mm_cvtepi8_epi32(_mm_loadu_si128((__m128i *)&state->m_goRicePar[start]));
      go_rice_tab = _mm_slli_epi32(go_rice_tab, 5);
      value = _mm_add_epi32(value, go_rice_tab);

      __m128i temp = _mm_add_epi32(coeff_frac_bits, _mm_i32gather_epi32(&g_goRiceBits[0][0], value, 4));
      rd_cost_a = _mm256_add_epi64(rd_cost_a, _mm256_cvtepi32_epi64(temp));
    } else {
      const int pqAs[4] = {0, 0, 3, 3};
      ALIGNED(32) int64_t rd_costs[4] = {0, 0, 0, 0};
      // AVX2 cannot be used so we have to loop the values normally
      for (int i = 0; i < 4; i++) {
        const int      state_offset = start + i;
        const int      pqA = pqAs[i];
        const int32_t* goRiceTab = g_goRiceBits[state->m_goRicePar[state_offset]];
        if (pqDataA->absLevel[pqA] < 4) {
          rd_costs[i] = state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqA]];
        } else {
          const coeff_t value = (pqDataA->absLevel[pqA] - 4) >> 1;
          rd_costs[i] += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqA] - (value << 1)] + goRiceTab[value < RICEMAX ? value : RICEMAX - 1];
        }
      }
      rd_cost_a = _mm256_add_epi64(rd_cost_a, _mm256_loadu_si256((__m256i const *)&rd_costs[0]));
    }

    // pqDataB, same stuff as for pqDataA
    if (pqDataA->absLevel[1] < 4 && pqDataA->absLevel[2] < 4) {
      __m128i offsets = _mm_set_epi32(18 + pqDataA->absLevel[1], 12 + pqDataA->absLevel[1], 6 + pqDataA->absLevel[2], 0 + pqDataA->absLevel[2]);
      __m128i coeff_frac_bits = _mm_i32gather_epi32(state->m_coeffFracBits[start], offsets, 4);
      __m256i ext_frac_bits = _mm256_cvtepi32_epi64(coeff_frac_bits);
      rd_cost_b = _mm256_add_epi64(rd_cost_b, ext_frac_bits);
    } else if (pqDataA->absLevel[1] >= 4 && pqDataA->absLevel[2] >= 4) {
      __m128i value = _mm_set_epi32((pqDataA->absLevel[1] - 4) >> 1, (pqDataA->absLevel[1] - 4) >> 1, (pqDataA->absLevel[2] - 4) >> 1, (pqDataA->absLevel[2] - 4) >> 1);

      __m128i offsets = _mm_set_epi32(18 + pqDataA->absLevel[1], 12 + pqDataA->absLevel[1], 6 + pqDataA->absLevel[2], 0 + pqDataA->absLevel[2]);
      __m128i t = _mm_slli_epi32(value, 1);
      offsets = _mm_sub_epi32(offsets, t);
      __m128i coeff_frac_bits = _mm_i32gather_epi32(state->m_coeffFracBits[start], offsets, 4);

      __m128i max_rice = _mm_set1_epi32(31);
      value = _mm_min_epi32(value, max_rice);
      __m128i go_rice_tab = _mm_cvtepi8_epi32(_mm_loadu_si128((__m128i*)&state->m_goRicePar[start]));
      go_rice_tab = _mm_slli_epi32(go_rice_tab, 5);
      value = _mm_add_epi32(value, go_rice_tab);

      __m128i temp = _mm_add_epi32(coeff_frac_bits, _mm_i32gather_epi32(&g_goRiceBits[0][0], value, 4));
      rd_cost_b = _mm256_add_epi64(rd_cost_b, _mm256_cvtepi32_epi64(temp));
    } else {
      const int pqBs[4] = {2, 2, 1, 1};
      int64_t rd_costs[4] = {0, 0, 0, 0}; 
      for (int i = 0; i < 4; i++) {
        const int      state_offset = start + i;
        const int      pqB = pqBs[i];
        const int32_t* goRiceTab = g_goRiceBits[state->m_goRicePar[state_offset]];
        if (pqDataA->absLevel[pqB] < 4) {
          rd_costs[i] = state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqB]];
        } else {
          const coeff_t value = (pqDataA->absLevel[pqB] - 4) >> 1;
          rd_costs[i] += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqB] - (value << 1)] + goRiceTab[value < RICEMAX ? value : RICEMAX - 1];
        }
      }
      rd_cost_b =
        _mm256_add_epi64(rd_cost_b, _mm256_loadu_si256((__m256i const *) & rd_costs[0]));
    }

    if (spt == SCAN_ISCSBB) {
      // This loads values such as that the values are
      // |State 0 Flag 0|State 0 Flag 1|State 1 Flag 0|State 1 Flag 1|State 2 Flag 0|State 2 Flag 1|State 3 Flag 0|State 3 Flag 1|
      // By setting the flag 1 bits to zero we get the flag 0 values as 64 bit integers (even) variable which we can be summed to the rd_cost
      // Flag 1 values can be shifted 32 to right and again we have 64 bit integeres holding the values (odd) which can be summed to the rd_cost
      __m256i original = _mm256_loadu_si256((__m256i const*)state->m_sigFracBits[start]);
      __m256i even      = _mm256_and_si256(original, _mm256_set1_epi64x(0xffffffff));
      __m256i odd  = _mm256_srli_epi64(original, 32);
      rd_cost_a = _mm256_add_epi64(rd_cost_a, odd);
      rd_cost_b = _mm256_add_epi64(rd_cost_b, odd);
      rd_cost_z = _mm256_add_epi64(rd_cost_z, even);
    } else if (spt == SCAN_SOCSBB) {
      __m256i original = _mm256_loadu_si256((__m256i const*)state->m_sigFracBits[start]);

      // Same here
      __m256i m_sigFracBits_0 = _mm256_and_si256(original, _mm256_set1_epi64x(0xffffffff));
      __m256i m_sigFracBits_1 = _mm256_srli_epi64(original, 32);

      original = _mm256_loadu_si256((__m256i const*)state->m_sbbFracBits[start]);
      __m256i m_sbbFracBits_1 = _mm256_srli_epi64(original, 32);

      rd_cost_a = _mm256_add_epi64(rd_cost_a, m_sbbFracBits_1);
      rd_cost_b = _mm256_add_epi64(rd_cost_b, m_sbbFracBits_1);
      rd_cost_z = _mm256_add_epi64(rd_cost_z, m_sbbFracBits_1);

      rd_cost_a = _mm256_add_epi64(rd_cost_a, m_sigFracBits_1);
      rd_cost_b = _mm256_add_epi64(rd_cost_b, m_sigFracBits_1);
      rd_cost_z = _mm256_add_epi64(rd_cost_z, m_sigFracBits_0);
    }
    else {
      int num_sig_sbb;
      memcpy(&num_sig_sbb, &state->m_numSigSbb[start], 4);
      // numSigSbb only has values 1 or zero, so if all 4 values are 1 the complete value is 0x01010101
      if (num_sig_sbb == 0x01010101) {
        __m256i original = _mm256_loadu_si256((__m256i const*)state->m_sigFracBits[start]);
      __m256i even      = _mm256_and_si256(original, _mm256_set1_epi64x(0xffffffff));
      __m256i odd  = _mm256_srli_epi64(original, 32);
        rd_cost_a = _mm256_add_epi64(rd_cost_a, odd);
        rd_cost_b = _mm256_add_epi64(rd_cost_b, odd);
        rd_cost_z = _mm256_add_epi64(rd_cost_z, even);     
      }
      else if (num_sig_sbb == 0) {
        rd_cost_z = _mm256_setr_epi64x(decisions->rdCost[0], decisions->rdCost[0], decisions->rdCost[3], decisions->rdCost[3]);
      }

      else {
        const int ALIGNED(32) pqAs[4] = {0, 0, 3, 3};
        _mm256_store_si256((__m256i*)temp_rd_cost_a, rd_cost_a);
        _mm256_store_si256((__m256i*)temp_rd_cost_b, rd_cost_b);
        _mm256_store_si256((__m256i*)temp_rd_cost_z, rd_cost_z);
        for (int i = 0; i < 4; i++) {
          const int state_offset = start + i;
          if (state->m_numSigSbb[state_offset]) {
            temp_rd_cost_a[i] += state->m_sigFracBits[state_offset][1];
            temp_rd_cost_b[i] += state->m_sigFracBits[state_offset][1];
            temp_rd_cost_z[i] += state->m_sigFracBits[state_offset][0];
          } else {
            temp_rd_cost_z[i] = decisions->rdCost[pqAs[i]];
          }
        }
        rd_cost_a = _mm256_loadu_si256((__m256i*)temp_rd_cost_a);
        rd_cost_b = _mm256_loadu_si256((__m256i*)temp_rd_cost_b);
        rd_cost_z = _mm256_loadu_si256((__m256i*)temp_rd_cost_z);
      }
    }
  } else if (state->all_lt_four) {
    __m128i scale_bits = _mm_set1_epi32(1 << SCALE_BITS);
    __m128i max_rice = _mm_set1_epi32(31);
    __m128i go_rice_zero = _mm_cvtepi8_epi32(_mm_loadu_si128((const __m128i*)&state->m_goRiceZero[start]));
    // RD cost A
    {
      __m128i pq_abs_a = _mm_set_epi32(pqDataA->absLevel[3], pqDataA->absLevel[3], pqDataA->absLevel[0], pqDataA->absLevel[0]);
      // Calculate mask for pqDataA->absLevel <= state->m_goRiceZero
      // The mask is reverse of the one that is used in the scalar code so the values are in other order in blendv
      __m128i cmp = _mm_cmpgt_epi32(pq_abs_a, go_rice_zero);

      // pqDataA->absLevel < RICEMAX ? pqDataA->absLevel : RICEMAX - 1
      __m128i go_rice_smaller = _mm_min_epi32(pq_abs_a, max_rice);

      // pqDataA->absLevel - 1
      __m128i other = _mm_sub_epi32(pq_abs_a, _mm_set1_epi32(1));

      __m128i selected = _mm_blendv_epi8(other, go_rice_smaller, cmp);

      // Again calculate the offset for the different go_rice_tabs
      __m128i go_rice_offset = _mm_cvtepi8_epi32(_mm_loadu_si128((__m128i*)&state->m_goRicePar[start]));
      go_rice_offset = _mm_slli_epi32(go_rice_offset, 5);

      __m128i offsets = _mm_add_epi32(selected, go_rice_offset);
      __m128i go_rice_tab = _mm_i32gather_epi32(&g_goRiceBits[0][0], offsets, 4);
      //(1 << SCALE_BITS) + goRiceTab[selected]
      __m128i temp = _mm_add_epi32(go_rice_tab, scale_bits);

      rd_cost_a = _mm256_add_epi64(rd_cost_a, _mm256_cvtepi32_epi64(temp));
    }
    // RD cost b, same as RD cost A
    {
      __m128i pq_abs_b = _mm_set_epi32(pqDataA->absLevel[1], pqDataA->absLevel[1], pqDataA->absLevel[2], pqDataA->absLevel[2]);
      __m128i cmp = _mm_cmpgt_epi32(pq_abs_b, go_rice_zero);

      __m128i go_rice_smaller = _mm_min_epi32(pq_abs_b, max_rice);

      __m128i other = _mm_sub_epi32(pq_abs_b, _mm_set1_epi32(1));

      __m128i selected = _mm_blendv_epi8(other, go_rice_smaller, cmp);


      __m128i go_rice_offset = _mm_cvtepi8_epi32(_mm_loadu_si128((__m128i*)&state->m_goRicePar[start]));
      go_rice_offset = _mm_slli_epi32(go_rice_offset, 5);

      __m128i offsets = _mm_add_epi32(selected, go_rice_offset);
      __m128i go_rice_tab = _mm_i32gather_epi32(&g_goRiceBits[0][0], offsets, 4);
      __m128i temp = _mm_add_epi32(go_rice_tab, scale_bits);

      rd_cost_b = _mm256_add_epi64(rd_cost_b, _mm256_cvtepi32_epi64(temp));
    }
    // RD cost Z
    {
      // This time the go_rice_tab is offset with only the go_rize_zero
      __m128i go_rice_offset = _mm_cvtepi8_epi32(_mm_loadu_si128((__m128i*)&state->m_goRicePar[start]));
      go_rice_offset = _mm_slli_epi32(go_rice_offset, 5);

      go_rice_offset = _mm_add_epi32(go_rice_offset, go_rice_zero);
      __m128i go_rice_tab = _mm_i32gather_epi32(&g_goRiceBits[0][0], go_rice_offset, 4);
      rd_cost_z = _mm256_add_epi64(rd_cost_z, _mm256_cvtepi32_epi64(go_rice_tab));
    }
  } else {
    const int pqAs[4] = {0, 0, 3, 3};
    const int pqBs[4] = {2, 2, 1, 1};
    const int decision_a[4] = {0, 2, 1, 3};
    for (int i = 0; i < 4; i++) {
      const int      state_offset = start + i;
      const int32_t* goRiceTab = g_goRiceBits[state->m_goRicePar[state_offset]];
      const int pqA = pqAs[i];
      const int pqB = pqBs[i];
      int64_t rdCostA = state->m_rdCost[state_offset] + pqDataA->deltaDist[pqA];
      int64_t rdCostB = state->m_rdCost[state_offset] + pqDataA->deltaDist[pqB];
      int64_t rdCostZ = state->m_rdCost[state_offset];
      if (state->m_remRegBins[state_offset] >= 4) {
        if (pqDataA->absLevel[pqA] < 4) {
          rdCostA += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqA]];
        } else {
          const coeff_t value = (pqDataA->absLevel[pqA] - 4) >> 1;
          rdCostA += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqA] - (value << 1)] + goRiceTab[value < RICEMAX ? value : RICEMAX - 1];
        }
        if (pqDataA->absLevel[pqB] < 4) {
          rdCostB += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqB]];
        } else {
          const coeff_t value = (pqDataA->absLevel[pqB] - 4) >> 1;
          rdCostB += state->m_coeffFracBits[state_offset][pqDataA->absLevel[pqB] - (value << 1)] + goRiceTab[value < RICEMAX ? value : RICEMAX - 1];
        }
        if (spt == SCAN_ISCSBB) {
          rdCostA += state->m_sigFracBits[state_offset][1];
          rdCostB += state->m_sigFracBits[state_offset][1];
          rdCostZ += state->m_sigFracBits[state_offset][0];
        } else if (spt == SCAN_SOCSBB) {
          rdCostA += state->m_sbbFracBits[state_offset][1] + state->m_sigFracBits[state_offset][1];
          rdCostB += state->m_sbbFracBits[state_offset][1] + state->m_sigFracBits[state_offset][1];
          rdCostZ += state->m_sbbFracBits[state_offset][1] + state->m_sigFracBits[state_offset][0];
        } else if (state->m_numSigSbb[state_offset]) {
          rdCostA += state->m_sigFracBits[state_offset][1];
          rdCostB += state->m_sigFracBits[state_offset][1];
          rdCostZ += state->m_sigFracBits[state_offset][0];
        } else {
          rdCostZ = decisions->rdCost[decision_a[i]];
        }
      } else {
        rdCostA += (1 << SCALE_BITS) + goRiceTab[pqDataA->absLevel[pqA] <= state->m_goRiceZero[state_offset] ? pqDataA->absLevel[pqA] - 1 : (pqDataA->absLevel[pqA] < RICEMAX ? pqDataA->absLevel[pqA] : RICEMAX - 1)];
        rdCostB += (1 << SCALE_BITS) + goRiceTab[pqDataA->absLevel[pqB] <= state->m_goRiceZero[state_offset] ? pqDataA->absLevel[pqB] - 1 : (pqDataA->absLevel[pqB] < RICEMAX ? pqDataA->absLevel[pqB] : RICEMAX - 1)];
        rdCostZ += goRiceTab[state->m_goRiceZero[state_offset]];
      }
      temp_rd_cost_a[i] = rdCostA;
      temp_rd_cost_b[i] = rdCostB;
      temp_rd_cost_z[i] = rdCostZ;
    }
    rd_cost_a = _mm256_loadu_si256((__m256i*)temp_rd_cost_a);
    rd_cost_b = _mm256_loadu_si256((__m256i*)temp_rd_cost_b);
    rd_cost_z = _mm256_loadu_si256((__m256i*)temp_rd_cost_z);
  }
  // Re order the cost so that cost of state 0 is in the first element state 1 in second etc
  rd_cost_a = _mm256_permute4x64_epi64(rd_cost_a, 216);
  rd_cost_b = _mm256_permute4x64_epi64(rd_cost_b, 141);
  rd_cost_z = _mm256_permute4x64_epi64(rd_cost_z, 216);
  __m256i rd_cost_decision = _mm256_load_si256((__m256i*)decisions->rdCost);

  __m256i decision_abs_coeff = _mm256_load_si256((__m256i*)decisions->absLevel);
  __m256i decision_prev_state = _mm256_load_si256((__m256i*)decisions->prevId);
  __m256i decision_data = _mm256_permute2x128_si256(decision_abs_coeff, decision_prev_state, 0x20);
  __m256i mask = _mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0);

  // Store data for all of the cost so that the lower 32 bits have coefficient magnitude and upper have the previous state
  decision_data = _mm256_permutevar8x32_epi32(decision_data, mask);
  __m256i a_data = _mm256_set_epi32(3, pqDataA->absLevel[3], 1, pqDataA->absLevel[0], 2, pqDataA->absLevel[3], 0, pqDataA->absLevel[0]);
  __m256i b_data = _mm256_set_epi32(2, pqDataA->absLevel[1], 0, pqDataA->absLevel[2], 3, pqDataA->absLevel[1], 1, pqDataA->absLevel[2]);
  __m256i z_data = _mm256_set_epi32(3, 0, 1, 0, 2, 0, 0, 0);

  __m256i a_vs_b = _mm256_cmpgt_epi64(rd_cost_a, rd_cost_b);
  __m256i cheaper_first = _mm256_blendv_epi8(rd_cost_a, rd_cost_b, a_vs_b);
  __m256i cheaper_first_data = _mm256_blendv_epi8(a_data, b_data, a_vs_b);

  __m256i z_vs_decision = _mm256_cmpgt_epi64(rd_cost_decision, rd_cost_z);
  __m256i cheaper_second = _mm256_blendv_epi8(rd_cost_decision, rd_cost_z, z_vs_decision);
  __m256i cheaper_second_data = _mm256_blendv_epi8(decision_data, z_data, z_vs_decision);

  __m256i final_decision = _mm256_cmpgt_epi64(cheaper_second, cheaper_first);
  __m256i final_rd_cost = _mm256_blendv_epi8(cheaper_second, cheaper_first, final_decision);
  __m256i final_data = _mm256_blendv_epi8(cheaper_second_data, cheaper_first_data, final_decision);

  _mm256_store_si256((__m256i*)decisions->rdCost, final_rd_cost);
  final_data = _mm256_permutevar8x32_epi32(final_data, _mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0));
  _mm256_storeu2_m128i((__m128i *)decisions->prevId, (__m128i *)decisions->absLevel, final_data);
}


static INLINE void checkRdCostSkipSbbZeroOut(
  Decision*                        decision,
  const all_depquant_states* const state,
  int                              decision_id,
  int                              skip_offset)
{
  int64_t rdCost = state->m_rdCost[decision_id + skip_offset] + state->m_sbbFracBits[decision_id + skip_offset][0];
  decision->rdCost[decision_id] = rdCost;
  decision->absLevel[decision_id] = 0;
  decision->prevId[decision_id] = 4 + state->m_stateId[decision_id + skip_offset];
}


static INLINE void checkRdCostSkipSbb(const all_depquant_states* const state, Decision * decisions, int decision_id, int skip_offset)
{
  int64_t rdCost = state->m_rdCost[skip_offset + decision_id] + state->m_sbbFracBits[skip_offset + decision_id][0];
  if (rdCost < decisions->rdCost[decision_id])
  {
    decisions->rdCost[decision_id] = rdCost;
    decisions->absLevel[decision_id] = 0;
    decisions->prevId[decision_id] = 4 + state->m_stateId[skip_offset + decision_id];
  }
}

static INLINE void checkRdCostStart(const depquant_state* const state, int32_t lastOffset, const PQData *pqData, Decision *decisions, int
                                    decision_id)
{
  int64_t rdCost = pqData->deltaDist[decision_id] + lastOffset;
  if (pqData->absLevel[decision_id] < 4) {
    rdCost += state->m_coeffFracBits[pqData->absLevel[decision_id]];
  }
  else {
    const coeff_t value = (pqData->absLevel[decision_id] - 4) >> 1;
    rdCost += state->m_coeffFracBits[pqData->absLevel[decision_id] - (value << 1)]
              + g_goRiceBits[state->m_goRicePar][value < RICEMAX ? value : RICEMAX - 1];
  }
  if (rdCost < decisions->rdCost[decision_id]) {
    decisions->rdCost[decision_id] = rdCost;
    decisions->absLevel[decision_id] = pqData->absLevel[decision_id];
    decisions->prevId[decision_id] = -1;
  }
}

static INLINE void preQuantCoeff(const quant_block * const qp, const coeff_t absCoeff, PQData* pqData, coeff_t quanCoeff)
{
  int64_t scaledOrg = (int64_t)(absCoeff) * quanCoeff;
  coeff_t  qIdx = MAX(1, (coeff_t)MIN(qp->m_maxQIdx, ((scaledOrg + qp->m_QAdd) >> qp->m_QShift)));
  int64_t scaledAdd = qIdx * qp->m_DistStepAdd - scaledOrg * qp->m_DistOrgFact;
  int index = qIdx & 3;
  pqData->deltaDist[index] = (scaledAdd * qIdx + qp->m_DistAdd) >> qp->m_DistShift;
  pqData->absLevel[index] = (++qIdx) >> 1;
  scaledAdd += qp->m_DistStepAdd;
  index = qIdx & 3;
  pqData->deltaDist[index] = (scaledAdd * qIdx + qp->m_DistAdd) >> qp->m_DistShift;
  pqData->absLevel[index] = (++qIdx) >> 1;
  scaledAdd += qp->m_DistStepAdd;
  index = qIdx & 3;
  pqData->deltaDist[index] = (scaledAdd * qIdx + qp->m_DistAdd) >> qp->m_DistShift;
  pqData->absLevel[index] = (++qIdx) >> 1;
  scaledAdd += qp->m_DistStepAdd;
  index = qIdx & 3;
  pqData->deltaDist[index] = (scaledAdd * qIdx + qp->m_DistAdd) >> qp->m_DistShift;
  pqData->absLevel[index] = (++qIdx) >> 1;
}


static const Decision startDec = { .rdCost = {INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2, INT64_MAX >> 2},
  .absLevel = {-1, -1, -1, -1, 0, 0, 0, 0}, .prevId = {-2, -2, -2, -2, 4, 5, 6, 7} };


static void xDecide(
  all_depquant_states* const all_states,
  depquant_state* const      m_startState,
  quant_block *              qp,
  const enum ScanPosType     spt,
  const coeff_t              absCoeff,
  const int                  lastOffset,
  Decision*                  decisions,
  bool                       zeroOut,
  coeff_t                    quanCoeff,
  const int                  skip_offset,
  const int                  prev_offset)
{
  memcpy(decisions, &startDec, sizeof(Decision));

  if (zeroOut) {
    if (spt == SCAN_EOCSBB) {
      checkRdCostSkipSbbZeroOut(decisions, all_states, 0, skip_offset);
      checkRdCostSkipSbbZeroOut(decisions, all_states, 1, skip_offset);
      checkRdCostSkipSbbZeroOut(decisions, all_states, 2, skip_offset);
      checkRdCostSkipSbbZeroOut(decisions, all_states, 3, skip_offset);
    }
    return;
  }

  PQData pqData;
  preQuantCoeff(qp, absCoeff, &pqData, quanCoeff);
  check_rd_costs_avx2(all_states, spt, &pqData, decisions, prev_offset);
  if (spt == SCAN_EOCSBB) {
    checkRdCostSkipSbb(all_states, decisions, 0, skip_offset);
    checkRdCostSkipSbb(all_states, decisions, 1, skip_offset);
    checkRdCostSkipSbb(all_states, decisions, 2, skip_offset);
    checkRdCostSkipSbb(all_states, decisions, 3, skip_offset);
  }

  checkRdCostStart(m_startState, lastOffset, &pqData, decisions, 0);
  checkRdCostStart(m_startState, lastOffset, &pqData, decisions, 2);
}


static void update_state_eos_avx2(context_store* ctxs, const uint32_t scan_pos, const uint32_t cg_pos,
                                  const uint32_t sigCtxOffsetNext, const uint32_t gtxCtxOffsetNext,
                                  const uint32_t width_in_sbb, const uint32_t height_in_sbb,
                                  const uint32_t next_sbb_right, const uint32_t next_sbb_below,
                                  const Decision* decisions)
{
  all_depquant_states* state = &ctxs->m_allStates;
  bool all_above_minus_two = true;
  bool all_between_zero_and_three = true;
  bool all_above_four = true;

  
  int state_offset = ctxs->m_curr_state_offset;
  __m256i rd_cost = _mm256_load_si256((__m256i const*)decisions->rdCost);
  _mm256_store_si256((__m256i *)& ctxs->m_allStates.m_rdCost[state_offset], rd_cost);
  for (int i = 0; i < 4; ++i) {
    all_above_minus_two &= decisions->prevId[i] > -2;
    all_between_zero_and_three &= decisions->prevId[i] >= 0 && decisions->prevId[i] < 4;
    all_above_four &= decisions->prevId[i] >= 4;
  }
  if (all_above_minus_two) {
    bool all_have_previous_state = true;
    __m128i prev_state;
    __m128i prev_state_no_offset;
    __m128i abs_level = _mm_load_si128((const __m128i*)decisions->absLevel);
    __m128i control = _mm_setr_epi8(0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12);
    if (all_above_four) {
      prev_state = _mm_set1_epi32(ctxs->m_skip_state_offset);
      prev_state_no_offset = _mm_sub_epi32(_mm_load_si128((const __m128i*)decisions->prevId), _mm_set1_epi32(4));
      prev_state = _mm_add_epi32(
        prev_state,
            prev_state_no_offset
      );
      memset(&state->m_numSigSbb[state_offset], 0, 4);
      memset(state->m_absLevels[state_offset >> 2], 0, 64 * sizeof(uint8_t));    
      
    } else if (all_between_zero_and_three) {
      prev_state_no_offset = _mm_load_si128((const __m128i*)decisions->prevId);
      prev_state = _mm_add_epi32(
        prev_state_no_offset,
        _mm_set1_epi32(ctxs->m_prev_state_offset)
      );
      // Set the high bytes to 0xff so that the shuffle will set them to zero and it won't cause problems with the min_epi32
      __m128i prev_state_with_ff_high_bytes = _mm_or_si128(prev_state, _mm_set1_epi32(0xffffff00));
      __m128i num_sig_sbb = _mm_load_si128((const __m128i*)state->m_numSigSbb);
      num_sig_sbb = _mm_shuffle_epi8(num_sig_sbb, prev_state_with_ff_high_bytes);
      num_sig_sbb = _mm_or_si128(
        num_sig_sbb,
        _mm_min_epi32(abs_level, _mm_set1_epi32(1))
      );

      num_sig_sbb = _mm_shuffle_epi8(num_sig_sbb, control);
      int num_sig_sbb_s = _mm_extract_epi32(num_sig_sbb, 0);
      memcpy(&state->m_numSigSbb[state_offset], &num_sig_sbb_s, 4);

      // Set this so that the temp_prev_state has the previous state set into the first 4 bytes and duplicated to the second 4 bytes
      __m128i temp_prev_state = _mm_shuffle_epi8(prev_state_no_offset, control);
      __m256i prev_state_256 = _mm256_castsi128_si256(temp_prev_state);
      // Duplicate the state all over the vector so that all 32 bytes hold the previous states
      prev_state_256 = _mm256_permute4x64_epi64(prev_state_256, 0);
      // Increment the second set by four, third by eight and fourth by twelve and repeat for the second lane
      __m256i temp_add = _mm256_setr_epi32(0, 0x04040404, 0x08080808, 0x0c0c0c0c, 0, 0x04040404, 0x08080808, 0x0c0c0c0c);
      prev_state_256 = _mm256_add_epi8(prev_state_256, temp_add);
      for (int i = 0; i < 64; i += (256 / (8 * sizeof(uint8_t)))) {
        __m256i data = _mm256_load_si256((__m256i*)&state->m_absLevels[ctxs->m_prev_state_offset >> 2][i]);
        data = _mm256_shuffle_epi8(data, prev_state_256);
        _mm256_store_si256((__m256i*)&state->m_absLevels[ctxs->m_curr_state_offset >> 2][i], data);
      }
    } else {
      // TODO: it would be possible to do the absLevels update with avx2 even here just would need to set the shuffle mask to
      // 0xff for the states that don't have previous state or the previous state is a skip state
      int prev_state_s[4] = {-1, -1, -1, -1};
      for (int i = 0; i < 4; ++i) {
        const int decision_id = i;
        const int curr_state_offset = state_offset + i;
        if (decisions->prevId[decision_id] >= 4) {
          prev_state_s[i] = ctxs->m_skip_state_offset + (decisions->prevId[decision_id] - 4);
          state->m_numSigSbb[curr_state_offset] = 0;
          for (int j = i; j < 64; j += 4) {
            state->m_absLevels[curr_state_offset >> 2][j] = 0;
          }
        } else if (decisions->prevId[decision_id] >= 0) {
          prev_state_s[i] = ctxs->m_prev_state_offset + decisions->prevId[decision_id];
          state->m_numSigSbb[curr_state_offset] = state->m_numSigSbb[prev_state_s[i]] || !!decisions->absLevel[decision_id];
          for (int j = 0; j < 64; j += 4) {
            state->m_absLevels[curr_state_offset >> 2][j + i] = state->m_absLevels[ctxs->m_prev_state_offset >> 2][j + decisions->prevId[decision_id]];
          }
        } else {
          state->m_numSigSbb[curr_state_offset] = 1;
          for (int j = i; j < 64; j += 4) {
            state->m_absLevels[curr_state_offset >> 2][j] = 0;
          }
          all_have_previous_state = false;
        }
      }
      prev_state = _mm_loadu_si128((__m128i const*)prev_state_s);
    }
    uint32_t level_offset = scan_pos & 15;
    __m128i  max_abs = _mm_min_epi32(abs_level, _mm_set1_epi32(255));
    max_abs = _mm_shuffle_epi8(max_abs, control);
    uint32_t packed_max_abs = _mm_extract_epi32(max_abs, 0);
    memcpy(&state->m_absLevels[state_offset >> 2][level_offset * 4], &packed_max_abs, 4);


    // Update common context
    __m128i last;
    {
      const uint32_t numSbb = width_in_sbb * height_in_sbb;
      common_context* cc = &ctxs->m_common_context;
      size_t   setCpSize = cc->m_nbInfo[scan_pos - 1].maxDist * sizeof(uint8_t);
      uint8_t* sbbFlags  = cc->m_allSbbCtx[cc->m_curr_sbb_ctx_offset].sbbFlags;
      uint8_t* levels   = cc->m_allSbbCtx[cc->m_curr_sbb_ctx_offset].levels + scan_pos * 4;
      uint8_t* levels_in = cc->m_allSbbCtx[cc->m_prev_sbb_ctx_offset].levels + scan_pos * 4;
      int      previous_state_array[4];
      _mm_storeu_si128((__m128i*)previous_state_array, prev_state);

      if (all_have_previous_state) {
        __m128i temp_p_state = _mm_shuffle_epi8(prev_state, control);
        // Similarly to how the abs level was done earlier set the previous state duplicated across the lane
        __m128i ref_sbb_ctx_offset = _mm_load_si128((__m128i*)ctxs->m_allStates.m_refSbbCtxId);
        ref_sbb_ctx_offset = _mm_shuffle_epi8(ref_sbb_ctx_offset, temp_p_state);
        // numSbb is two or four, in case it is one this function is never called
        if (numSbb <= 4) {
          __m128i incremented_ref_sbb_ctx_offset = _mm_add_epi8(
            ref_sbb_ctx_offset,
            _mm_setr_epi8(0, 0, 0, 0, 4, 4, 4, 4, 8, 8, 8, 8, 12, 12, 12, 12)
          );
          // In case the ref_sbb_ctx is minus one the values need to be set to zero, which is achieved by
          // first finding which states have the minus one and then the blend is used after the load to
          // set the corresponding values to zero
          __m128i blend_mask = _mm_cmpeq_epi8(ref_sbb_ctx_offset, _mm_set1_epi32(0xffffffff));
          __m128i sbb_flags = _mm_loadu_si128((__m128i*)cc->m_allSbbCtx[cc->m_prev_sbb_ctx_offset].sbbFlags);
          sbb_flags = _mm_shuffle_epi8(sbb_flags, incremented_ref_sbb_ctx_offset);
          sbb_flags = _mm_blendv_epi8(sbb_flags, _mm_set1_epi64x(0), blend_mask);
          if (numSbb == 2) {
            uint64_t temp = _mm_extract_epi64(sbb_flags, 0);
            memcpy(sbbFlags, &temp, 8);
          } else {
            _mm_storeu_si128((__m128i*)sbbFlags, sbb_flags);
          }
        } else {
          __m256i extended_ref_state = _mm256_zextsi128_si256(ref_sbb_ctx_offset);
          extended_ref_state = _mm256_permute4x64_epi64(extended_ref_state, 0);
          __m256i inc_ref_state = _mm256_add_epi8(
            extended_ref_state,
            _mm256_setr_epi32(0, 0x04040404, 0x08080808, 0x0c0c0c0c,0, 0x04040404, 0x08080808, 0x0c0c0c0c)
          );
          // Unlike the case for two or four sbb, the blendv is used to set the shuffle mask to -1 so that
          // the shuffle will set the values to zero. Its better to do this way here so that the blendv is
          // not called in the loop, and the other is done the otherway because I implemented it first
          // and only realized afterwards that this order is better
          __m256i blend_mask = _mm256_cmpeq_epi8(extended_ref_state, _mm256_set1_epi32(0xffffffff));
          inc_ref_state = _mm256_blendv_epi8(inc_ref_state, _mm256_set1_epi32(0xffffffff), blend_mask);
          for (int i = 0; i < numSbb * 4; i += 32) {
            __m256i sbb_flags = _mm256_loadu_si256((__m256i*)(&cc->m_allSbbCtx[cc->m_prev_sbb_ctx_offset].sbbFlags[i]));
            sbb_flags = _mm256_shuffle_epi8(sbb_flags, inc_ref_state);
            _mm256_store_si256((__m256i*)&sbbFlags[i], sbb_flags);
          }
        }
        // The first 16 variables will be loaded from the previous state so this can be started from 16
        int levels_start = 16;
        // Do avx2 optimized version for the amount that is divisible by 8 (four states of 8 1-byte values)
        const uint64_t limit        = setCpSize & ~(8 - 1);
        if (levels_start < limit) {
          // Overall this is the same to the numSbb > 4
          __m256i extended_ref_state = _mm256_zextsi128_si256(ref_sbb_ctx_offset);
          extended_ref_state = _mm256_permute4x64_epi64(extended_ref_state, 0);
          __m256i inc_ref_state = _mm256_add_epi8(
            extended_ref_state,
            _mm256_setr_epi32(0, 0x04040404, 0x08080808, 0x0c0c0c0c,0, 0x04040404, 0x08080808, 0x0c0c0c0c)
          );
          __m256i blend_mask = _mm256_cmpeq_epi8(extended_ref_state, _mm256_set1_epi32(0xffffffff));
          inc_ref_state = _mm256_blendv_epi8(inc_ref_state, _mm256_set1_epi32(0xffffffff), blend_mask);
          for (; levels_start < limit; levels_start += 8) {
             __m256i levels_v = _mm256_loadu_si256((__m256i*)(&levels_in[levels_start * 4]));
            levels_v = _mm256_shuffle_epi8(levels_v, inc_ref_state);
             _mm256_store_si256((__m256i*)&levels[levels_start * 4], levels_v);
          }
        }
        uint8_t ref_sbb[4];
        int     temp_sbb_ref = _mm_extract_epi32(ref_sbb_ctx_offset, 0);
        memcpy(ref_sbb, &temp_sbb_ref, 4);
        // Do the excess that is not divisible by 8
        for (;levels_start < setCpSize; ++levels_start) {
          uint8_t new_values[4];
          new_values[0] = ref_sbb[0] != 0xff ? levels_in[levels_start * 4 + ref_sbb[0]] : 0;
          new_values[1] = ref_sbb[1] != 0xff ? levels_in[levels_start * 4 + ref_sbb[1]] : 0;
          new_values[2] = ref_sbb[2] != 0xff ? levels_in[levels_start * 4 + ref_sbb[2]] : 0;
          new_values[3] = ref_sbb[3] != 0xff ? levels_in[levels_start * 4 + ref_sbb[3]] : 0;
          memcpy(&levels[levels_start * 4], new_values, 4);
        }

      }
      else {
        //TODO: This could also be done using avx2 just need to check for both wheter the previous state
        // is minus one and that if the ref_sbb_ctx_id is minus one. 
        for (int curr_state = 0; curr_state < 4; ++curr_state) {
          const int p_state = previous_state_array[curr_state];
          if (p_state != -1 && ctxs->m_allStates.m_refSbbCtxId[p_state] >= 0) {
            const int prev_sbb = ctxs->m_allStates.m_refSbbCtxId[p_state];
            for (int i = 0; i < numSbb; ++i) {
              sbbFlags[i * 4 + curr_state] = cc->m_allSbbCtx[cc->m_prev_sbb_ctx_offset].sbbFlags[i * 4 + prev_sbb];
            }
            for (int i = 16; i < setCpSize; ++i) {
              levels[i * 4 + curr_state] = levels_in[i * 4 + prev_sbb];
            }
          } else {
            for (int i = 0; i < numSbb; ++i) {
              sbbFlags[i * 4 + curr_state] = 0;
            }
            for (int i = 16; i < setCpSize; ++i) {
              levels[ i * 4 + curr_state] = 0;
            }
          }
        }
      }
      memcpy(levels, ctxs->m_allStates.m_absLevels[state_offset / 4], 64);
      memcpy(&sbbFlags[cg_pos * 4], &ctxs->m_allStates.m_numSigSbb[state_offset], 4);
      
      __m128i sbb_right = next_sbb_right ?
          _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)&cc->m_allSbbCtx[cc->m_curr_sbb_ctx_offset].sbbFlags[next_sbb_right * 4])) :
          _mm_set1_epi32(0);
      
      __m128i sbb_below = next_sbb_below ?
        _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)&cc->m_allSbbCtx[cc->m_curr_sbb_ctx_offset].sbbFlags[next_sbb_below * 4])) :
        _mm_set1_epi32(0);

      __m128i sig_sbb = _mm_or_si128(sbb_right, sbb_below);
      sig_sbb = _mm_min_epi32(sig_sbb, _mm_set1_epi32(1));
      // Gather is not necessary here put it would require at least five operation to do the same thing
      // so the performance gain in my opinion is not worth the readability loss
      __m256i sbb_frac_bits = _mm256_i32gather_epi64((const long long int *)cc->m_sbbFlagBits[0], sig_sbb, 8);
      _mm256_store_si256((__m256i*)state->m_sbbFracBits[state_offset], sbb_frac_bits);

      memset(&state->m_numSigSbb[state_offset], 0, 4);
      memset(&state->m_goRicePar[state_offset], 0, 4);

      uint8_t states[4] = {0, 1, 2, 3};
      memcpy(&state->m_refSbbCtxId[state_offset], states, 4);
      if (all_have_previous_state) {
        __m128i rem_reg_bins = _mm_i32gather_epi32(state->m_remRegBins, prev_state, 4);
        _mm_store_si128((__m128i*) & state->m_remRegBins[state_offset], rem_reg_bins);
      } else {
        const int temp = (state->effWidth * state->effHeight * 28) / 16;
        for (int i = 0; i < 4; ++i) {
          if (previous_state_array[i] != -1) {
            state->m_remRegBins[i + state_offset] = state->m_remRegBins[previous_state_array[i]];
          } else {
            state->m_remRegBins[i + state_offset] = temp;
          }
        }
      }
      
      const int        scanBeg = scan_pos - 16;
      const NbInfoOut* nbOut = cc->m_nbInfo + scanBeg;
      const uint8_t*   absLevels = cc->m_allSbbCtx[cc->m_curr_sbb_ctx_offset].levels + scanBeg * 4;

      __m128i          ones = _mm_set1_epi32(1);
      __m128i         fours = _mm_set1_epi32(4);
      __m256i          all[4];
      uint64_t         temp[4];

      for (int id = 0; id < 16; id++, nbOut++) {
        if (nbOut->num == 0) {
          temp[id % 4] = 0;
          if (id % 4 == 3) {
            all[id / 4] = _mm256_loadu_si256((__m256i const*)temp);
          }
          continue;
        }
        __m128i sum_abs = _mm_set1_epi32(0);
        __m128i sum_abs_1 = _mm_set1_epi32(0);
        __m128i sum_num = _mm_set1_epi32(0);
        switch (nbOut->num) {
        case 5:
          {
            __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&absLevels[nbOut->outPos[4] * 4])));
            sum_abs = _mm_add_epi32(sum_abs, t);
            sum_num   = _mm_add_epi32(sum_num, _mm_min_epi32(t, ones));
            __m128i min_t = _mm_min_epi32(
              t,
              _mm_add_epi32(
                fours,
                _mm_and_si128(t, ones)
              )
            );
            sum_abs_1 = _mm_add_epi32(sum_abs_1, min_t);
          }
        case 4: {
            __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&absLevels[nbOut->outPos[3] * 4])));
            sum_abs = _mm_add_epi32(sum_abs, t);
            sum_num   = _mm_add_epi32(sum_num, _mm_min_epi32(t, ones));
            __m128i min_t = _mm_min_epi32(
              t,
              _mm_add_epi32(
                fours,
                _mm_and_si128(t, ones)));
            sum_abs_1 = _mm_add_epi32(sum_abs_1, min_t);
        }
        case 3: {
            __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&absLevels[nbOut->outPos[2] * 4])));
            sum_abs = _mm_add_epi32(sum_abs, t);
            sum_num   = _mm_add_epi32(sum_num, _mm_min_epi32(t, ones));
            __m128i min_t = _mm_min_epi32(
              t,
              _mm_add_epi32(
                fours,
                _mm_and_si128(t, ones)));
            sum_abs_1 = _mm_add_epi32(sum_abs_1, min_t);
        }
        case 2: {
            __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&absLevels[nbOut->outPos[1] * 4])));
            sum_abs = _mm_add_epi32(sum_abs, t);
            sum_num   = _mm_add_epi32(sum_num, _mm_min_epi32(t, ones));
            __m128i min_t = _mm_min_epi32(
              t,
              _mm_add_epi32(
                fours,
                _mm_and_si128(t, ones)));
            sum_abs_1 = _mm_add_epi32(sum_abs_1, min_t);
        }
        case 1: {
            __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&absLevels[nbOut->outPos[0] * 4])));
            sum_abs = _mm_add_epi32(sum_abs, t);
            sum_num = _mm_add_epi32(sum_num, _mm_min_epi32(t, ones));
            __m128i min_t = _mm_min_epi32(
              t,
              _mm_add_epi32(
                fours,
                _mm_and_si128(t, ones)));
            sum_abs_1 = _mm_add_epi32(sum_abs_1, min_t);
        }
            break;
        default:
          assert(0);
        }
        sum_abs_1 = _mm_slli_epi32(sum_abs_1, 3);
        sum_abs = _mm_slli_epi32(_mm_min_epi32(_mm_set1_epi32(127), sum_abs), 8);
        __m128i template_ctx_init = _mm_add_epi32(sum_num, sum_abs);
        template_ctx_init = _mm_add_epi32(template_ctx_init, sum_abs_1);
        __m128i shuffle_mask = _mm_setr_epi8(0, 1, 4, 5, 8, 9, 12, 13, 0, 0, 0, 0, 0, 0, 0, 0);
        __m128i shuffled_template_ctx_init = _mm_shuffle_epi8(template_ctx_init, shuffle_mask);
        temp[id % 4] = _mm_extract_epi64(shuffled_template_ctx_init, 0);
        if (id % 4 == 3) {
          all[id / 4] = _mm256_loadu_si256((__m256i const*)temp);
          last = template_ctx_init;
        }
      }
      
      _mm256_storeu_si256((__m256i*)(&state->m_ctxInit[state_offset >> 2][0]), all[0]);
      _mm256_storeu_si256((__m256i*)(&state->m_ctxInit[state_offset >> 2][16]), all[1]);
      _mm256_storeu_si256((__m256i*)(&state->m_ctxInit[state_offset >> 2][32]), all[2]);
      _mm256_storeu_si256((__m256i*)(&state->m_ctxInit[state_offset >> 2][48]), all[3]);
      
      memset(state->m_absLevels[state_offset >> 2], 0, 16 * 4);      
    }
    // End update common context

    __m128i sum_num = _mm_and_si128(last, _mm_set1_epi32(7));
    __m128i sum_abs1 = _mm_and_si128(
      _mm_srli_epi32(last, 3),
      _mm_set1_epi32(31));

    __m128i sum_abs_min = _mm_min_epi32(
      _mm_set1_epi32(3),
      _mm_srli_epi32(
        _mm_add_epi32(sum_abs1, _mm_set1_epi32(1)),
        1));

    __m128i offsets = _mm_set_epi32(12 * 3, 12 * 2, 12 * 1, 12 * 0);
    offsets = _mm_add_epi32(offsets, _mm_set1_epi32(sigCtxOffsetNext));
    offsets         = _mm_add_epi32(offsets, sum_abs_min);
    __m256i sig_frac_bits = _mm256_i32gather_epi64((long long const*)&state->m_sigFracBitsArray[state_offset][0][0], offsets, 8);
    _mm256_store_si256((__m256i*)&state->m_sigFracBits[state_offset][0], sig_frac_bits);


    __m128i sum_gt1 = _mm_sub_epi32(sum_abs1, sum_num);
    __m128i min_gt1 = _mm_min_epi32(sum_gt1, _mm_set1_epi32(4));
    uint32_t sum_gt1_s[4];
    _mm_storeu_si128((__m128i*)sum_gt1_s, min_gt1);
    // These are 192 bits so no benefit from using avx2
    for (int i = 0; i < 4; ++i) {
      memcpy(state->m_coeffFracBits[state_offset + i], state->m_gtxFracBitsArray[sum_gt1_s[i] + gtxCtxOffsetNext], sizeof(state->m_coeffFracBits[0]));
    }
  }
  else {
    for (int i = 0; i < 4; i++) {
      uvg_dep_quant_update_state_eos(
        ctxs,
        scan_pos,
        cg_pos,
        sigCtxOffsetNext,
        gtxCtxOffsetNext,
        width_in_sbb,
        height_in_sbb,
        next_sbb_right,
        next_sbb_below,
        decisions,
        i);
    }
  }
}

static INLINE void update_states_avx2(
  context_store*  ctxs,
  int             numIPos,
  const uint32_t  scan_pos,
  const Decision* decisions,
  const uint32_t  sigCtxOffsetNext,
  const uint32_t  gtxCtxOffsetNext,
  const NbInfoSbb next_nb_info_ssb,
  const int       baseLevel,
  const bool      extRiceFlag)
{
  all_depquant_states* state = &ctxs->m_allStates;

  bool all_non_negative = true;
  bool all_above_minus_two = true;
  bool all_minus_one = true;
  for (int i = 0; i < 4; ++i) {
    all_non_negative &= decisions->prevId[i] >= 0;
    all_above_minus_two &= decisions->prevId[i] > -2;
    all_minus_one &= decisions->prevId[i] == -1;
  }
  int state_offset = ctxs->m_curr_state_offset;
  __m256i rd_cost = _mm256_load_si256((__m256i const*)decisions->rdCost);
  _mm256_store_si256((__m256i *)& ctxs->m_allStates.m_rdCost[state_offset], rd_cost);
  if (all_above_minus_two) {

    bool    rem_reg_all_gte_4 = true;
    bool    rem_reg_all_lt4 = true;
    __m128i control = _mm_setr_epi8(0, 4, 8, 12, 0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1);

    __m128i abs_level = _mm_load_si128((__m128i const*)decisions->absLevel);
    if (all_non_negative) {
      __m128i prv_states_o  = _mm_load_si128((__m128i const*)decisions->prevId);
      __m128i prev_offset = _mm_set1_epi32(ctxs->m_prev_state_offset);
      __m128i prv_states     = _mm_add_epi32(prv_states_o, prev_offset);
      __m128i shuffled_prev_states = _mm_shuffle_epi8(prv_states, control);

      // sig_sbb values matter only whether they are one or zero so make sure that they stay at one or zero
      // which allows some optimizations when handling the values in update_state_eos_avx2
      __m128i sig_sbb   = _mm_load_si128((__m128i const*)state->m_numSigSbb);
      sig_sbb = _mm_shuffle_epi8(sig_sbb, shuffled_prev_states);
      __m128i has_coeff = _mm_min_epi32(abs_level, _mm_set1_epi32(1));
      has_coeff         = _mm_shuffle_epi8(has_coeff, control);
      sig_sbb           = _mm_or_si128(sig_sbb, has_coeff);
      int sig_sbb_i = _mm_extract_epi32(sig_sbb, 0);
      memcpy(&state->m_numSigSbb[state_offset], &sig_sbb_i, 4);

      // These following two are jus shuffled and then extracted the 4 bytes that store the values
      __m128i ref_sbb_ctx_idx = _mm_load_si128((__m128i const*)state->m_refSbbCtxId);
      ref_sbb_ctx_idx = _mm_shuffle_epi8(ref_sbb_ctx_idx, shuffled_prev_states);
      int ref_sbb_ctx = _mm_extract_epi32(ref_sbb_ctx_idx, 0);
      memcpy(&state->m_refSbbCtxId[state_offset], &ref_sbb_ctx, 4);
      
      __m128i go_rice_par = _mm_load_si128((__m128i const*)state->m_goRicePar);
      go_rice_par = _mm_shuffle_epi8(go_rice_par, shuffled_prev_states);
      int go_rice_par_i = _mm_extract_epi32(go_rice_par, 0);
      memcpy(&state->m_goRicePar[state_offset], &go_rice_par_i, 4);

      // Again gather is not necessary but it is easier to read and shouldn't have too large of a performance hit
      // Should be true for all gathers here
      __m256i sbb_frac_bits = _mm256_i32gather_epi64((const long long *)state->m_sbbFracBits[0], prv_states, 8);
      _mm256_store_si256((__m256i*)&state->m_sbbFracBits[state_offset][0], sbb_frac_bits);

      // Next three lines: state->m_remRegBins = prvState->m_remRegBins - 1;
      __m128i rem_reg_bins = _mm_i32gather_epi32(state->m_remRegBins, prv_states, 4);
      __m128i ones = _mm_set1_epi32(1);
      rem_reg_bins = _mm_sub_epi32(rem_reg_bins, ones);

      __m128i reg_bins_sub = _mm_set1_epi32(0);
      // Next two lines: (decision->absLevel < 2 ? (unsigned)decision->absLevel : 3)
      __m128i abs_level_smaller_than_two = _mm_cmplt_epi32(abs_level, _mm_set1_epi32(2));
      __m128i secondary = _mm_blendv_epi8(_mm_set1_epi32(3), abs_level, abs_level_smaller_than_two);

      // Depending on whether the rem_reg_bins are smaller than four or not,
      // the reg_bins_sub is either 0 or result of the above operation
      __m128i rem_reg_bins_smaller_than_four = _mm_cmplt_epi32(rem_reg_bins, _mm_set1_epi32(4));
      reg_bins_sub = _mm_blendv_epi8(secondary, reg_bins_sub, rem_reg_bins_smaller_than_four);
      rem_reg_bins = _mm_sub_epi32(rem_reg_bins, reg_bins_sub);
      _mm_store_si128((__m128i*)&state->m_remRegBins[state_offset], rem_reg_bins);

      // Save whether all rem_reg_bins are smaller than four or not and gte 4 as these
      // are needed in multiple places
      __m128i mask = _mm_cmpgt_epi32(rem_reg_bins, _mm_set1_epi32(3)); 
      int     bit_mask = _mm_movemask_epi8(mask);           
      rem_reg_all_gte_4 = (bit_mask == 0xFFFF);
      mask = _mm_cmplt_epi32(rem_reg_bins, _mm_set1_epi32(4));
      bit_mask = _mm_movemask_epi8(mask); 
      rem_reg_all_lt4 = (bit_mask == 0xFFFF);

      // This is the same as in update_state_eos_avx2
      __m128i temp_prev_state = _mm_shuffle_epi8(prv_states_o, control);
      __m256i prev_state_256  = _mm256_castsi128_si256(temp_prev_state);
      prev_state_256          = _mm256_permute4x64_epi64(prev_state_256, 0);
      __m256i temp_add        = _mm256_setr_epi32(
        0,
        0x04040404,
        0x08080808,
        0x0c0c0c0c,
        0,
        0x04040404,
        0x08080808,
        0x0c0c0c0c);
      prev_state_256 = _mm256_add_epi8(prev_state_256, temp_add);
      for (int i = 0; i < 64; i += (256 / (8 * sizeof(uint8_t)))) {
        __m256i data = _mm256_load_si256((__m256i*)&state->m_absLevels[ctxs->m_prev_state_offset >> 2][i]);
        data = _mm256_shuffle_epi8(data, prev_state_256);
        _mm256_store_si256((__m256i*)&state->m_absLevels[ctxs->m_curr_state_offset >> 2][i], data);
      }

      // This is overall the same as absLevels but since the ctx values are two bytes all of the
      // masks have to account for that
      __m256i prev_state_full = _mm256_load_si256((__m256i const*)decisions->prevId);
      __m256i shuffle_mask = _mm256_setr_epi8(0, 0, 4, 4,8, 8, 12, 12, 0, 0, 4, 4, 8, 8, 12, 12,0, 0, 0, 0,0, 0, 0, 0,16, 16, 16, 16, 16, 16, 16, 16);
      prev_state_full = _mm256_shuffle_epi8(prev_state_full, shuffle_mask);
      prev_state_full  = _mm256_permute4x64_epi64(prev_state_full, 0);
      prev_state_full  = _mm256_slli_epi16(prev_state_full, 1);
      temp_add = _mm256_setr_epi8(
        0, 1, 0, 1, 0, 1, 0, 1,
        8, 9, 8, 9, 8, 9, 8, 9, 
        16, 17, 16, 17, 16, 17, 16, 17,
        24, 25, 24, 25, 24, 25, 24, 25);
      prev_state_full = _mm256_add_epi8(prev_state_full, temp_add);

      for (int i = 0; i < 64; i += (256 / (8 * sizeof(uint16_t)))) {
         __m256i data   = _mm256_load_si256((__m256i*)(&state->m_ctxInit[(ctxs->m_prev_state_offset >> 2)][i]));
        data = _mm256_shuffle_epi8(data, prev_state_full);
        _mm256_store_si256((__m256i*)(&state->m_ctxInit[(state_offset >> 2)][i]), data);
      }
    }
    else if (all_minus_one) {
      memset(&state->m_numSigSbb[state_offset], 1, 4);
      memset(&state->m_refSbbCtxId[state_offset], -1, 4);

      const int a = (state->effWidth * state->effHeight * 28) / 16;

      __m128i   rem_reg_bins = _mm_set1_epi32(a);
      __m128i   sub = _mm_blendv_epi8(
        _mm_set1_epi32(3),
        abs_level,
        _mm_cmplt_epi32(abs_level, _mm_set1_epi32(2))
      );
      rem_reg_bins = _mm_sub_epi32(rem_reg_bins, sub);
      _mm_store_si128((__m128i*) & state->m_remRegBins[state_offset], rem_reg_bins);

      __m128i mask = _mm_cmpgt_epi32(rem_reg_bins, _mm_set1_epi32(3));
      int     bit_mask = _mm_movemask_epi8(mask);
      rem_reg_all_gte_4 = (bit_mask == 0xFFFF);
      mask = _mm_cmplt_epi32(rem_reg_bins, _mm_set1_epi32(4));
      bit_mask = _mm_movemask_epi8(mask);
      rem_reg_all_lt4 = (bit_mask == 0xFFFF);
      
      memset(state->m_absLevels[state_offset >> 2], 0, 16 * sizeof(uint8_t) * 4);
      memset(state->m_ctxInit[state_offset >> 2], 0, 16 * sizeof(uint16_t) * 4);
      
    }
    else {
      for (int i = 0; i< 4; ++i) {
        const int decision_id = i;
        const int state_id = state_offset + i;
        if (decisions->prevId[decision_id] >= 0) {
          const int prvState = ctxs->m_prev_state_offset + decisions->prevId[decision_id];
          state->m_numSigSbb[state_id] = (state->m_numSigSbb[prvState]) || !!decisions->absLevel[decision_id];
          state->m_refSbbCtxId[state_id] = state->m_refSbbCtxId[prvState];
          state->m_sbbFracBits[state_id][0] = state->m_sbbFracBits[prvState][0];
          state->m_sbbFracBits[state_id][1] = state->m_sbbFracBits[prvState][1];
          state->m_remRegBins[state_id] = state->m_remRegBins[prvState] - 1;
          state->m_goRicePar[state_id] = state->m_goRicePar[prvState];
          if (state->m_remRegBins[state_id] >= 4) {
            state->m_remRegBins[state_id] -= (decisions->absLevel[decision_id] < 2 ? (unsigned)decisions->absLevel[decision_id] : 3);
          }
        } else {
          state->m_numSigSbb[state_id] = 1;
          state->m_refSbbCtxId[state_id] = -1;
          int ctxBinSampleRatio = 28;
          state->m_remRegBins[state_id] = (state->effWidth * state->effHeight * ctxBinSampleRatio) / 16 - (decisions->absLevel[decision_id] < 2 ? (unsigned)decisions->absLevel[decision_id] : 3);
        }
        rem_reg_all_gte_4 &= state->m_remRegBins[state_id] >= 4;
        rem_reg_all_lt4 &= state->m_remRegBins[state_id] < 4;
      }
      {
        // Same as for the all_non_negative but use blendv to set the shuffle mask to -1 for the states that do not have previous state
        __m256i prev_state_full = _mm256_load_si256((__m256i const*)decisions->prevId);
        __m256i shuffle_mask = _mm256_setr_epi8(0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
        prev_state_full = _mm256_shuffle_epi8(prev_state_full, shuffle_mask);
        prev_state_full  = _mm256_permute4x64_epi64(prev_state_full, 0);
        __m256i temp_add        = _mm256_setr_epi32(
          0,
          0x04040404,
          0x08080808,
          0x0c0c0c0c,
          0,
          0x04040404,
          0x08080808,
          0x0c0c0c0c);
        __m256i comp_mask = _mm256_cmpeq_epi8(prev_state_full, _mm256_set1_epi64x(-1));
        prev_state_full = _mm256_add_epi8(prev_state_full, temp_add);
        prev_state_full = _mm256_blendv_epi8(prev_state_full, _mm256_set1_epi64x(-1), comp_mask);
        for (int i = 0; i < 64; i += (256 / (8 * sizeof(uint8_t)))) {
          __m256i data = _mm256_load_si256((__m256i*)&state->m_absLevels[ctxs->m_prev_state_offset >> 2][i]);
          data = _mm256_shuffle_epi8(data, prev_state_full);
          _mm256_store_si256((__m256i*)&state->m_absLevels[ctxs->m_curr_state_offset >> 2][i], data);
        }
      }

      {
        __m256i prev_state_full = _mm256_load_si256((__m256i const*)decisions->prevId);
        __m256i shuffle_mask = _mm256_setr_epi8(0, 0, 4, 4,8, 8, 12, 12, 0, 0, 4, 4, 8, 8, 12, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
        prev_state_full = _mm256_shuffle_epi8(prev_state_full, shuffle_mask);
        prev_state_full  = _mm256_permute4x64_epi64(prev_state_full, 0);
        __m256i comp_mask = _mm256_cmpeq_epi8(prev_state_full, _mm256_set1_epi64x(-1));
        prev_state_full  = _mm256_slli_epi16(prev_state_full, 1);
        __m256i temp_add = _mm256_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 8, 9, 8, 9, 8, 9, 8, 9, 16, 17, 16, 17,16, 17,16, 17, 24, 25,24,25,24,25,24,25);

        prev_state_full = _mm256_add_epi8(prev_state_full, temp_add);
        prev_state_full = _mm256_blendv_epi8(prev_state_full, _mm256_set1_epi64x(-1), comp_mask);

        for (int i = 0; i < 64; i += (256 / 8 / sizeof(uint16_t))) {
          __m256i data   = _mm256_load_si256((__m256i*)(&state->m_ctxInit[(ctxs->m_prev_state_offset >> 2)][i]));
          data = _mm256_shuffle_epi8(data, prev_state_full);
          _mm256_store_si256((__m256i*)(&state->m_ctxInit[(state_offset >> 2)][i]), data);
        }
      }
    }
    uint32_t level_offset   = scan_pos & 15;
    __m128i  max_abs        = _mm_min_epi32(abs_level, _mm_set1_epi32(255));
    max_abs                 = _mm_shuffle_epi8(max_abs, control);
    uint32_t packed_max_abs = _mm_extract_epi32(max_abs, 0);
    memcpy(&state->m_absLevels[state_offset >> 2][level_offset * 4], &packed_max_abs,4);

    state->all_gte_four = rem_reg_all_gte_4;
    state->all_lt_four = rem_reg_all_lt4;

    if (rem_reg_all_gte_4) {
      const __m128i  ones = _mm_set1_epi32(1);
      const uint32_t tinit_offset = MIN(level_offset - 1u, 15u);
      __m128i   tinit = _mm_loadu_si128((__m128i*)(&state->m_ctxInit[state_offset >> 2][tinit_offset * 4]));
      tinit = _mm_cvtepu16_epi32(tinit); 
      __m128i sum_abs1 = _mm_and_si128(_mm_srli_epi32(tinit, 3), _mm_set1_epi32(31));
      __m128i sum_num = _mm_and_si128(tinit, _mm_set1_epi32(7));

      uint8_t* levels = (uint8_t*)state->m_absLevels[state_offset >> 2];
      switch (numIPos) {
      case 5:
        {
          __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&levels[next_nb_info_ssb.inPos[4] * 4])));
          __m128i min_arg = _mm_min_epi32(
            _mm_add_epi32(_mm_set1_epi32(4), _mm_and_si128(t, ones)),
            t
          );
          sum_abs1 = _mm_add_epi32(
            sum_abs1,
            min_arg
          );
          sum_num = _mm_add_epi32(
            sum_num,
            _mm_min_epi32(t, ones));
        }
      case 4:
        {
          __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&levels[next_nb_info_ssb.inPos[3] * 4])));
          __m128i min_arg = _mm_min_epi32(
            _mm_add_epi32(_mm_set1_epi32(4), _mm_and_si128(t, ones)),
            t
          );
          sum_abs1 = _mm_add_epi32(
            sum_abs1,
            min_arg
          );
          sum_num = _mm_add_epi32(sum_num, _mm_min_epi32(t, ones));
        }
      case 3:
        {
          __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&levels[next_nb_info_ssb.inPos[2] * 4])));
          __m128i min_arg = _mm_min_epi32(
            _mm_add_epi32(_mm_set1_epi32(4), _mm_and_si128(t, ones)),
            t
          );
          sum_abs1 = _mm_add_epi32(
            sum_abs1,
            min_arg
          );
          sum_num = _mm_add_epi32(sum_num, _mm_min_epi32(t, ones));
        }
      case 2:
        {
          __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&levels[next_nb_info_ssb.inPos[1] * 4])));
        __m128i min_arg = _mm_min_epi32(
              _mm_add_epi32(_mm_set1_epi32(4), _mm_and_si128(t, ones)),
              t
            );
          sum_abs1 = _mm_add_epi32(
            sum_abs1,
            min_arg
          );
          sum_num = _mm_add_epi32(sum_num, _mm_min_epi32(t, ones));
        }
      case 1: {
          __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&levels[next_nb_info_ssb.inPos[0] * 4])));
          __m128i min_arg = _mm_min_epi32(
            _mm_add_epi32(_mm_set1_epi32(4), _mm_and_si128(t, ones)),
            t
          );
          sum_abs1 = _mm_add_epi32(
            sum_abs1,
            min_arg
            );
          sum_num = _mm_add_epi32(sum_num, _mm_min_epi32(t, ones));
        } break;
      default:
          assert(0);
      }
      __m128i sum_gt1 = _mm_sub_epi32(sum_abs1, sum_num);
      __m128i  offsets = _mm_set_epi32(12 * 3, 12 * 2, 12 * 1, 12 * 0);
      offsets = _mm_add_epi32(offsets, _mm_set1_epi32(sigCtxOffsetNext));
      __m128i temp = _mm_min_epi32(
        _mm_srli_epi32(_mm_add_epi32(sum_abs1, ones), 1),
        _mm_set1_epi32(3));
      offsets = _mm_add_epi32(offsets, temp);
      __m256i sig_frac_bits = _mm256_i32gather_epi64((const long long *)state->m_sigFracBitsArray[state_offset][0], offsets, 8);
      _mm256_store_si256((__m256i*)&state->m_sigFracBits[state_offset][0], sig_frac_bits);

      sum_gt1 = _mm_min_epi32(sum_gt1, _mm_set1_epi32(4));
      sum_gt1 = _mm_add_epi32(sum_gt1, _mm_set1_epi32(gtxCtxOffsetNext));
      uint32_t sum_gt1_s[4];
      _mm_storeu_si128((__m128i*)sum_gt1_s, sum_gt1);
      for (int i = 0; i < 4; ++i) {
        memcpy(state->m_coeffFracBits[state_offset + i], state->m_gtxFracBitsArray[sum_gt1_s[i]], sizeof(state->m_coeffFracBits[0]));
      }

      __m128i sum_abs = _mm_srli_epi32(tinit, 8);
      sum_abs = _mm_min_epi32(sum_abs, _mm_set1_epi32(255));
      switch (numIPos) {
        case 5:
          {
          __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&levels[next_nb_info_ssb.inPos[4] * 4])));
          sum_abs = _mm_add_epi32(t, sum_abs);
          }
        case 4:
          {
          __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&levels[next_nb_info_ssb.inPos[3] * 4])));
          sum_abs = _mm_add_epi32(t, sum_abs);
          }
        case 3:
          {
          __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&levels[next_nb_info_ssb.inPos[2] * 4])));
          sum_abs = _mm_add_epi32(t, sum_abs);
          }
        case 2:
          {
          __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&levels[next_nb_info_ssb.inPos[1] * 4])));
          sum_abs = _mm_add_epi32(t, sum_abs);
          }
        case 1:
          {
          __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&levels[next_nb_info_ssb.inPos[0] * 4])));
          sum_abs = _mm_add_epi32(t, sum_abs);
          } break;
        default:
          assert(0);
      }
      if (extRiceFlag) {
        assert(0 && "Not implemented for avx2");
      } else {
        // int sumAll = MAX(MIN(31, (int)sumAbs - 4 * 5), 0);
        __m128i sum_all = _mm_max_epi32(
          _mm_min_epi32(
            _mm_set1_epi32(31),
            _mm_sub_epi32(sum_abs, _mm_set1_epi32(20))),
          _mm_set1_epi32(0));
        __m128i temp = _mm_i32gather_epi32(g_goRiceParsCoeff, sum_all, 4);
        __m128i control = _mm_setr_epi8(0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
        __m128i go_rice_par = _mm_shuffle_epi8(temp, control);
        int     go_rice_par_i = _mm_extract_epi32(go_rice_par, 0);
        memcpy(&state->m_goRicePar[state_offset], &go_rice_par_i, 4);
      }
    }

    else if (rem_reg_all_lt4) {
      uint8_t*       levels = (uint8_t*)state->m_absLevels[state_offset >> 2];
      const uint32_t tinit_offset = MIN(level_offset - 1u, 15u);
      __m128i   tinit = _mm_loadu_si128((__m128i*)(&state->m_ctxInit[state_offset >> 2][tinit_offset * 4]));
      tinit = _mm_cvtepu16_epi32(tinit); 
      __m128i sum_abs = _mm_srli_epi32(tinit, 8);
      sum_abs         = _mm_min_epi32(sum_abs, _mm_set1_epi32(255));
      switch (numIPos) {
        case 5:
          {
          __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&levels[next_nb_info_ssb.inPos[4] * 4])));
          sum_abs = _mm_add_epi32(t, sum_abs);
          }
        case 4:
          {
          __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&levels[next_nb_info_ssb.inPos[3] * 4])));
          sum_abs = _mm_add_epi32(t, sum_abs);
          }
        case 3:
          {
          __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&levels[next_nb_info_ssb.inPos[2] * 4])));
          sum_abs = _mm_add_epi32(t, sum_abs);
          }
        case 2:
          {
          __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&levels[next_nb_info_ssb.inPos[1] * 4])));
          sum_abs = _mm_add_epi32(t, sum_abs);
          }
        case 1:
          {
          __m128i t = _mm_cvtepu8_epi32(_mm_loadu_si128((__m128i*)(&levels[next_nb_info_ssb.inPos[0] * 4])));
          sum_abs = _mm_add_epi32(t, sum_abs);
          } break;
        default:
          assert(0);
      }
      if (extRiceFlag) {
        assert(0 && "Not implemented for avx2");
      } else {
        __m128i sum_all = _mm_min_epi32(_mm_set1_epi32(31), sum_abs);
        __m128i temp = _mm_i32gather_epi32(g_goRiceParsCoeff, sum_all, 4);
        __m128i control = _mm_setr_epi8(0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
        __m128i go_rice_par = _mm_shuffle_epi8(temp, control);
        int     go_rice_par_i = _mm_extract_epi32(go_rice_par, 0);
        memcpy(&state->m_goRicePar[state_offset], &go_rice_par_i, 4);

        // This cannot be vectorized because there is no way to dynamically shift values
        for (int i = 0; i < 4; ++i) {
          state->m_goRiceZero[state_offset + i] = (i < 2 ? 1 : 2) << state->m_goRicePar[state_offset + i];          
        }
      }

    }
    else {
      for (int i = 0; i < 4; ++i) {
        const int state_id = state_offset + i;
        uint8_t*  levels = (uint8_t*)(state->m_absLevels[state_offset >> 2]);
        if (state->m_remRegBins[state_id] >= 4) {
          coeff_t tinit = state->m_ctxInit[state_offset >> 2][((scan_pos - 1) & 15) * 4 + i];
          coeff_t sumAbs1 = (tinit >> 3) & 31;
          coeff_t sumNum = tinit & 7;
#define UPDATE(k)                                  \
  {                                                \
    coeff_t t = levels[next_nb_info_ssb.inPos[k] * 4 + i]; \
    sumAbs1 += MIN(4 + (t & 1), t);                \
    sumNum += !!t;                                 \
  }
          switch (numIPos) {
            case 5: UPDATE(4);
            case 4: UPDATE(3);
            case 3: UPDATE(2);
            case 2: UPDATE(1);
            case 1: UPDATE(0); break;
            default: assert(0);
          }
#undef UPDATE
          coeff_t sumGt1 = sumAbs1 - sumNum;
          state->m_sigFracBits[state_id][0] = state->m_sigFracBitsArray[state_id][sigCtxOffsetNext + MIN((sumAbs1 + 1) >> 1, 3)][0];
          state->m_sigFracBits[state_id][1] = state->m_sigFracBitsArray[state_id][sigCtxOffsetNext + MIN((sumAbs1 + 1) >> 1, 3)][1];
          memcpy(state->m_coeffFracBits[state_id], state->m_gtxFracBitsArray[gtxCtxOffsetNext + (sumGt1 < 4 ? sumGt1 : 4)], sizeof(state->m_coeffFracBits[0]));


          coeff_t sumAbs = state->m_ctxInit[state_offset >> 2][((scan_pos - 1) & 15) * 4 + i] >> 8;
#define UPDATE(k)                                  \
  {                                                \
    coeff_t t = levels[next_nb_info_ssb.inPos[k] * 4 + i]; \
    sumAbs += t;                                   \
  }
          switch (numIPos) {
            case 5: UPDATE(4);
            case 4: UPDATE(3);
            case 3: UPDATE(2);
            case 2: UPDATE(1);
            case 1: UPDATE(0); break;
            default: assert(0);
          }
#undef UPDATE
          if (extRiceFlag) {
            assert(0 && "Not implemented for avx2");
          } else {
            int sumAll = MAX(MIN(31, (int)sumAbs - 4 * 5), 0);
            state->m_goRicePar[state_id] = g_goRiceParsCoeff[sumAll];
          }
        } else {
          coeff_t sumAbs = (state->m_ctxInit[state_offset >> 2][((scan_pos - 1) & 15) * 4 + i]) >> 8;
#define UPDATE(k)                                  \
  {                                                \
    coeff_t t = levels[next_nb_info_ssb.inPos[k] * 4 + i]; \
    sumAbs += t;                                   \
  }
          switch (numIPos) {
            case 5: UPDATE(4);
            case 4: UPDATE(3);
            case 3: UPDATE(2);
            case 2: UPDATE(1);
            case 1: UPDATE(0); break;
            default: assert(0);
          }
#undef UPDATE
          if (extRiceFlag) {
            assert(0 && "Not implemented for avx2");
          } else {
            sumAbs = MIN(31, sumAbs);
            state->m_goRicePar[state_id] = g_goRiceParsCoeff[sumAbs];
          }
          state->m_goRiceZero[state_id] = ((state_id & 3) < 2 ? 1 : 2) << state->m_goRicePar[state_id];
        }
      }
    }
  } else {
    for (int i = 0; i < 4; ++i) {
      state->all_gte_four = true;
      state->all_lt_four = true;
      uvg_dep_quant_update_state(
        ctxs,
        numIPos,
        scan_pos,
        decisions,
        sigCtxOffsetNext,
        gtxCtxOffsetNext,
        next_nb_info_ssb,
        baseLevel,
        extRiceFlag,
        i);
    }
  }
}

void uvg_dep_quant_decide_and_update_avx2(
  rate_estimator_t*                         re,
  context_store*                          ctxs,
  struct dep_quant_scan_info const* const scan_info,
  const coeff_t                           absCoeff,
  const uint32_t                          scan_pos,
  const uint32_t                          width_in_sbb,
  const uint32_t                          height_in_sbb,
  const NbInfoSbb                         next_nb_info_ssb,
  bool                                    zeroOut,
  coeff_t                                 quantCoeff,
  const uint32_t                          effWidth,
  const uint32_t                          effHeight,
  bool                                    is_chroma)
{
  Decision* decisions = &ctxs->m_trellis[scan_pos];
  SWAP(ctxs->m_curr_state_offset, ctxs->m_prev_state_offset, int);

  enum ScanPosType spt = 0;
  if ((scan_pos & 15) == 15 && scan_pos > 16 && scan_pos < effHeight * effWidth - 1)
  {
    spt = SCAN_SOCSBB;
  }
  else if ((scan_pos & 15) == 0 && scan_pos > 0 && scan_pos < effHeight * effWidth - 16)
  {
    spt = SCAN_EOCSBB;
  }

  xDecide(&ctxs->m_allStates, &ctxs->m_startState, ctxs->m_quant, spt, absCoeff, re->m_lastBitsX[scan_info->pos_x] + re->m_lastBitsY[scan_info->pos_y], decisions, zeroOut, quantCoeff,ctxs->m_skip_state_offset, ctxs->m_prev_state_offset);
  decisions->zero_out = zeroOut;

  if (scan_pos) {
    if (!(scan_pos & 15)) {
      SWAP(ctxs->m_common_context.m_curr_sbb_ctx_offset, ctxs->m_common_context.m_prev_sbb_ctx_offset, int);
      update_state_eos_avx2(ctxs, scan_pos, scan_info->cg_pos, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], width_in_sbb, height_in_sbb, scan_info->next_sbb_right, scan_info->next_sbb_below, decisions);
      memcpy(decisions->prevId + 4, decisions->prevId, 4 * sizeof(int32_t));
      memcpy(decisions->absLevel + 4, decisions->absLevel, 4 * sizeof(int32_t));
      memcpy(decisions->rdCost + 4, decisions->rdCost, 4 * sizeof(int64_t));
    } else if (!zeroOut) {
      update_states_avx2(ctxs, next_nb_info_ssb.num, scan_pos, decisions, scan_info->sig_ctx_offset[is_chroma], scan_info->gtx_ctx_offset[is_chroma], next_nb_info_ssb, 4, false);
    }

    if (spt == SCAN_SOCSBB) {
      SWAP(ctxs->m_skip_state_offset, ctxs->m_prev_state_offset, int);
    }
  }
}


void uvg_find_first_non_zero_avx2(const coeff_t* srcCoeff, const bool enableScalingLists, const context_store * const dep_quant_context, const uint32_t* const scan, const int32_t* q_coeff, int* firstTestPos, const int width, const int height)
{
  const int default_quant_coeff = dep_quant_context->m_quant->m_QScale;
  const int32_t thres  = dep_quant_context->m_quant->m_thresLast;
  int temp = *firstTestPos;
  if (enableScalingLists) {
    for (; temp >= 0; (temp)--) {
      coeff_t thresTmp = thres / (4 * q_coeff[scan[(temp)]]);
      if (abs(srcCoeff[scan[(temp)]]) > thresTmp) {
        break;
      }
    }
  } else {
    coeff_t thresTmp = thres / (4 * default_quant_coeff);
    if (temp >= 16 && height >= 4) {
      __m256i th = _mm256_set1_epi16(thresTmp);
      temp -= 15;
      for (; temp >= 0; temp -= 16) {
        __m256i sbb_data;
        if (width <= 4) {
          sbb_data = _mm256_loadu_si256((__m256i const*)&srcCoeff[scan[temp]]);
        } else if (width == 8) {
          uint32_t i     = scan[temp];
          __m256i  first = _mm256_loadu_si256((__m256i const*)&srcCoeff[i]);
          __m256i  second = _mm256_loadu_si256((__m256i const*)&srcCoeff[i+ 12]);
          sbb_data       = _mm256_blend_epi32(first, second, 204);
        } else {
          int16_t temp_d[16];
          uint32_t i = scan[temp];
          memcpy(temp_d, &srcCoeff[i], 8);
          i += width;
          memcpy(temp_d + 4, &srcCoeff[i], 8);
          i += width;
          memcpy(temp_d + 8, &srcCoeff[i], 8);
          i += width;
          memcpy(temp_d + 12, &srcCoeff[i], 8);

          sbb_data = _mm256_loadu_si256((__m256i const*)temp_d);
        }
        sbb_data = _mm256_abs_epi16(sbb_data);

        __m256i a = _mm256_cmpgt_epi16(sbb_data, th);
        if (!_mm256_testz_si256(a, a))
        {
          if (temp >= 0) {
            temp += 15;
          }
          break;
        }
      }
    }
    for (;temp >= 0; temp--) {
      if (abs(srcCoeff[scan[(temp)]]) > thresTmp) {
        break;
      }
    }
  }

  *firstTestPos = temp;
}


#endif //COMPILE_INTEL_AVX2 && defined X86_64

int uvg_strategy_register_depquant_avx2(void* opaque, uint8_t bitdepth)
{
  bool success = true;

#if COMPILE_INTEL_AVX2 && defined X86_64
  success &= uvg_strategyselector_register(opaque, "dep_quant_decide_and_update", "avx2", 40, &uvg_dep_quant_decide_and_update_avx2);
  success &= uvg_strategyselector_register(opaque, "find_first_non_zero_coeff", "avx2", 40, &uvg_find_first_non_zero_avx2);
#endif //COMPILE_INTEL_AVX2 && defined X86_64

  return success;
}
