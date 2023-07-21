#ifndef STRATEGIES_DEPQUANT_H_
#define STRATEGIES_DEPQUANT_H_
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
 * \ingroup Optimization
 * \file
 * Interface for sao functions.
 */

#include "encoder.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "uvg266.h"
#include "dep_quant.h"


// Declare function pointers.
typedef int(dep_quant_decide_and_update_func)(
  rate_estimator_t*                       re,
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
  bool                                    is_chroma);

typedef void (find_first_non_zero_coeff_func)(
  const coeff_t*             srcCoeff,
  const bool                 enableScalingLists,
  const context_store* const dep_quant_context,
  const uint32_t* const      scan,
  const int32_t*             q_coeff,
  int*                       firstTestPos,
  int                        width,
  int                        height);


// Declare function pointers.
extern dep_quant_decide_and_update_func* uvg_dep_quant_decide_and_update;
extern find_first_non_zero_coeff_func* uvg_find_first_non_zero_coeff;

int uvg_strategy_register_depquant(void* opaque, uint8_t bitdepth);


#define STRATEGIES_DEPQUANT_EXPORTS \
  {"dep_quant_decide_and_update", (void**)&uvg_dep_quant_decide_and_update}, \
  {"find_first_non_zero_coeff", (void**)&uvg_find_first_non_zero_coeff}, \



#endif //STRATEGIES_DEPQUANT_H_
