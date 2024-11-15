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

#include "strategies/strategies-alf.h"
#include "strategies/sse41/alf-sse41.h"
#if defined(__AVX512F__)
#include "strategies/avx2/alf-avx2.h"
#endif
#include "strategies/generic/alf-generic.h"
#include "strategyselector.h"


// Define function pointers.
alf_derive_classification_blk_func* uvg_alf_derive_classification_blk;
alf_filter_5x5_blk_func* uvg_alf_filter_5x5_blk;
alf_filter_7x7_blk_func* uvg_alf_filter_7x7_blk;
alf_get_blk_stats_func* uvg_alf_get_blk_stats;

int uvg_strategy_register_alf(void* opaque, uint8_t bitdepth) {
  bool success = true;

  success &= uvg_strategy_register_alf_generic(opaque, bitdepth);

  if (uvg_g_hardware_flags.intel_flags.sse41) {
    success &= uvg_strategy_register_alf_sse41(opaque, bitdepth);
  }
#if defined(__AVX512F__)
  if (uvg_g_hardware_flags.intel_flags.avx2) {
    success &= uvg_strategy_register_alf_avx2(opaque, bitdepth);
  }
#endif
  return success;
}
