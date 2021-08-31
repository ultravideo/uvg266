/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2021 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * Kvazaar is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#include "strategies/strategies-alf.h"
#include "strategies/sse41/alf-sse41.h"
#include "strategies/avx2/alf-avx2.h"
#include "strategies/generic/alf-generic.h"
#include "strategyselector.h"


// Define function pointers.
alf_derive_classification_blk_func* kvz_alf_derive_classification_blk;
alf_filter_5x5_blk_func* kvz_alf_filter_5x5_blk;
alf_filter_7x7_blk_func* kvz_alf_filter_7x7_blk;
alf_get_blk_stats_func* kvz_alf_get_blk_stats;

int kvz_strategy_register_alf(void* opaque, uint8_t bitdepth) {
  bool success = true;

  success &= kvz_strategy_register_alf_generic(opaque, bitdepth);

  if (kvz_g_hardware_flags.intel_flags.sse41) {
    success &= kvz_strategy_register_alf_sse41(opaque, bitdepth);
  }
  if (kvz_g_hardware_flags.intel_flags.avx2) {
    success &= kvz_strategy_register_alf_avx2(opaque, bitdepth);
  }

  return success;
}