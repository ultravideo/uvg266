#pragma once
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

/**
 * \ingroup Optimization
 * \file
 * Interface for sao functions.
 */

#include "encoder.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "alf.h"


// Declare function pointers.
typedef void (alf_derive_classification_blk_func)(encoder_state_t * const state,
  const int shift,
  const int n_height,
  const int n_width,
  const int blk_pos_x,
  const int blk_pos_y,
  const int blk_dst_x,
  const int blk_dst_y,
  const int vb_ctu_height,
  int vb_pos);

// Declare function pointers.
extern alf_derive_classification_blk_func * alf_derive_classification_blk;

int kvz_strategy_register_alf(void* opaque, uint8_t bitdepth);


#define STRATEGIES_ALF_EXPORTS \
  {"alf_derive_classification_blk", (void**) &alf_derive_classification_blk}, \
 

