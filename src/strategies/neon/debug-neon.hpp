#ifndef STRATEGIES_DEBUG_NEON_HPP_
#define STRATEGIES_DEBUG_NEON_HPP_
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
 * Arm Neon optimizations.
 */

#include "config-neon.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>

// ---------------------------------------------------------
// NEON_ERROR_MSG
// ---------------------------------------------------------

#define NEON_ERROR_MSG_FULL(...)\
  fprintf(stderr, "%s:%i: %s: ", __FILE__, __LINE__, __func__);\
  fprintf(stderr, __VA_ARGS__)

#define NEON_ERROR_MSG_BRIEF(...)\
  fprintf(stderr, __VA_ARGS__)

#if NEON_ENABLE_DEBUG
  #define NEON_ERROR_MSG NEON_ERROR_MSG_FULL
#else
  #define NEON_ERROR_MSG NEON_ERROR_MSG_BRIEF
#endif

// ---------------------------------------------------------
// NEON_UNREACHABLE
// ---------------------------------------------------------

#if NEON_ENABLE_DEBUG
  #define NEON_UNREACHABLE(...)\
    NEON_ERROR_MSG(__VA_ARGS__);\
    abort()
#else
  #if defined(__has_builtin) && __has_builtin(__builtin_unreachable)
    #define NEON_UNREACHABLE(...)\
      __builtin_unreachable()
  #else
    #define NEON_UNREACHABLE(...)\
      NEON_ERROR_MSG(__VA_ARGS__);\
      abort()
  #endif
#endif

// ---------------------------------------------------------
// NEON_ASSERT, NEON_ALWAYS_ASSERT
// ---------------------------------------------------------

#define NEON_ACTUALLY_ASSERT(expr)\
  if (!(expr)) {\
    NEON_ERROR_MSG("Assertion failed: '%s'\n", #expr);\
    abort();\
  } else {}

#define NEON_SKIP_ASSERT(expr)\
  ((void)(expr))

#if NEON_ENABLE_DEBUG
  #define NEON_ASSERT NEON_ACTUALLY_ASSERT
#else
  #define NEON_ASSERT NEON_SKIP_ASSERT
#endif

#define NEON_ALWAYS_ASSERT NEON_ACTUALLY_ASSERT

// ---------------------------------------------------------

#endif // STRATEGIES_DEBUG_NEON_HPP_
