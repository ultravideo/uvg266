############################################################
# Helpers
############################################################

function(uvg_try_compile file flags result)
  try_compile(
    success
    "${CMAKE_BINARY_DIR}"
    "${file}"
    COMPILE_DEFINITIONS "${flags}"
  )
  if(success)
    set(${result} 1 PARENT_SCOPE)
  else()
    set(${result} 0 PARENT_SCOPE)
  endif()
endfunction()

############################################################
# x86-64
############################################################

macro(uvg_check_have_x86_64)
  uvg_try_compile("${CMAKE_SOURCE_DIR}/cmake/checks/arch/have_x86_64.c" "" UVG_HAVE_X86_64)
endmacro()

macro(uvg_check_have_sse2)
  if(MSVC)
    set(UVG_CFLAGS_SSE2 "")
  else()
    set(UVG_CFLAGS_SSE2 "-msse2")
  endif()
  uvg_try_compile("${CMAKE_SOURCE_DIR}/cmake/checks/arch/have_x86_64_sse2.c" "${UVG_CFLAGS_SSE2}" UVG_HAVE_X86_64_SSE2)
endmacro()

macro(uvg_check_have_sse41)
  if(MSVC)
    set(UVG_CFLAGS_SSE41 "")
  else()
    set(UVG_CFLAGS_SSE41 "-msse4.1")
  endif()
  uvg_try_compile("${CMAKE_SOURCE_DIR}/cmake/checks/arch/have_x86_64_sse41.c" "${UVG_CFLAGS_SSE41}" UVG_HAVE_X86_64_SSE41)
endmacro()

macro(uvg_check_have_sse42)
  if(MSVC)
    set(UVG_CFLAGS_SSE42 "")
  else()
    set(UVG_CFLAGS_SSE42 "-msse4.2")
  endif()
  uvg_try_compile("${CMAKE_SOURCE_DIR}/cmake/checks/arch/have_x86_64_sse42.c" "${UVG_CFLAGS_SSE42}" UVG_HAVE_X86_64_SSE42)
endmacro()

macro(uvg_check_have_avx2)
  if(MSVC)
    set(UVG_CFLAGS_AVX2 "/arch:AVX2")
  else()
    set(UVG_CFLAGS_AVX2 "-mavx2 -mbmi -mpopcnt -mlzcnt -mbmi2")
  endif()
  uvg_try_compile("${CMAKE_SOURCE_DIR}/cmake/checks/arch/have_x86_64_avx2.c" "${UVG_CFLAGS_AVX2}" UVG_HAVE_X86_64_AVX2)
endmacro()

macro(uvg_maybe_enable_sse2)
  if(UVG_ENABLE_SSE2)
    message(STATUS "Detecting SSE2 support")
    uvg_check_have_sse2()
    if(UVG_HAVE_X86_64_SSE2)
      message(STATUS "Detecting SSE2 support - yes")
    else()
      message(STATUS "Detecting SSE2 support - no")
    endif()
  endif()
endmacro()

macro(uvg_maybe_enable_sse41)
  if(UVG_ENABLE_SSE41)
    message(STATUS "Detecting SSE4.1 support")
    uvg_check_have_sse41()
    if(UVG_HAVE_X86_64_SSE41)
      message(STATUS "Detecting SSE4.1 support - yes")
    else()
      message(STATUS "Detecting SSE4.1 support - no")
    endif()
  endif()
endmacro()

macro(uvg_maybe_enable_sse42)
  if(UVG_ENABLE_SSE42)
    message(STATUS "Detecting SSE4.2 support")
    uvg_check_have_sse42()
    if(UVG_HAVE_X86_64_SSE42)
      message(STATUS "Detecting SSE4.2 support - yes")
    else()
      message(STATUS "Detecting SSE4.2 support - no")
    endif()
  endif()
endmacro()

macro(uvg_maybe_enable_avx2)
  if(UVG_ENABLE_AVX2)
    message(STATUS "Detecting AVX2 support")
    uvg_check_have_avx2()
    if(UVG_HAVE_X86_64_AVX2)
      message(STATUS "Detecting AVX2 support - yes")
    else()
      message(STATUS "Detecting AVX2 support - no")
    endif()
  endif()
endmacro()

############################################################
# PowerPC
############################################################

macro(uvg_check_have_powerpc)
  uvg_try_compile("${CMAKE_SOURCE_DIR}/cmake/checks/arch/have_powerpc.c" "" UVG_HAVE_POWERPC)
endmacro()

macro(uvg_check_have_altivec)
  set(UVG_CFLAGS_ALTIVEC "")
  uvg_try_compile("${CMAKE_SOURCE_DIR}/cmake/checks/arch/have_powerpc_altivec.c" "${UVG_CFLAGS_ALTIVEC}" UVG_HAVE_POWERPC_ALTIVEC)
endmacro()

macro(uvg_maybe_enable_altivec)
  if(UVG_ENABLE_ALTIVEC)
    message(STATUS "Detecting AltiVec support")
    uvg_check_have_altivec()
    if(UVG_HAVE_POWERPC_ALTIVEC)
      message(STATUS "Detecting AltiVec support - yes")
    else()
      message(STATUS "Detecting AltiVec support - no")
    endif()
  endif()
endmacro()

############################################################
# RISC-V
############################################################

macro(uvg_check_have_riscv)
  uvg_try_compile("${CMAKE_SOURCE_DIR}/cmake/checks/arch/have_riscv.c" "" UVG_HAVE_RISCV)
endmacro()

macro(uvg_check_have_riscv_64)
  uvg_try_compile("${CMAKE_SOURCE_DIR}/cmake/checks/arch/have_riscv_64.c" "" UVG_HAVE_RISCV_64)
endmacro()

############################################################
# Generic
############################################################

macro(uvg_detect_cpu_features)
  # These variables must always be set to some value because
  # they are used in build_config.h. It would be nice to figure
  # out a more robust way of handling this:

  set(UVG_HAVE_X86_64 0)
  set(UVG_HAVE_X86_64_SSE2 0)
  set(UVG_HAVE_X86_64_SSE41 0)
  set(UVG_HAVE_X86_64_SSE42 0)
  set(UVG_HAVE_X86_64_AVX2 0)

  set(UVG_HAVE_POWERPC 0)
  set(UVG_HAVE_POWERPC_ALTIVEC 0)

  set(UVG_HAVE_RISCV 0)
  set(UVG_HAVE_RISCV_64 0)

  uvg_check_have_x86_64()
  uvg_check_have_powerpc()
  uvg_check_have_riscv()
  uvg_check_have_riscv_64()

  if(UVG_HAVE_X86_64)
    message(STATUS "Targeting x86-64")
    uvg_maybe_enable_sse2()
    uvg_maybe_enable_sse41()
    uvg_maybe_enable_sse42()
    uvg_maybe_enable_avx2()
  elseif(UVG_HAVE_POWERPC)
    message(STATUS "Targeting PowerPC")
    uvg_maybe_enable_altivec()
  elseif(UVG_HAVE_RISCV_64)
    message(STATUS "Targeting RISC-V 64-bit")
  elseif(UVG_HAVE_RISCV)
    message(STATUS "Targeting RISC-V")
  endif()
endmacro()
