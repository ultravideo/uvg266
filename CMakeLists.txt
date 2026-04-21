cmake_minimum_required(VERSION 3.22)

include(cmake/build_options.cmake)
include(cmake/ccache.cmake)
include(cmake/compiler_options.cmake)
include(cmake/cpu.cmake)
include(cmake/git.cmake)
include(cmake/project_settings.cmake)
include(cmake/system.cmake)
include(cmake/try_compile.cmake)
include(cmake/utilities.cmake)

if(NOT DEFINED CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 11)
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

project(uvg266
  LANGUAGES C CXX
  HOMEPAGE_URL https://github.com/ultravideo/uvg266
  DESCRIPTION "An open-source VVC encoder licensed under 3-clause BSD"
  VERSION 0.8.1
)

uvg_option_enum(UVG_BIT_DEPTH "Supported bit depth" "8;10" "8")
uvg_option_bool(UVG_BUILD_EXECUTABLE "Build the uvg266 executable" ON)
uvg_option_bool(UVG_BUILD_SHARED "Build libuvg266 as a shared library" ON)
uvg_option_bool(UVG_BUILD_STATIC "Build libuvg266 as a static library" ON)
uvg_option_bool(UVG_BUILD_TESTS "Build tests" ON)
uvg_option_bool(UVG_ENABLE_ALTIVEC "Enable AltiVec optimizations (if supported)" ON)
uvg_option_bool(UVG_ENABLE_AVX2 "Enable AVX2 optimizations (if supported)" ON)
uvg_option_bool(UVG_ENABLE_SSE2 "Enable SSE2 optimizations (if supported)" ON)
uvg_option_bool(UVG_ENABLE_SSE41 "Enable SSE4.1 optimizations (if supported)" ON)
uvg_option_bool(UVG_ENABLE_SSE42 "Enable SSE4.2 optimizations (if supported)" ON)
uvg_option_bool(UVG_ENABLE_CCACHE "Use ccache (if available)" ON)
uvg_option_bool(UVG_SYNC_SUBMODULES "Ensure git submodules are up-to-date" ON)

uvg_maybe_override_project_version()
uvg_maybe_update_submodules()
uvg_setup_project_settings()
uvg_validate_build_options()
uvg_setup_compiler_options()
uvg_detect_system_features()
uvg_detect_cpu_features()
uvg_maybe_use_ccache()

add_subdirectory(src)

if(UVG_BUILD_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()
