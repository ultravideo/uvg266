############################################################
# Generate files
############################################################

string(TIMESTAMP UVG_BUILD_DATE %Y-%m-%d)
set(UVG_COMPILER_STRING "${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/build_config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/build_config.h"
  @ONLY
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/version.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/version.h"
  @ONLY
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/uvg266.pc.in"
  "${CMAKE_CURRENT_BINARY_DIR}/uvg266.pc"
  @ONLY
)

############################################################
# Target: libuvg266 -> uvg266::libuvg266_object
############################################################

add_library(libuvg266 OBJECT
  alf.c
  bitstream.c
  cabac.c
  cfg.c
  checkpoint.c
  constraint.c
  context.c
  cu.c
  debug.c
  dep_quant.c
  encode_coding_tree.c
  encoder_state-bitstream.c
  encoder_state-ctors_dtors.c
  encoder_state-geometry.c
  encoder.c
  encoderstate.c
  encoding_resume.c
  fast_coeff_cost.c
  filter.c
  hashmap.c
  image.c
  imagelist.c
  input_frame_buffer.c
  inter.c
  intra.c
  ml_classifier_intra_depth_pred.c
  ml_intra_cu_depth_pred.c
  nal.c
  rate_control.c
  rdo.c
  reshape.c
  sao.c
  scalinglist.c
  search_ibc.c
  search_inter.c
  search_intra.c
  search.c
  strategyselector.c
  tables.c
  threadqueue.c
  transform.c
  uvg266.c
  videoframe.c
  extras/libmd5.c
  strategies/strategies-alf.c
  strategies/strategies-dct.c
  strategies/strategies-depquant.c
  strategies/strategies-encode.c
  strategies/strategies-intra.c
  strategies/strategies-ipol.c
  strategies/strategies-nal.c
  strategies/strategies-picture.c
  strategies/strategies-quant.c
  strategies/strategies-sao.c
)

foreach(source
  strategies/generic/alf-generic.c
  strategies/generic/dct-generic.c
  strategies/generic/depquant-generic.c
  strategies/generic/encode_coding_tree-generic.c
  strategies/generic/intra-generic.c
  strategies/generic/ipol-generic.c
  strategies/generic/nal-generic.c
  strategies/generic/picture-generic.c
  strategies/generic/quant-generic.c
  strategies/generic/sao-generic.c
)
  target_sources(libuvg266 PRIVATE ${source})
  set_property(SOURCE ${source} APPEND_STRING PROPERTY COMPILE_FLAGS " ${UVG_CFLAGS_GENERIC} ${CFLAGS_GENERIC} $ENV{CFLAGS_GENERIC}")
endforeach()

if(UVG_HAVE_POWERPC_ALTIVEC)
  foreach(source
    strategies/altivec/picture-altivec.c
  )
    target_sources(libuvg266 PRIVATE ${source})
    set_property(SOURCE ${source} APPEND_STRING PROPERTY COMPILE_FLAGS " ${UVG_CFLAGS_ALTIVEC} ${CFLAGS_ALTIVEC} $ENV{CFLAGS_ALTIVEC}")
  endforeach()
endif()

if(UVG_HAVE_X86_64_AVX2)
  foreach(source
    strategies/avx2/alf-avx2.c
    strategies/avx2/dct-avx2.c
    strategies/avx2/depquant-avx2.c
    # FIXME: Not yet converted to VVC and fails to build.
    #strategies/avx2/encode_coding_tree-avx2.c
    strategies/avx2/intra-avx2.c
    strategies/avx2/ipol-avx2.c
    strategies/avx2/picture-avx2.c
    strategies/avx2/quant-avx2.c
    strategies/avx2/sao-avx2.c
  )
    target_sources(libuvg266 PRIVATE ${source})
    set_property(SOURCE ${source} APPEND_STRING PROPERTY COMPILE_FLAGS " ${UVG_CFLAGS_AVX2} ${CFLAGS_AVX2} $ENV{CFLAGS_AVX2}")
  endforeach()
endif()

if(UVG_HAVE_X86_64_SSE2)
  foreach(source
    strategies/sse2/picture-sse2.c
  )
    target_sources(libuvg266 PRIVATE ${source})
    set_property(SOURCE ${source} APPEND_STRING PROPERTY COMPILE_FLAGS " ${UVG_CFLAGS_SSE2} ${CFLAGS_SSE2} $ENV{CFLAGS_SSE2}")
  endforeach()
endif()

if(UVG_HAVE_X86_64_SSE41)
  foreach(source
    strategies/sse41/alf-sse41.c
    strategies/sse41/picture-sse41.c
  )
    target_sources(libuvg266 PRIVATE ${source})
    set_property(SOURCE ${source} APPEND_STRING PROPERTY COMPILE_FLAGS " ${UVG_CFLAGS_SSE41} ${CFLAGS_SSE41} $ENV{CFLAGS_SSE41}")
  endforeach()
endif()

if(UVG_HAVE_X86_64_SSE42)
  foreach(source
    strategies/sse42/picture-sse42.c
  )
    target_sources(libuvg266 PRIVATE ${source})
    set_property(SOURCE ${source} APPEND_STRING PROPERTY COMPILE_FLAGS " ${UVG_CFLAGS_SSE42} ${CFLAGS_SSE42} $ENV{CFLAGS_SSE42}")
  endforeach()
endif()

if(MSVC)
  foreach(source
    threadwrapper/src/pthread.cpp
    threadwrapper/src/semaphore.cpp
  )
    target_sources(libuvg266 PRIVATE ${source})
  endforeach()
  target_include_directories(libuvg266 SYSTEM PRIVATE threadwrapper/include)
else()
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(Threads REQUIRED)
  target_link_libraries(libuvg266 PRIVATE Threads::Threads)
  target_link_libraries(libuvg266 PRIVATE m)
endif()

if(UVG_BUILD_SHARED)
  set_target_properties(libuvg266 PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_compile_definitions(libuvg266 PUBLIC PIC)
endif()

target_include_directories(libuvg266
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}"
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_compile_definitions(libuvg266 PRIVATE UVG_DLL_EXPORTS)
set_target_properties(libuvg266 PROPERTIES EXPORT_NAME libuvg266_object)

if(UVG_BIT_DEPTH STREQUAL "8")
  target_compile_definitions(libuvg266 PUBLIC UVG_BIT_DEPTH=8)
elseif(UVG_BIT_DEPTH STREQUAL "10")
  target_compile_definitions(libuvg266 PUBLIC UVG_BIT_DEPTH=10)
else()
  message(FATAL_ERROR "Unexpected bit depth: ${UVG_BIT_DEPTH}")
endif()

############################################################
# Target: libuvg266_shared -> uvg266::libuvg266
############################################################

add_library(libuvg266_shared SHARED)
target_link_libraries(libuvg266_shared PUBLIC libuvg266)
target_compile_definitions(libuvg266_shared PUBLIC UVG_DLL_EXPORTS)

set_target_properties(libuvg266_shared PROPERTIES
  EXPORT_NAME libuvg266
)

if(WIN32)
  set_target_properties(libuvg266_shared PROPERTIES
    ARCHIVE_OUTPUT_NAME libuvg266-stub
    RUNTIME_OUTPUT_NAME libuvg266
  )
else()
  set_target_properties(libuvg266_shared PROPERTIES
    OUTPUT_NAME uvg266
  )
endif()

############################################################
# Target: libuvg266_static -> uvg266::libuvg266_static
############################################################

add_library(libuvg266_static STATIC)
target_link_libraries(libuvg266_static PUBLIC libuvg266)
target_compile_definitions(libuvg266_static PUBLIC UVG_STATIC_LIB)

set_target_properties(libuvg266_static PROPERTIES
  EXPORT_NAME libuvg266_static
)

if(WIN32)
  set_target_properties(libuvg266_static PROPERTIES
    OUTPUT_NAME libuvg266
  )
else()
  set_target_properties(libuvg266_static PROPERTIES
    OUTPUT_NAME uvg266
  )
endif()

############################################################
# Target: uvg266 -> uvg266::uvg266enc
############################################################

add_executable(uvg266
  cli.c
  encmain.c
  yuv_io.c
)

if(MSVC)
  target_sources(uvg266
    PRIVATE
      extras/getopt.c
      threadwrapper/src/pthread.cpp
      threadwrapper/src/semaphore.cpp
  )
  target_include_directories(uvg266 SYSTEM PRIVATE extras threadwrapper/include)
else()
  target_link_libraries(uvg266 PRIVATE Threads::Threads)
endif()

set_target_properties(uvg266 PROPERTIES OUTPUT_NAME uvg266)
target_include_directories(uvg266 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

if(UVG_BUILD_SHARED)
  target_link_libraries(uvg266 PRIVATE libuvg266_shared)
  set_target_properties(uvg266 PROPERTIES INSTALL_RPATH_USE_LINK_PATH ON)
  if(APPLE)
    set_target_properties(uvg266 PROPERTIES INSTALL_RPATH "@loader_path;@loader_path/../${CMAKE_INSTALL_LIBDIR}")
  elseif(UNIX)
    set_target_properties(uvg266 PROPERTIES INSTALL_RPATH "$ORIGIN:$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
  endif()
else()
  target_link_libraries(uvg266 PRIVATE libuvg266_static)
endif()

############################################################
# Install rules
############################################################

if(CMAKE_SKIP_INSTALL_RULES)
  return()
endif()

if(UVG_BUILD_EXECUTABLE)
  install(
    TARGETS uvg266
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  )
endif()

if(UVG_BUILD_SHARED)
  install(
    TARGETS libuvg266_shared
    EXPORT uvg266Targets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  )
endif()

if(UVG_BUILD_STATIC)
  install(
    TARGETS libuvg266_static
    EXPORT uvg266Targets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  )
endif()

install(
  TARGETS libuvg266
  EXPORT uvg266Targets
)

install(
  FILES uvg266.h
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/uvg266"
)

install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/uvg266.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
)

install(
  FILES ../doc/uvg266.1
  DESTINATION "${CMAKE_INSTALL_MANDIR}/man1"
)

install(
  EXPORT uvg266Targets
  FILE uvg266Targets.cmake
  NAMESPACE uvg266::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/uvg266
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
  "${CMAKE_SOURCE_DIR}/cmake/uvg266Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/uvg266Config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/uvg266
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/uvg266ConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/uvg266Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/uvg266ConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/uvg266
)
