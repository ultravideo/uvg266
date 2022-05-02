cmake_minimum_required(VERSION 3.12)

project(uvg266
LANGUAGES C CXX
HOMEPAGE_URL https://github.com/ultravideo/uvg266
DESCRIPTION "An open-source VVC encoder licensed under 3-clause BSD"
VERSION 0.2.3 )

# Apply dynamic info to the config files
configure_file("${PROJECT_SOURCE_DIR}/src/uvg266.pc.in" "${PROJECT_SOURCE_DIR}/src/uvg266.pc" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/src/version.h.in" "${PROJECT_SOURCE_DIR}/src/version.h" @ONLY)

option(BUILD_SHARED_LIBS "Build using shared uvg266 library" ON)

option(BUILD_TESTS "Build tests" ON)


find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(WARNING "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/greatest/greatest.h")
    message(WARNING "The submodule greatest was not loaded, some tests may fail")
endif()

# Add all sources in src/ base
file(GLOB LIB_SOURCES RELATIVE ${PROJECT_SOURCE_DIR} "src/*.h" "src/*.c")

# We don't want CLI main in the library
list(REMOVE_ITEM LIB_SOURCES "src/encmain.c" "src/cli.c" "src/cli.h" "src/yuv_io.c" "src/yuv_io.h")

# Add also all the strategies
file(GLOB_RECURSE LIB_SOURCES_STRATEGIES RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/*.c")

# ToDo: do something with encode_coding_tree-avx2, currently not converted to VVC
list(REMOVE_ITEM LIB_SOURCES_STRATEGIES "src/strategies/avx2/encode_coding_tree-avx2.c")

list(APPEND LIB_SOURCES ${LIB_SOURCES_STRATEGIES})

# We also need the libmd5
list(APPEND LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/extras/libmd5.c)

add_definitions(-DUVG_DLL_EXPORTS)

if(BUILD_SHARED_LIBS)
  add_definitions(-DPIC)
endif()

# For visual studio / windows we also need our own pthread implementation and getopt
if(MSVC)
  list(APPEND LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/extras/getopt.c ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/pthread.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/semaphore.cpp)
  add_definitions(-DWIN32_LEAN_AND_MEAN -D_WIN32 -DWIN32 -DWIN64)
endif()

if(BUILD_SHARED_LIBS)
  list( APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" "./" "../lib" )  
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  add_library(uvg266 SHARED ${LIB_SOURCES})
else()
  add_library(uvg266 STATIC ${LIB_SOURCES})
  if(MSVC) # Fix a linking problem with visual studio when the library is the same name as the binary
    set_target_properties(uvg266 PROPERTIES OUTPUT_NAME libuvg266)
  endif()
  
endif()

target_include_directories(uvg266 PUBLIC src)
target_include_directories(uvg266 PUBLIC src/extras)
target_include_directories(uvg266 PUBLIC src/strategies)

file(GLOB LIB_SOURCES_STRATEGIES_AVX2 RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/avx2/*.c")
file(GLOB LIB_SOURCES_STRATEGIES_SSE41 RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/sse41/*.c")

set(CLI_SOURCES "src/encmain.c" "src/cli.c" "src/cli.h" "src/yuv_io.c" "src/yuv_io.h")

# Add the getopt and pthread for visual studio
if(MSVC) 
  list(APPEND CLI_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/extras/getopt.c ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/pthread.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/semaphore.cpp)
endif()

add_executable(uvg266-bin ${CLI_SOURCES})

target_link_libraries(uvg266-bin PUBLIC uvg266)

set_target_properties(uvg266-bin PROPERTIES OUTPUT_NAME uvg266)
set_target_properties(uvg266-bin PROPERTIES RUNTIME_OUTPUT_NAME uvg266)


if(MSVC)
  target_include_directories(uvg266 PUBLIC src/threadwrapper/include)
  set_property( SOURCE ${LIB_SOURCES_STRATEGIES_AVX2} APPEND PROPERTY COMPILE_FLAGS "/arch:AVX2" )
else()
  list(APPEND ALLOW_AVX2 "x86_64" "AMD64")
  if(${CMAKE_SYSTEM_PROCESSOR} IN_LIST ALLOW_AVX2) 
    set_property( SOURCE ${LIB_SOURCES_STRATEGIES_AVX2} APPEND PROPERTY COMPILE_FLAGS "-mavx2 -mbmi -mpopcnt -mlzcnt -mbmi2" )
    set_property( SOURCE ${LIB_SOURCES_STRATEGIES_SSE41} APPEND PROPERTY COMPILE_FLAGS "-msse4.1" )    
  endif()
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(Threads REQUIRED)
  target_link_libraries(uvg266 PUBLIC Threads::Threads)

  include(CheckLibraryExists)

  CHECK_LIBRARY_EXISTS(m sin "" HAVE_LIB_M)

  if (HAVE_LIB_M)
      set(EXTRA_LIBS ${EXTRA_LIBS} m)
  endif (HAVE_LIB_M)

  target_link_libraries(uvg266-bin PUBLIC ${EXTRA_LIBS})
endif()

# Source grouping

# Some basic structuring of the files based on previous visual studio project files
file(GLOB SOURCE_GROUP_BITSTREAM RELATIVE ${PROJECT_SOURCE_DIR} "src/encode_coding_tree.*" "src/encoder_state-bitstream.*" "src/nal.*")
file(GLOB SOURCE_GROUP_CABAC RELATIVE ${PROJECT_SOURCE_DIR} "src/bitstream.*" "src/cabac.*" "src/context.*")
file(GLOB SOURCE_GROUP_COMPRESSION RELATIVE ${PROJECT_SOURCE_DIR} "src/search*" "src/rdo.*" "src/fast_coeff*")
file(GLOB SOURCE_GROUP_CONSTRAINT RELATIVE ${PROJECT_SOURCE_DIR} "src/constraint.*" "src/ml_*")
file(GLOB SOURCE_GROUP_CONTROL RELATIVE ${PROJECT_SOURCE_DIR} "src/cfg.*" "src/encoder.*" "src/encoder_state-c*" "src/encoder_state-g*" "src/encoderstate*" "src/gop.*" "src/input_frame_buffer.*" "src/uvg266*" "src/rate_control.*" "src/mip_data.h")
file(GLOB SOURCE_GROUP_DATA_STRUCTURES RELATIVE ${PROJECT_SOURCE_DIR} "src/cu.*" "src/image.*" "src/imagelist.*" "src/videoframe.*")
file(GLOB SOURCE_GROUP_EXTRAS RELATIVE ${PROJECT_SOURCE_DIR} "src/extras/*.h" "src/extras/*.c")
file(GLOB_RECURSE SOURCE_GROUP_STRATEGIES RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/*.h" "src/strategies/*.c")
file(GLOB SOURCE_GROUP_RECON RELATIVE ${PROJECT_SOURCE_DIR} "src/alf.*" "src/filter.*" "src/inter.*" "src/intra.*" "src/reshape.*" "src/sao.*" "src/scalinglist.*" "src/tables.*" "src/transform.*")
file(GLOB SOURCE_GROUP_THREADING RELATIVE ${PROJECT_SOURCE_DIR} "src/threadqueue.*" "src/threads.*")
file(GLOB_RECURSE SOURCE_GROUP_THREADWRAPPER RELATIVE ${PROJECT_SOURCE_DIR} "src/threadwrapper/*.cpp" "src/threadwrapper/*.h")
file(GLOB SOURCE_GROUP_TOPLEVEL RELATIVE ${PROJECT_SOURCE_DIR} "src/debug.*" "src/global.h" "src/version.h" "src/uvg_math.h" "src/checkpoint.*")

source_group( "Bitstream" FILES  ${SOURCE_GROUP_BITSTREAM})
source_group( "CABAC" FILES ${SOURCE_GROUP_CABAC})
source_group( "Compression" FILES ${SOURCE_GROUP_COMPRESSION})
source_group( "Constraint" FILES ${SOURCE_GROUP_CONSTRAINT})
source_group( "Control" FILES ${SOURCE_GROUP_CONTROL})
source_group( "Data Structures" FILES ${SOURCE_GROUP_DATA_STRUCTURES})
source_group( "Extras" FILES ${SOURCE_GROUP_EXTRAS})

# Handle the strategies directory structure better in visual studio
if(MSVC)
  foreach(source IN LISTS SOURCE_GROUP_STRATEGIES)
      get_filename_component(source_path "${source}" PATH)
      string(REPLACE "src/" "" source_path_msvc "${source_path}")
      string(REPLACE "/" "\\" source_path_msvc "${source_path_msvc}")
      source_group("Optimization\\${source_path_msvc}" FILES "${source}")
  endforeach()
else()
  source_group( "Optimization" FILES ${SOURCE_GROUP_STRATEGIES})
endif()
source_group( "Optimization" FILES "src/strategyselector.c" "src/strategyselector.h")

source_group( "Reconstruction" FILES ${SOURCE_GROUP_RECON})
source_group( "Threading" FILES ${SOURCE_GROUP_THREADING})
source_group( "Threadwrapper" FILES ${SOURCE_GROUP_THREADWRAPPER})
source_group( "" FILES ${SOURCE_GROUP_TOPLEVEL})

# INSTALL

# ToDo: make configurable

install(FILES ${PROJECT_SOURCE_DIR}/src/uvg266.pc DESTINATION ${CMAKE_INSTALL_PREFIX}/share/pkgconfig)
install(TARGETS uvg266-bin DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
install(TARGETS uvg266 DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
if(BUILD_SHARED_LIBS) # Just add the lib to the bin directory for now
  if(MSVC)
    install(TARGETS uvg266 DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
  endif()
endif()
install(FILES ${PROJECT_SOURCE_DIR}/src/uvg266.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
install(FILES ${PROJECT_SOURCE_DIR}/doc/uvg266.1 DESTINATION ${CMAKE_INSTALL_PREFIX}/share/man)

# DIST

file(GLOB_RECURSE DIST_SOURCES RELATIVE ${PROJECT_SOURCE_DIR} "*.c" "*.h" "*.sh" "*.py" "*.md" "*.yml" "*.in" "placeholder.txt" "CMakeLists.txt" "doc/*" "greatest/*")
list(APPEND DIST_SOURCES ".clang-format" ".gitignore" ".gitmodules" "tests/tsan_suppressions.txt" ".travis-install.bash" "CREDITS" "Dockerfile" "docs.doxy" "LICENSE" "LICENSE.EXT.greatest")

add_custom_target(dist
                COMMAND tar -zcvf "${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz" --transform 's,^,${PROJECT_NAME}-${PROJECT_VERSION}/,' -- ${DIST_SOURCES}
                COMMENT "Make distribution ${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz"
                WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
            )

# TESTS

enable_testing()

if(MSVC OR MINGW OR MSYS)
  if(BUILD_SHARED_LIBS) 
    set(BUILD_TESTS OFF)
    message(INFO " Disable test building, fails in MSVC/MINGW/MSYS2 when building shared binaries")
  endif()
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/greatest/greatest.h" AND BUILD_TESTS)
  add_subdirectory( "tests/" )
  add_test( NAME Test_uvg266 COMMAND uvg266_tests )
endif()

if(NOT DEFINED MSVC)

  # Disable tests in XFAIL_TESTS env, usage: `XFAIL_TESTS="test_intra;test_external_symbols" cmake .`
  list(APPEND XFAIL "off")
  if(DEFINED ENV{XFAIL_TESTS})
    list(APPEND XFAIL $ENV{XFAIL_TESTS})
  endif()  

  if(NOT "test_external_symbols" IN_LIST XFAIL)
    add_test( NAME test_external_symbols COMMAND ${PROJECT_SOURCE_DIR}/tests/test_external_symbols.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()
  if(NOT "test_intra" IN_LIST XFAIL)
    add_test( NAME test_intra COMMAND ${PROJECT_SOURCE_DIR}/tests/test_intra.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()
  if(NOT "test_invalid_input" IN_LIST XFAIL)
    add_test( NAME test_invalid_input COMMAND ${PROJECT_SOURCE_DIR}/tests/test_invalid_input.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()
  if(NOT "test_monochrome" IN_LIST XFAIL)
    add_test( NAME test_monochrome COMMAND ${PROJECT_SOURCE_DIR}/tests/test_monochrome.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()
  if(NOT "test_lmcs" IN_LIST XFAIL)
    add_test( NAME test_lmcs COMMAND ${PROJECT_SOURCE_DIR}/tests/test_lmcs.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()
  if(NOT "test_tools" IN_LIST XFAIL)
    add_test( NAME test_tools COMMAND ${PROJECT_SOURCE_DIR}/tests/test_tools.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()
  if(NOT "test_slices" IN_LIST XFAIL)
    add_test( NAME test_slices COMMAND ${PROJECT_SOURCE_DIR}/tests/test_slices.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()
endif()