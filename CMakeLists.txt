cmake_minimum_required(VERSION 3.12)

project(uvg266
LANGUAGES C CXX
HOMEPAGE_URL https://github.com/ultravideo/uvg266
DESCRIPTION "An open-source VVC encoder licensed under 3-clause BSD"
VERSION 0.8.1 )

option(BUILD_SHARED_LIBS "Build using shared uvg266 library" ON)

option(BUILD_TESTS "Build tests" ON)

include(GNUInstallDirs) #Helps to define correct distro specific install directories

set(UVG266_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}" CACHE PATH "uvg266 library install path")
set(UVG266_INSTALL_BINDIR "${CMAKE_INSTALL_BINDIR}" CACHE PATH "uvg266 binary install path")
set(UVG266_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}" CACHE PATH "uvg266 include install path")
set(UVG266_INSTALL_MANDIR "${CMAKE_INSTALL_MANDIR}/man1" CACHE PATH "uvg266 manual page file install path")

# https://www.kitware.com/cmake-and-the-default-build-type/
# Set a default build type if none was specified
set(UVG_DEFAULT_BUILD_TYPE "RelWithDebInfo")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type specified, setting to '${UVG_DEFAULT_BUILD_TYPE}'.")
  set(CMAKE_BUILD_TYPE "${UVG_DEFAULT_BUILD_TYPE}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


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
    # Check git hash and fetch tag
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_HEAD_OK
                    OUTPUT_VARIABLE GIT_HEAD)
    if(GIT_HEAD_OK EQUAL "0")
        string(SUBSTRING ${GIT_HEAD} 0 30 GIT_TAG_LONG)    
        execute_process(COMMAND ${GIT_EXECUTABLE} name-rev --tags --name-only ${GIT_TAG_LONG}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_TAG_OK
                        OUTPUT_VARIABLE GIT_TAG)
        string(SUBSTRING ${GIT_TAG} 0 9 GIT_TAG_STRIP)
      
        # If tag is not defined, add part of the commit hash to the version
        if(GIT_TAG_OK EQUAL "0" AND GIT_TAG_STRIP STREQUAL "undefined")
          string(SUBSTRING ${GIT_HEAD} 0 7 GIT_TAG_SHORT)
          set(PROJECT_VERSION ${PROJECT_VERSION}-${GIT_TAG_SHORT})
          message(INFO " No tag detected, version changed to ${PROJECT_VERSION}")
        endif()
    endif()    
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/greatest/greatest.h")
    message(WARNING "The submodule greatest was not loaded, some tests may fail")
endif()

# Grab <year>-<month>-<day> timestamp for debug purposes
string(TIMESTAMP CMAKE_BUILD_DATE %Y-%m-%d)

set(UVG_COMPILER_VERSION "${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")

if(MSVC)
    if(MSVC_VERSION LESS 1800)
        set(UVG_COMPILER_VERSION "VS")
    elseif(MSVC_VERSION LESS 1900)
        set(UVG_COMPILER_VERSION "VS2013")
    elseif(MSVC_VERSION LESS 1910)
        set(UVG_COMPILER_VERSION "VS2015")
    elseif(MSVC_VERSION LESS 1920)
        set(UVG_COMPILER_VERSION "VS2017")
    elseif(MSVC_VERSION LESS 1930)
        set(UVG_COMPILER_VERSION "VS2019")
    else()
        set(UVG_COMPILER_VERSION "VS2022")
    endif()
endif()

# Set compiler info to print at runtime
set(UVG_COMPILER_STRING "${UVG_COMPILER_VERSION}")

# Apply dynamic info to the config files
configure_file("${PROJECT_SOURCE_DIR}/src/uvg266.pc.in" "${PROJECT_SOURCE_DIR}/src/uvg266.pc" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/src/version.h.in" "${PROJECT_SOURCE_DIR}/src/version.h" @ONLY)

# Add all sources in src/ base
file(GLOB LIB_SOURCES RELATIVE ${PROJECT_SOURCE_DIR} "src/*.h" "src/*.c")

# We don't want CLI main in the library
list(REMOVE_ITEM LIB_SOURCES "src/encmain.c" "src/cli.c" "src/cli.h" "src/yuv_io.c" "src/yuv_io.h")

# Add also all the strategies
file(GLOB_RECURSE LIB_SOURCES_STRATEGIES RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/*.h" "src/strategies/*.c")

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
  add_definitions(-DWIN32_LEAN_AND_MEAN -D_WIN32 -DWIN32 -DWIN64 -D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_DEPRECATE)
endif()

if(BUILD_SHARED_LIBS)
  list( APPEND CMAKE_INSTALL_RPATH "${UVG266_INSTALL_LIBDIR}" "./" "../lib" )
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
file(GLOB LIB_SOURCES_STRATEGIES_SSE42 RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/sse42/*.c")

set(CLI_SOURCES "src/encmain.c" "src/cli.c" "src/cli.h" "src/yuv_io.c" "src/yuv_io.h")

# Add the getopt and pthread for visual studio
if(MSVC) 
  list(APPEND CLI_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/extras/getopt.c ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/pthread.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/semaphore.cpp)
  
  string(REPLACE "/Zi" "/ZI" CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
  string(REPLACE "/Zi" "/ZI" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
  string(REPLACE "/Zi" "/ZI" CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")
  string(REPLACE "/Zi" "/ZI" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

  target_link_options(uvg266 PUBLIC "/INCREMENTAL")
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
    set_property( SOURCE ${LIB_SOURCES_STRATEGIES_SSE42} APPEND PROPERTY COMPILE_FLAGS "-msse4.2" )
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
file(GLOB SOURCE_GROUP_DATA_STRUCTURES RELATIVE ${PROJECT_SOURCE_DIR} "src/cu.*" "src/image.*" "src/imagelist.*" "src/videoframe.*" "src/hashmap.*")
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

install(FILES ${PROJECT_SOURCE_DIR}/src/uvg266.pc DESTINATION ${UVG266_INSTALL_LIBDIR}/pkgconfig)
install(TARGETS uvg266-bin DESTINATION ${UVG266_INSTALL_BINDIR})
install(TARGETS uvg266 ARCHIVE DESTINATION "${UVG266_INSTALL_LIBDIR}" LIBRARY DESTINATION "${UVG266_INSTALL_LIBDIR}" RUNTIME DESTINATION "${UVG266_INSTALL_BINDIR}")
install(FILES ${PROJECT_SOURCE_DIR}/src/uvg266.h DESTINATION ${UVG266_INSTALL_INCLUDEDIR})
install(FILES ${PROJECT_SOURCE_DIR}/doc/uvg266.1 DESTINATION ${UVG266_INSTALL_MANDIR})

IF(UNIX)
# DIST

set(GIT_LS_TREE_OK "1")

# By default grab the list of files in the git repo
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    execute_process(COMMAND ${GIT_EXECUTABLE} ls-tree --name-only -r HEAD
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_LS_TREE_OK
                    OUTPUT_VARIABLE GIT_LS_TREE)
    if(GIT_LS_TREE_OK EQUAL "0")
        string(REGEX REPLACE "\n" ";" GIT_LS_TREE "${GIT_LS_TREE}")
        string(REGEX REPLACE "\r" "" GIT_LS_TREE "${GIT_LS_TREE}")
        list(APPEND DIST_SOURCES ${GIT_LS_TREE})
    endif()
endif()
if(NOT GIT_LS_TREE_OK EQUAL "0")
    file(GLOB_RECURSE DIST_SOURCES RELATIVE ${PROJECT_SOURCE_DIR} "src/*.c" "src/*.h" "tests/*.sh" "tools/*.sh" "tools/*.py" ".github/*.yml" "src/*.in" "placeholder.txt" "CMakeLists.txt" "doc/*" "examples/*" "rdcost-weight-tool/*" "greatest/*.h" "greatest/*.md")
    list(APPEND DIST_SOURCES ".clang-format" ".gitignore" ".gitmodules" "tests/tsan_suppressions.txt" ".travis-install.bash" "CREDITS" "Dockerfile" "docs.doxy" ".gitlab-ci.yml" "LICENSE" "LICENSE.EXT.greatest" "README.md")
endif()

add_custom_target(dist
                COMMAND echo \"Writing log to ${PROJECT_SOURCE_DIR}/dist.log\"
                 && tar -zcvf "${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz" --transform 's,^,${PROJECT_NAME}-${PROJECT_VERSION}/,' -- ${DIST_SOURCES} > dist.log 2>&1 || { echo \"\\033[0;31mfailed to pack ${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz, check ${PROJECT_SOURCE_DIR}/dist.log.\\033\[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                COMMENT "Make distribution ${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz"
                WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
                BYPRODUCTS ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz
            )

# DISTCHECK

set(TEMP_DISTCHECK_DIR "_distcheck")

add_custom_target(distcheck
                  COMMAND echo \"Writing log to ${PROJECT_SOURCE_DIR}/distcheck.log\"
                  && cd ${PROJECT_SOURCE_DIR}
                  && mkdir -p ${TEMP_DISTCHECK_DIR}
                  && cd ${TEMP_DISTCHECK_DIR}
                  && tar -zxf ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz > ${PROJECT_SOURCE_DIR}/distcheck.log || { echo \"\\033[0;31mfailed to unpack ${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz.\\033\[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                  && echo \"\\033[0;32mFile unpack ok\\033[m\"
                  && cd ${PROJECT_NAME}-${PROJECT_VERSION}
                  && mkdir -p build
                  && cd build                  
                  && cmake -DCMAKE_INSTALL_PREFIX=./ -DBUILD_SHARED_LIBS=OFF -G "Unix Makefiles" .. >> ${PROJECT_SOURCE_DIR}/distcheck.log || { echo \"\\033[0;31mcmake failed to configure.\\033[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                  && echo \"\\033[0;32mCMake configure ok\\033[m\"
                  && make -j >> ${PROJECT_SOURCE_DIR}/distcheck.log || { echo \"\\033[0;31mmake failed.\\033[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                  && echo \"\\033[0;32mMake ok\\033[m\"
                  # Full tests might be too demanding to run, enable with parameter?
                  #&& make test || (echo \"\\e[0;31mmake test failed.\\033[m\" && false)
                  && tests/uvg266_tests >> ${PROJECT_SOURCE_DIR}/distcheck.log 2>&1 || { echo \"\\033[0;31mtests failed.\\033[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                  && echo \"\\033[0;32mTests ok\\033[m\"
                  && make install >> ${PROJECT_SOURCE_DIR}/distcheck.log || { echo \"\\033[0;31mmake install failed.\\033[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                  && echo \"\\033[0;32mInstall ok\\033[m\"
                  && ${CMAKE_INSTALL_BINDIR}/uvg266 --help >> ${PROJECT_SOURCE_DIR}/distcheck.log || { echo \"\\033[0;31muvg266 binary failed to run.\\033[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                  && echo \"\\033[0;32m${CMAKE_INSTALL_BINDIR}/uvg266 ok\\033[m\"
                  && make clean >> ${PROJECT_SOURCE_DIR}/distcheck.log || { echo \"\\033[0;31mmake clean failed.\\033[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                  && echo \"\\033[0;32mmake clean ok\\033[m\"
                  && cd ${PROJECT_SOURCE_DIR}
                  && rm -rf "${PROJECT_SOURCE_DIR}/${TEMP_DISTCHECK_DIR}"
                  && echo \"\\033[0;32m==============================================================\\033[m\"
                  && echo \"\\033[0;32m${PROJECT_NAME}-${PROJECT_VERSION} archives ready for distribution:\\033[m\"
                  && echo \"\\033[0;32m${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz\\033[m\"
                  && echo \"\\033[0;32m==============================================================\\033[m\"
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  DEPENDS ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz
                  COMMENT "Checking ${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz.."
                  )
endif() #Unix

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
  if(NOT "test_mtt" IN_LIST XFAIL)
    add_test( NAME test_mtt COMMAND ${PROJECT_SOURCE_DIR}/tests/test_mtt.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
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
  if(NOT "test_cabac_state" IN_LIST XFAIL)
    add_test( NAME test_cabac_state COMMAND ${PROJECT_SOURCE_DIR}/tests/test_cabac_state.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()
endif()
