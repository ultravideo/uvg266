############################################################
# Tests implemented with greatest
############################################################

function(add_greatest_tests)
  add_executable(run-greatest
    coeff_sum_tests.c
    dct_tests.c
    # FIXME
    #inter_recon_bipred_tests.c
    intra_sad_tests.c
    mts_tests.c
    mv_cand_tests.c
    sad_tests.c
    satd_tests.c
    speed_tests.c
    test_strategies.c
    tests_main.c
  )
  target_link_libraries(run-greatest PRIVATE libuvg266_static)
  target_include_directories(run-greatest PRIVATE ${PROJECT_SOURCE_DIR})
  if(MSVC)
    target_include_directories(run-greatest SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/src/threadwrapper/include)
  endif()
  add_test(NAME "Test_uvg266" COMMAND run-greatest)
endfunction()

if(EXISTS "${PROJECT_SOURCE_DIR}/greatest/greatest.h")
  add_greatest_tests()
else()
  message(WARNING "Couldn't find greatest.h - some tests will not be built")
endif()

############################################################
# Tests implemented with shell scripts
############################################################

function(add_shell_tests)
  # Disable tests in XFAIL_TESTS env, usage: `XFAIL_TESTS="test_intra;test_external_symbols" cmake .`
  list(APPEND XFAIL "off")
  if(DEFINED ENV{XFAIL_TESTS})
    list(APPEND XFAIL $ENV{XFAIL_TESTS})
  endif()

  foreach(TEST
    test_external_symbols
    test_mtt_intra
    test_mtt_inter
    test_mtt_both
    test_intra
    test_invalid_input
    test_monochrome
    test_lmcs
    test_tools
    test_slices
    test_cabac_state
  )
    if(NOT "${TEST}" IN_LIST XFAIL)
      add_test(
        NAME "${TEST}"
        COMMAND ${PROJECT_SOURCE_DIR}/tests/${TEST}.sh
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests
      )
    endif()
  endforeach()
endfunction()

if(UNIX)
  add_shell_tests()
endif()

############################################################
# Target: run-catch2
############################################################

function(try_build_tests)
  find_package(Catch2 QUIET)

  if(NOT Catch2_FOUND)
    message(WARNING "Catch2 not found - tests will not be built")
    return()
  endif()

  add_executable(run-catch2
    run-catch2.cpp
    strategies/transpose.test.cpp
  )
  target_link_libraries(run-catch2
    PRIVATE
      libuvg266_static
      Catch2::Catch2
  )
  target_compile_definitions(run-catch2 PRIVATE UVG_NO_BENCHMARKS=1)

  include(Catch)
  catch_discover_tests(run-catch2)
endfunction()

try_build_tests()

############################################################
# Target: run-benchmark
############################################################

function(try_build_benchmarks)
  find_package(benchmark QUIET)
  find_package(CLI11 QUIET)

  if(NOT benchmark_FOUND)
    message(WARNING "Google Benchmark not found - benchmarks will not be built")
    return()
  endif()

  if(NOT CLI11_FOUND)
    message(WARNING "CLI11 not found - benchmarks will not be built")
    return()
  endif()

  add_executable(run-benchmark
    run-benchmark.cpp
    strategies/transpose.test.cpp
  )
  target_link_libraries(run-benchmark
    PRIVATE
      libuvg266_static
      benchmark::benchmark
      CLI11::CLI11
  )
  target_compile_definitions(run-benchmark PRIVATE UVG_NO_TESTS=1)
endfunction()

try_build_benchmarks()
