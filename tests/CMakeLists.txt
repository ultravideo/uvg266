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
