macro(uvg_msvc_enable_incremental_link)
  foreach(flags
    CMAKE_C_FLAGS_DEBUG
    CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_RELWITHDEBINFO
  )
    string(REPLACE "/Zi" "" ${flags} "${${flags}}")
    string(REPLACE "/ZI" "" ${flags} "${${flags}}")
    string(APPEND ${flags} " /ZI")
  endforeach()

  foreach(flags
    CMAKE_EXE_LINKER_FLAGS_DEBUG
    CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO
    CMAKE_SHARED_LINKER_FLAGS_DEBUG
    CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO
  )
    string(APPEND ${flags} " /INCREMENTAL")
  endforeach()
endmacro()

macro(uvg_msvc_disable_crt_warnings)
  add_compile_definitions(
    _CRT_NONSTDC_NO_DEPRECATE
    _CRT_SECURE_NO_WARNINGS
  )
endmacro()

macro(uvg_enable_better_diagnostics)
  if(CMAKE_C_COMPILER_ID MATCHES ".*Clang")
    add_compile_options(
      $<$<COMPILE_LANGUAGE:C>:-fcolor-diagnostics>
      $<$<COMPILE_LANGUAGE:CXX>:-fcolor-diagnostics>
    )
  elseif(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    add_compile_options(
      $<$<COMPILE_LANGUAGE:C>:-fdiagnostics-color=always>
      $<$<COMPILE_LANGUAGE:CXX>:-fdiagnostics-color=always>
    )
  elseif(CMAKE_C_COMPILER_ID STREQUAL "MSVC" AND MSVC_VERSION GREATER 1900)
    add_compile_options(
      $<$<COMPILE_LANGUAGE:C>:/diagnostics:column>
      $<$<COMPILE_LANGUAGE:CXX>:/diagnostics:column>
    )
  endif()
endmacro()

macro(uvg_setup_compiler_options)
  if(MSVC)
    uvg_msvc_enable_incremental_link()
    uvg_msvc_disable_crt_warnings()
  endif()

  uvg_enable_better_diagnostics()
endmacro()
