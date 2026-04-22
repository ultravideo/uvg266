macro(uvg_try_compile file flags result)
  try_compile(
    success
    "${CMAKE_BINARY_DIR}"
    "${file}"
    COMPILE_DEFINITIONS "${flags}"
  )
  if(success)
    set(${result} 1)
  else()
    set(${result} 0)
  endif()
endmacro()
