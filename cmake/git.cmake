function(uvg_find_git)
  if(DEFINED UVG_GIT_EXECUTABLE)
    return()
  endif()

  message(STATUS "Looking for git")
  find_program(UVG_GIT_EXECUTABLE git)
  if(UVG_GIT_EXECUTABLE)
    message(STATUS "Looking for git - ${UVG_GIT_EXECUTABLE}")
    set(UVG_GIT_EXECUTABLE "${UVG_GIT_EXECUTABLE}" PARENT_SCOPE)
  else()
    message(STATUS "Looking for git - not found")
  endif()
endfunction()

macro(uvg_maybe_override_project_version)
  if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    uvg_find_git()
    if(UVG_GIT_EXECUTABLE)
      execute_process(
        COMMAND ${UVG_GIT_EXECUTABLE} describe --tags --dirty
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
      )

      if(result EQUAL 0 AND output MATCHES "v[0-9]+\\.[0-9]+\\.[0-9]+.*")
        # The output starts with something that looks like a version tag,
        # so strip the leading 'v' and use that as the project version.
        string(STRIP "${output}" output)
        string(SUBSTRING "${output}" 1 -1 output)
        set(PROJECT_VERSION "${output}")
        message(STATUS "Project version -> ${PROJECT_VERSION}")
      endif()
    endif()
  endif()
endmacro()

function(uvg_maybe_update_submodules)
  if(NOT EXISTS "${CMAKE_SOURCE_DIR}/.git" OR NOT UVG_SYNC_SUBMODULES)
    return()
  endif()

  uvg_find_git()
  if(NOT UVG_GIT_EXECUTABLE)
    return()
  endif()

  message(STATUS "Updating submodules")
  execute_process(
    COMMAND ${UVG_GIT_EXECUTABLE} submodule update --init --recursive
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(WARNING "Failed to update submodules!")
  endif()
endfunction()
