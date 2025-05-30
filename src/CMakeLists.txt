# Gtest
if (BUILD_CI_TESTS)
  set(BUILD_TESTS ON)
endif ()
if (BUILD_TESTS)
  add_subdirectory(gtest)
endif ()

# Required libraries
add_subdirectory(dagmc)

# UWUW
if (BUILD_MCNP5 OR BUILD_MCNP6 OR BUILD_GEANT4 OR BUILD_FLUKA OR BUILD_CI_TESTS)
  set(BUILD_UWUW ON)
endif ()
if (BUILD_UWUW)
  add_subdirectory(uwuw)

  if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(COMMAND ${GIT_EXECUTABLE} "submodule" "update" "--init" "--recursive"
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL 0)
      message(FATAL_ERROR "git submodule update --init --recursive failed with \
        ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()

  # Check to see if submodules exist (by checking one)
  if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/pyne/pyne/readme.rst")
    message(FATAL_ERROR "The git submodules were not downloaded! GIT_SUBMODULE was \
      turned off or failed. Please update submodules and try again.")
  endif()
  add_subdirectory(pyne)
endif ()

# dagtally library (required for DAG-MCNP and CI)
if (BUILD_MCNP5 OR BUILD_MCNP6 OR BUILD_CI_TESTS)
  set(BUILD_TALLY ON)
endif ()
if (BUILD_TALLY)
  add_subdirectory(tally)
endif ()

# Tools
if (BUILD_BUILD_OBB)
  add_subdirectory(build_obb)
endif ()
if (BUILD_MAKE_WATERTIGHT)
  add_subdirectory(make_watertight)
endif ()

# Physics code interfaces
if (BUILD_MCNP5 OR BUILD_MCNP6 OR BUILD_CI_TESTS)
  add_subdirectory(mcnp)
endif ()
if (BUILD_GEANT4)
  add_subdirectory(geant4)
endif ()
if (BUILD_FLUKA)
  add_subdirectory(fluka)
endif ()

if (BUILD_OVERLAP_CHECK)
  add_subdirectory(overlap_check)
endif ()
