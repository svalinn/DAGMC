# Figure out which physics codes to build
if (BUILD_ALL)
  set(BUILD_MCNP5 ON)
  set(BUILD_MCNP6 ON)
  set(BUILD_GEANT4 ON)
  set(BUILD_FLUKA ON)
endif ()

# DAG-MCNP and CI require the dagtally library
if (BUILD_MCNP5 OR BUILD_MCNP6 OR BUILD_CI_TESTS)
  set(BUILD_TALLY ON)
endif ()

# Required libraries
add_subdirectory(gtest)
add_subdirectory(dagmc)
add_subdirectory(pyne)
add_subdirectory(uwuw)

# Non-required libraries
if (BUILD_TALLY)
  add_subdirectory(tally)
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

# Tools
add_subdirectory(build_obb)
add_subdirectory(make_watertight)
