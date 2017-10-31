project(DAGMC)
cmake_minimum_required(VERSION 2.8)
enable_language(Fortran)

# Make the scripts in the "cmake" directory available to CMake
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(DAGMC_macros)

# Setup
dagmc_set_build_type()
dagmc_setup_build_options()

# Make sure we have HDF5 and MOAB
find_package(HDF5_DAGMC REQUIRED)
find_package(MOAB REQUIRED)

# More setup
dagmc_setup_flags()
dagmc_setup_rpath()

add_subdirectory(gtest)
add_subdirectory(pyne)
add_subdirectory(dagmc)
add_subdirectory(uwuw)
add_subdirectory(tools)

if (BUILD_TALLY)
  add_subdirectory(tally)
endif ()
if (BUILD_MCNP5 OR BUILD_MCNP6 OR BUILD_CI_TESTS)
  add_subdirectory(mcnp)
endif ()
if (BUILD_GEANT4)
  add_subdirectory(geant4)
endif ()
if (BUILD_FLUKA)
  add_subdirectory(fluka)
endif ()
