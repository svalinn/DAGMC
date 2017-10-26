project(DAGMC)

cmake_minimum_required(VERSION 2.8)

enable_language(Fortran)

# Make the scripts in the "cmake" directory available to CMake
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(DAGMC_macros)

# Setup
dagmc_set_build_type()
dagmc_setup_compile_flags()
dagmc_setup_link_flags()

# Make sure we have HDF5
find_package(HDF5 REQUIRED)
message(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
include_directories(${HDF5_INCLUDE_DIRS})

# Make sure we have MOAB
find_package(MOAB REQUIRED)
include_directories(${MOAB_INCLUDE_DIRS})

dagmc_setup_build_options()

add_subdirectory(gtest)
add_subdirectory(pyne)
add_subdirectory(dagmc)
add_subdirectory(uwuw)
add_subdirectory(tools)

if (BUILD_TALLY)
  add_subdirectory(tally)
endif ()

if (BUILD_MCNP5 OR BUILD_MCNP6)
  add_subdirectory(mcnp)
endif ()

if (BUILD_GEANT4)
  add_subdirectory(geant4)
endif ()

if (BUILD_FLUKA)
  add_subdirectory(fluka)
endif ()

# Build tests that wouldn't normally be built for the CI
if (BUILD_CI_TESTS)
  if (NOT BUILD_MCNP5 AND NOT BUILD_MCNP6)
    add_subdirectory(mcnp)
  endif ()
endif ()
