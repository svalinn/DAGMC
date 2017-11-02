project(DAGMC)
cmake_minimum_required(VERSION 2.8)
enable_language(Fortran)

# Make the scripts in the "cmake" directory available to CMake
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(DAGMC_macros)

dagmc_setup_build()
dagmc_setup_options()

find_package(HDF5_DAGMC REQUIRED)
find_package(MOAB REQUIRED)

dagmc_setup_flags()
dagmc_setup_rpath()

if (BUILD_TESTS)
  add_subdirectory(gtest)
endif ()

add_subdirectory(src)

message("")
