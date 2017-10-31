project(DAGMC)
cmake_minimum_required(VERSION 2.8)
enable_language(Fortran)

# Make the scripts in the "cmake" directory available to CMake
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(DAGMC_macros)

set(INSTALL_BIN_DIR     bin)
set(INSTALL_LIB_DIR     lib)
set(INSTALL_INCLUDE_DIR include)
set(INSTALL_TESTS_DIR   tests)
set(INSTALL_TOOLS_DIR   tools)
set(INSTALL_SHARE_DIR   share)

dagmc_set_build_type()
dagmc_setup_build_options()

find_package(HDF5_DAGMC REQUIRED)
find_package(MOAB REQUIRED)

dagmc_setup_flags()

if (BUILD_TESTS)
  add_subdirectory(gtest)
endif ()

add_subdirectory(src)
