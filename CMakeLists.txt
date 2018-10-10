project(DAGMC)
cmake_minimum_required(VERSION 2.8)
enable_language(Fortran)

# Make the scripts in the "cmake" directory available to CMake
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(DAGMC_macros)

dagmc_setup_build()
dagmc_setup_options()

find_package(MOAB REQUIRED)

dagmc_setup_flags()

if (BUILD_TESTS)
  enable_testing()
endif ()

add_subdirectory(src)

dagmc_make_configure_file()

message("")
