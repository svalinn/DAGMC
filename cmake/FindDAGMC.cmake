# Try to find DAGMC
#
# Once done this will define
#
#  DAGMC_FOUND - system has DAGMC
#  DAGMC_INCLUDE_DIRS - the DAGMC include directory
#  DAGMC_LIBRARIES - Link these to use DAGMC

find_path(DAGMC_CMAKE_CONFIG NAMES DAGMCConfig.cmake
          HINTS ${DAGMC_ROOT}
          PATHS ENV LD_LIBRARY_PATH
          PATH_SUFFIXES lib Lib cmake cmake/DAGMC
          NO_DEFAULT_PATH)

message(STATUS "Found DAGMC in ${DAGMC_CMAKE_CONFIG}")

include(${DAGMC_CMAKE_CONFIG}/DAGMCConfig.cmake)
