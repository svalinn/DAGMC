# Try to find DAGMC
#
# Once done this will define
#
#  DAGMC_FOUND        - indicates that the package has been found by CMake
#  DAGMC_INCLUDE_DIRS - include directories for installed DAGMC headers
#  DAGMC_LIBRARIES    - location of installed DAGMC libraries
#  DAGMC_LINK_LIBS    - standard libraries installed with DAGMC, use to link applications against DAGMC

find_path(DAGMC_CMAKE_CONFIG NAMES DAGMCConfig.cmake
          HINTS ${DAGMC_ROOT}
          PATHS ENV LD_LIBRARY_PATH
          PATH_SUFFIXES lib Lib cmake cmake/DAGMC
          NO_DEFAULT_PATH)

message(STATUS "Found DAGMC in ${DAGMC_CMAKE_CONFIG}")

include(${DAGMC_CMAKE_CONFIG}/DAGMCConfig.cmake)
