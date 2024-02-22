# Try to find dd
#
# Once done this will define
#
#  dd_FOUND        - indicates that the package has been found by CMake
#  dd_INCLUDE_DIRS - include directories for installed dd headers
#  dd_LIBRARY_DIRS - location of installed dd libraries
#  dd_LIBRARIES    - set of libraries installed with dd, use to link applications against dd
#  dd_VERSION      - version of installed dd

find_path(dd_CMAKE_CONFIG NAMES ddConfig.cmake
          HINTS ${dd_ROOT} $ENV{dd_ROOT}
          PATHS ENV LD_LIBRARY_PATH
          PATHS ${DOUBLE_DOWN_DIR}
          PATH_SUFFIXES lib Lib cmake lib/cmake/ lib/cmake/dd
          NO_DEFAULT_PATH)

if (dd_CMAKE_CONFIG)
    message(STATUS "Found dd in ${dd_CMAKE_CONFIG}")
    include(${dd_CMAKE_CONFIG}/ddConfig.cmake)
    include(${dd_CMAKE_CONFIG}/ddConfigVersion.cmake)
    set(dd_VERSION ${PACKAGE_VERSION})
    message(STATUS "Found dd version ${dd_VERSION}")
else()
    message(FATAL_ERROR "Could not find dd")
endif()
