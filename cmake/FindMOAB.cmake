# This leverages HDF5_macro.cmake in order to find HDF5 libraries.

message("")

# Find MOAB cmake config file
# Only used to determine the location of the HDF5 with which MOAB was built
set(MOAB_SEARCH_DIRS)
file(GLOB MOAB_SEARCH_DIRS ${MOAB_SEARCH_DIRS} "${MOAB_DIR}/lib*/cmake/MOAB")
string(REPLACE "\n" ";" MOAB_SEARCH_DIRS "${MOAB_SEARCH_DIRS}")
find_path(MOAB_CMAKE_CONFIG
  NAMES MOABConfig.cmake
  PATHS ${MOAB_SEARCH_DIRS}
  NO_DEFAULT_PATH
)
if (MOAB_CMAKE_CONFIG)
  set(MOAB_CMAKE_CONFIG ${MOAB_CMAKE_CONFIG}/MOABConfig.cmake)
  message(STATUS "MOAB_CMAKE_CONFIG: ${MOAB_CMAKE_CONFIG}")
else ()
  message(FATAL_ERROR "Could not find MOAB. Set -DMOAB_DIR=<MOAB_DIR> when running cmake or use the $MOAB_DIR environment variable.")
endif ()

# Find HDF5
include(${MOAB_CMAKE_CONFIG})
include(HDF5_macro)
find_set_HDF5()

if(MSVC)
    set(BUILD_STATIC_LIBS TRUE)
    set(BUILD_SHARED_LIBS OFF)
endif()
# Find MOAB library (shared)
if (BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX})
  find_library(MOAB_LIBRARIES_SHARED
    NAMES MOAB
    HINTS ${MOAB_LIBRARY_DIRS}
    NO_DEFAULT_PATH
  )
  list(APPEND MOAB_LIBRARIES_SHARED)
endif ()

# Find MOAB library (static)
if (BUILD_STATIC_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
  find_library(MOAB_LIBRARIES_STATIC
    NAMES MOAB
    HINTS ${MOAB_LIBRARY_DIRS}
    NO_DEFAULT_PATH
  )
  list(APPEND MOAB_LIBRARIES_STATIC)
endif ()

message(STATUS "MOAB_INCLUDE_DIRS: ${MOAB_INCLUDE_DIRS}")
message(STATUS "MOAB_LIBRARY_DIRS: ${MOAB_LIBRARY_DIRS}")
message(STATUS "MOAB_LIBRARIES_SHARED: ${MOAB_LIBRARIES_SHARED}")
message(STATUS "MOAB_LIBRARIES_STATIC: ${MOAB_LIBRARIES_STATIC}")

if (MOAB_INCLUDE_DIRS AND (MOAB_LIBRARIES_SHARED OR NOT BUILD_SHARED_LIBS) AND
    (MOAB_LIBRARIES_STATIC OR NOT BUILD_STATIC_LIBS))
  message(STATUS "Found MOAB")
else ()
  message(FATAL_ERROR "Could not find MOAB")
endif ()

include_directories(${MOAB_INCLUDE_DIRS})
