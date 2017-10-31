# Find MOAB include directory
find_path(MOAB_INCLUDE_DIRS
  NAMES MBiMesh.hpp
  HINTS ${MOAB_ROOT}
  HINTS ${MOAB_DIR}
  PATHS ENV MOAB_ROOT
  PATHS ENV MOAB_DIR
  PATHS ENV LD_LIBRARY_PATH
  PATH_SUFFIXES . include Include ../include ../Include
  NO_DEFAULT_PATH
)
get_filename_component(MOAB_INCLUDE_DIRS ${MOAB_INCLUDE_DIRS} ABSOLUTE)

# Find MOAB library (shared)
set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
find_library(MOAB_LIBRARIES_SHARED
  NAMES MOAB
  HINTS ${MOAB_ROOT}
  HINTS ${MOAB_DIR}
  PATHS ENV MOAB_ROOT
  PATHS ENV MOAB_DIR
  PATHS ENV LD_LIBRARY_PATH
  PATH_SUFFIXES . lib Lib ../lib ../Lib
  NO_DEFAULT_PATH
)
get_filename_component(MOAB_LIBRARIES_SHARED ${MOAB_LIBRARIES_SHARED} ABSOLUTE)

# Find MOAB library (static)
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
find_library(MOAB_LIBRARIES_STATIC
  NAMES MOAB
  HINTS ${MOAB_ROOT}
  HINTS ${MOAB_DIR}
  PATHS ENV MOAB_ROOT
  PATHS ENV MOAB_DIR
  PATHS ENV LD_LIBRARY_PATH
  PATH_SUFFIXES . lib Lib ../lib ../Lib
  NO_DEFAULT_PATH
)
get_filename_component(MOAB_LIBRARIES_STATIC ${MOAB_LIBRARIES_STATIC} ABSOLUTE)

message(STATUS "MOAB_INCLUDE_DIRS: ${MOAB_INCLUDE_DIRS}")
message(STATUS "MOAB_LIBRARIES_SHARED: ${MOAB_LIBRARIES_SHARED}")
message(STATUS "MOAB_LIBRARIES_STATIC: ${MOAB_LIBRARIES_STATIC}")

if (MOAB_INCLUDE_DIRS AND MOAB_LIBRARIES_SHARED AND MOAB_LIBRARIES_STATIC)
  message(STATUS "Found MOAB")
else ()
  message(FATAL_ERROR "Could not find MOAB")
endif ()

include_directories(${MOAB_INCLUDE_DIRS})
