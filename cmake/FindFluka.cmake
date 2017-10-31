message("")

find_path(FLUKA_LIBRARIES
  NAMES libflukahp.a
  HINTS ${FLUKA_ROOT}
  HINTS ${FLUKA_DIR}
  PATHS ENV FLUKA_ROOT
  PATHS ENV FLUKA_DIR
  PATHS ENV LD_LIBRARY_PATH
  NO_DEFAULT_PATH
)
get_filename_component(FLUKA_LIBRARIES ${FLUKA_LIBRARIES} ABSOLUTE)
set(FLUKA_LIBRARIES ${FLUKA_LIBRARIES}/libflukahp.a gfortran)

message(STATUS "FLUKA_LIBRARIES: ${FLUKA_LIBRARIES}")

if (FLUKA_LIBRARIES)
  message(STATUS "Found Fluka")
else ()
  messate(FATAL_ERROR "Could not find Fluka")
endif ()

message("")
