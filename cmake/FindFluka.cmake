find_path(FLUKA_LIBRARIES
  NAMES libflukahp.a
  HINTS ${FLUKA_DIR}
  PATHS ENV FLUKA_DIR
  NO_DEFAULT_PATH
)
if (FLUKA_LIBRARIES)
  get_filename_component(FLUKA_LIBRARIES ${FLUKA_LIBRARIES} ABSOLUTE)
endif ()

set(FLUKA_LIBRARIES ${FLUKA_LIBRARIES}/libflukahp.a gfortran)

message(STATUS "FLUKA_LIBRARIES: ${FLUKA_LIBRARIES}")

if (FLUKA_LIBRARIES)
  message(STATUS "Found Fluka")
else ()
  message(FATAL_ERROR "Could not find Fluka")
endif ()
