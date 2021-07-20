find_path(FLUKA_LIBRARIES
  NAMES libfluka.a
  HINTS ${FLUKA_DIR} ${FLUKA_DIR}/lib
  PATHS ENV FLUKA_DIR
  NO_DEFAULT_PATH
)
if (FLUKA_LIBRARIES)
  get_filename_component(FLUKA_LIBRARIES ${FLUKA_LIBRARIES} ABSOLUTE)
endif ()

set(FLUKA_LIBRARIES ${FLUKA_LIBRARIES}/libfluka.a gfortran)

message(STATUS "FLUKA_LIBRARIES: ${FLUKA_LIBRARIES}")

if (FLUKA_LIBRARIES)
  message(STATUS "Found Fluka")
else ()
  message(FATAL_ERROR "Could not find Fluka")
endif ()
