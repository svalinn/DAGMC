find_path(FLUKA_LIBDIR 
    NAMES libfluka.a libflukahp.a
    HINTS ${FLUKA_DIR}
    PATHS ENV FLUKA_DIR
    NO_DEFAULT_PATH
)

if (NOT FLUKA_LIBDIR)
    message(FATAL_ERROR "Could not find Fluka. Please set FLUKA_DIR or ENV{FLUKA_DIR}.")
endif()

if (EXISTS "${FLUKA_LIBDIR}/libfluka.a")
    # FLUKA.CERN
    set(FLUKA_LIBRARIES 
        ${FLUKA_LIBDIR}/libfluka.a
        ${FLUKA_LIBDIR}/libgeometry.a
        ${FLUKA_LIBDIR}/libdata.a
        ${FLUKA_LIBDIR}/libmath.a
        ${FLUKA_LIBDIR}/libtool.a
        gfortran
        z
        pthread
    )
    message(STATUS "Found FLUKA installation)")

elseif (EXISTS "${FLUKA_LIBDIR}/libflukahp.a")
    # Old format: Single high-precision library
    set(FLUKA_LIBRARIES 
        ${FLUKA_LIBDIR}/libflukahp.a 
        gfortran
    )
    message(STATUS "FLUKA installation")

else()
    # Failsafe in case the directory was found but files are missing
    message(FATAL_ERROR "Found FLUKA_LIBDIR at ${FLUKA_LIBDIR}, but required libraries are missing.")
endif()

message(STATUS "FLUKA_LIBRARIES: ${FLUKA_LIBRARIES}")
