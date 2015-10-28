# Source file to build fluka lib and mainfludag
SET ( FLUKA_SRCDIR      ${CMAKE_CURRENT_SOURCE_DIR} )

# set the source files
SET ( FLUKA_SRC_FILES   "${FLUKA_SRCDIR}/fluka_funcs.cpp" )
SET ( FLUKA_PUB_HEADERS "${FLUKA_SRCDIR}/fluka_funcs.h" )


# include fluka lib as static archive
ADD_LIBRARY(fludag STATIC ${FLUKA_SRC_FILES})

# set the incs
INCLUDE_DIRECTORIES( ${HDF5_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${CODE_SRC_DIR} )


