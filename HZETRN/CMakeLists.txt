#############################################################3
# HZETRN Build
#############################################################3
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET ( CODE_SRCDIR      "HZETRN" )

# SET ( CODE_SRC_FILES   "${CODE_SRCDIR}/hzetrn_paths.cpp" )
# SET ( CODE_PUB_HEADERS "${CODE_SRCDIR}/fluka_funcs.h" )

SET ( CODE_CONFIG       )

######################################################
# To add this project's files to dagmciface, uncomment
# SET ( DI_PUB_HEADERS ${DI_PUB_HEADERS}
#                    ${CODE_PUB_HEADERS}
#                    ${CODE_CONFIG}
#     PARENT_SCOPE)

######################################################
# To add this project's files to dagmciface, uncomment
# SET ( DI_SRC_FILES ${DI_SRC_FILES} 
#                    ${CODE_SRC_FILES} 
#      PARENT_SCOPE)

# INCLUDE_DIRECTORIES( ${HDF5_INCLUDE_DIRS} )
# INCLUDE_DIRECTORIES( ${CODE_SRC_DIR} )

SET(HZETRN_LIBS
     ${DAGMC_LIBRARIES} 
  )

ADD_EXECUTABLE ( hzetrn hzetrn_paths.cc )
INSTALL( TARGETS hzetrn
  RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
# Order very important!
TARGET_LINK_LIBRARIES ( hzetrn ${HZETRN_LIBS} )

# ADD_SUBDIRECTORY(tests)

