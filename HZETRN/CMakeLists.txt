#############################################################3
# HZETRN Build
#############################################################3
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
# MESSAGE ( STATUS "Including FLUKA in dagmciface library.")

SET ( CODE_SRCDIR      "HZETRN" )

# SET ( CODE_SRC_FILES   "${CODE_SRCDIR}/hzetrn_paths.cpp" )
# SET ( CODE_PUB_HEADERS "${CODE_SRCDIR}/fluka_funcs.h" )

SET ( CODE_CONFIG       )

SET ( DI_PUB_HEADERS ${DI_PUB_HEADERS}
                     ${CODE_PUB_HEADERS}
                     ${CODE_CONFIG}
      PARENT_SCOPE)

SET ( DI_SRC_FILES ${DI_SRC_FILES} 
                   ${CODE_SRC_FILES} 
      PARENT_SCOPE)

# INCLUDE_DIRECTORIES( ${HDF5_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${CODE_SRC_DIR} )

# ENABLE_LANGUAGE( Fortran )

SET(HZETRNS_LIBS
     # dagmciface
     # ${HDF5_LIBRARIES}
     ${DAGMC_LIBRARIES} 
     # gfortran
  )

# If making mainfludag, define FLUPRO)
# add_custom_command (TARGET mainfludag
# PRE_BUILD COMMAND 

# ADD_EXECUTABLE ( mainhzetrn mainFluDAG.cpp )
ADD_EXECUTABLE ( hzetrn hzetrn_paths.cc )
INSTALL( TARGETS  
  RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
# Order very important!
TARGET_LINK_LIBRARIES ( hzetrn dagmciface ${DAGMC_LIBS} )

# ADD_SUBDIRECTORY(tests)

