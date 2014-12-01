CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
MESSAGE ( STATUS "Including FLUKA in dagmciface library.")

SET ( CODE_SRCDIR      "fluka" )

SET ( CODE_SRC_FILES   "${CODE_SRCDIR}/fluka_funcs.cpp" )
SET ( CODE_PUB_HEADERS "${CODE_SRCDIR}/fluka_funcs.h" )

SET ( CODE_CONFIG       )

SET ( DI_PUB_HEADERS ${DI_PUB_HEADERS}
                     ${CODE_PUB_HEADERS}
                     ${CODE_CONFIG}
      PARENT_SCOPE)

SET ( DI_SRC_FILES ${DI_SRC_FILES} 
                   ${CODE_SRC_FILES} 
      PARENT_SCOPE)

#INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/${CODE_SRCDIR})

INCLUDE_DIRECTORIES( ${HDF5_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${CODE_SRC_DIR} )

ENABLE_LANGUAGE( Fortran )

SET(FLUDAG_LIBS
     dagmciface
     ${HDF5_hdf5_hl_LIBRARY} ${HDF5_hdf5_LIBRARY}
     ${PYNE_LIBS_DIR}/libpyne.so
     ${DAGMC_LIBRARIES} 
     gfortran
  )

# If making mainfludag, define FLUPRO)
# add_custom_command (TARGET mainfludag
# PRE_BUILD COMMAND 

ADD_EXECUTABLE ( mainfludag mainFluDAG.cpp )
INSTALL( TARGETS mainfludag 
  RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
# Order very important!
TARGET_LINK_LIBRARIES ( mainfludag ${FLUDAG_LIBS} ${FLUKA_DIR}/libflukahp.a gfortran)

ADD_SUBDIRECTORY(tests)

