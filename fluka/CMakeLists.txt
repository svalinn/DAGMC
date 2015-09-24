# Source file to build fluka lib and mainfludag
SET ( FLUKA_SRCDIR      ${CMAKE_CURRENT_SOURCE_DIR} )

# build the static lib
ADD_SUBDIRECTORY(lib)

# set the incs
INCLUDE_DIRECTORIES( ${HDF5_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${CODE_SRC_DIR} )

#enable fortran
ENABLE_LANGUAGE( Fortran )

# set the libraries
SET(FLUDAG_LIBS
     fludag
     dagmciface
     ${HDF5_hdf5_hl_LIBRARY} ${HDF5_hdf5_LIBRARY}
     pyne_dagmc
     ${DAG_LINK_LIBS} 
     gfortran
  )

# If we have chosen to have a static exectuable
IF(STATIC_LIB)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static")
ENDIF(STATIC_LIB)

INCLUDE_DIRECTORIES(mainfludag PUBLIC lib)

# the main exec
ADD_EXECUTABLE ( mainfludag mainFluDAG.cpp )
INSTALL( TARGETS mainfludag RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
# Order very important!
TARGET_LINK_LIBRARIES ( mainfludag ${FLUDAG_LIBS} ${FLUKA_DIR}/libflukahp.a gfortran)

ADD_SUBDIRECTORY(tests)

