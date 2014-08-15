######################################################################
#
# CMakeLists.txt for dagmc Geant4
#
######################################################################
PROJECT( DAGSOLID )

MESSAGE ( STATUS "Building the DagSolid Library for  GEANT4." )
SET ( DAGSOLID_SRC  "DagSolid.cc" )
SET ( DAGSOLID_HEADERS "DagSolid.hh" )

# Geant4 include dirs
INCLUDE_DIRECTORIES( ${INCLUDE_DIRECTORIES} ${GEANT4_DIR}/include/Geant4 )

# Create DagSolid library
IF ( STATIC_LIB )
  ADD_LIBRARY( dagsolid STATIC ${DAGSOLID_SRC} )
ELSE ( STATIC_LIB )
  ADD_LIBRARY( dagsolid SHARED ${DAGSOLID_SRC} )
ENDIF ( STATIC_LIB )

# set the headers to install
SET_TARGET_PROPERTIES( dagsolid PROPERTIES PUBLIC_HEADER "${DAGSOLID_HEADERS}" )

# Put the library in Lib
INSTALL ( TARGETS dagsolid 
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}" )

# make the library
TARGET_LINK_LIBRARIES( dagsolid )
ADD_SUBDIRECTORY( tests )




