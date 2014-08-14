######################################################################
#
# CMakeLists.txt for dagmc Geant4
#
######################################################################


MESSAGE ( STATUS "Building the DagSolid Library for  GEANT4." )
SET ( DAGSOLID_SRCDIR  "DagSolid.cc" )
SET ( DAGSOLID_HEADERS "DagSolid.hh" )

#INCLUDE_DIRECTORIES( "${DAG_GEANT_SRCDIR}" )

# Geant4 include dirs
INCLUDE_DIRECTORIES( ${INCLUDE_DIRECTORIES} ${GEANT4_DIR}/include/Geant4 )

# Create DagSolid library
IF ( STATIC_LIB )
  ADD_LIBRARY( dagsolid STATIC ${DAGSOLID_SRCDIR} )
ELSE ( STATIC_LIB )
  ADD_LIBRARY( dagsolid SHARED ${DAGSOLID_SRCDIR} )
ENDIF ( STATIC_LIB )

# set the headers to install
SET_TARGET_PROPERTIES( dagsolid PROPERTIES PUBLIC_HEADER "${DAGSOLID_HEADERS}" )

# Put the library in Lib
INSTALL ( TARGETS dagsolid 
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}" )

# make the library
TARGET_LINK_LIBRARIES(dagsolid)




