######################################################################
#
# CMakeLists.txt for dagmc Geant4
#
######################################################################


MESSAGE ( STATUS "Building the DagSolid Library for  GEANT4." )
SET ( DAG_GEANT_SRCDIR "geant4" )
SET ( DAGMCIFACE_GEANT_SRCDIR "${DAG_GEANT_SRCDIR}/DagSolid.cc" )
SET ( DAGMCIFACE_GEANT_HEADERS "${DAG_GEANT_SRCDIR}/DagSolid.hh" )

SET ( DAGSOLID_LIB ${DAGMCIFACE_GEANT_SRCDIR} )
# all headers
SET( dagmciface_headers ${dagmciface_headers} ${DAGMCIFACE_GEANT_HEADERS} )

INCLUDE_DIRECTORIES( "${DAG_GEANT_SRCDIR}" )

# Geant4 include dirs
INCLUDE_DIRECTORIES( ${INCLUDE_DIRECTORIES} ${GEANT4_DIR}/include/Geant4 )

# Create DagSolid library                                                                                                        
IF ( STATIC_LIB )
  ADD_LIBRARY( dagsolid STATIC ${DAGMCIFACE_GEANT_SRCDIR} )
ELSE ( STATIC_LIB )
  ADD_LIBRARY( dagsolid SHARED ${DAGMCIFACE_GEANT_SRCDIR} )
ENDIF ( STATIC_LIB )

# set the headers to install
SET_TARGET_PROPERTIES( dagsolid PROPERTIES PUBLIC_HEADER "${DAGMCIFACE_GEANT_HEADERS}" )

# Put the library in Lib
INSTALL ( TARGETS dagsolid 
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}" )

# put the header in inclucde
TARGET_LINK_LIBRARIES(dagsolid)




