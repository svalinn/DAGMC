######################################################################
#
# CMakeLists.txt for dagmc Geant4
#
######################################################################


MESSAGE ( STATUS "Building the DagSolid Library for  GEANT4.")
SET ( DAG_GEANT_SRCDIR "geant4" )
SET ( DAGMCIFACE_GEANT_SRCDIR "${DAG_GEANT_SRCDIR}/DagSolid.cc" )
SET ( DAGMCIFACE_GEANT_HEADERS "${DAG_GEANT_SRCDIR}/DagSolid.hh")

SET ( DAGSOLID_LIB ${DAGMCIFACE_GEANT_SRCDIR} )
# all headers
SET( dagmciface_headers ${dagmciface_headers} ${DAGMCIFACE_GEANT_HEADERS} )

INCLUDE_DIRECTORIES( "${DAG_GEANT_SRCDIR}" )

# Geant4 include dirs
INCLUDE_DIRECTORIES( ${INCLUDE_DIRECTORIES} ${GEANT4_DIR}/include/Geant4)

# Create DAGMC Tally library                                                                                                        
ADD_LIBRARY(dagsolid SHARED ${DAGMCIFACE_GEANT_SRCDIR} )

TARGET_LINK_LIBRARIES(dagsolid)



