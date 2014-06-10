######################################################################
#
# CMakeLists.txt for dagmc Geant4
#
######################################################################

MESSAGE ( STATUS "Including GEANT4 in dagmciface library.")
SET ( DAG_GEANT_SRCDIR "geant4/src" )
SET ( DAGMCIFACE_GEANT_SRC "${DAG_GEANT_SRCDIR}/DagSolid.cc" )

SET ( DAGMCIFACE_SRC ${DAGMCIFACE_SRC} ${DAGMCIFACE_GEANT_SRC} )
INCLUDE_DIRECTORIES( "${DAG_GEANT_SRCDIR}" )

# Geant4 include dirs
INCLUDE_DIRECTORIES( ${INCLUDE_DIRECTORIES} ${GEANT4_DIR}/include/Geant4)
