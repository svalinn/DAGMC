####################################################
#
# CMakeLists.txt for top-level DAGMC repository
#
#####################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT (DAGMC)
SET ( DAGMC_VERSION_MAJOR 1 )
SET ( DAGMC_VERSION_MINOR 0 )

FIND_PACKAGE ( MOAB REQUIRED ) 

###################################
#
# DAGMC PHYSICS INTERFACE LIBRARY
#
###################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MCNP5
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SET ( DAG_MCNP_SRCDIR "mcnp5/src" )
SET ( DAGMCIFACE_MCNP_SRC "${DAG_MCNP_SRCDIR}/mcnp_funcs.cpp" )

SET ( DAGMCIFACE_SRC ${DAGMCIFACE_SRC} ${DAGMCIFACE_MCNP_SRC} )
INCLUDE_DIRECTORIES( "${DAG_MCNP_SRCDIR}")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FLUKA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SET ( DAG_FLUKA_SRCDIR "fluka/src" )
SET ( DAGMCIFACE_FLUKA_SRC "${DAG_FLUKA_SRCDIR}/fluka_funcs.cpp" )

SET ( DAGMCIFACE_SRC ${DAGMCIFACE_SRC} ${DAGMCIFACE_FLUKA_SRC} )
INCLUDE_DIRECTORIES( "${DAG_FLUKA_SRCDIR}" )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GEANT4
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SET ( DAG_GEANT_SRCDIR "geant4/src" )
SET ( DAGMCIFACE_GEANT_SRC "${DAG_GEANT_SRCDIR}/DagSolid.cc" )

SET ( DAGMCIFACE_SRC ${DAGMCIFACE_SRC} ${DAGMCIFACE_GEANT_SRC} )
INCLUDE_DIRETORIES( "${DAG_GEANT_SRCDIR}" )

ADD_LIBRARY(dagmciface ${DAGMCIFACE_SRC} )

###################################
#
# DAGMC TALLY LIBRARY
#
###################################

ADD_SUBDIRECTORY(tally)
