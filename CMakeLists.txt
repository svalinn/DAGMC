####################################################
#
# CMakeLists.txt for top-level DAGMC repository
#
#####################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT (DAGMC)
SET ( DAGMC_VERSION_MAJOR 1 )
SET ( DAGMC_VERSION_MINOR 0 )

# FIND_PACKAGE ( MOAB REQUIRED ) 

###################################
#
# DAGMC PHYSICS INTERFACE LIBRARY
#
###################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MCNP5
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IF ( BUILD_MCNP5 OR BUILD_ALL )
   MESSAGE ( STATUS "Including MCNP5 in dagmciface library.")
   SET ( DAG_MCNP_SRCDIR "mcnp5/src" )
   SET ( DAGMCIFACE_MCNP_SRC "${DAG_MCNP_SRCDIR}/mcnp_funcs.cpp" )

   SET ( DAGMCIFACE_SRC ${DAGMCIFACE_SRC} ${DAGMCIFACE_MCNP_SRC} )
   INCLUDE_DIRECTORIES( "${DAG_MCNP_SRCDIR}")
ENDIF ( BUILD_MCNP5 OR BUILD_ALL )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FLUKA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IF ( BUILD_FLUKA OR BUILD_ALL )
   MESSAGE ( STATUS "Including FLUKA in dagmciface library.")
   SET ( DAG_FLUKA_SRCDIR "fluka/src" )
   SET ( DAGMCIFACE_FLUKA_SRC "${DAG_FLUKA_SRCDIR}/fluka_funcs.cpp" )
 
   SET ( DAGMCIFACE_SRC ${DAGMCIFACE_SRC} ${DAGMCIFACE_FLUKA_SRC} )
   INCLUDE_DIRECTORIES( "${DAG_FLUKA_SRCDIR}" )
ENDIF ( BUILD_FLUKA OR BUILD_ALL )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GEANT4
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IF ( BUILD_GEANT OR BUILD_ALL )
 INCLUDE(geant4/CMakeLists.txt)
ENDIF ( BUILD_GEANT OR BUILD_ALL )

###################################
#
# DAGMC TALLY LIBRARY
#
###################################

IF ( BUILD_TALLY )
   MESSAGE ( STATUS "Building dagmctally library." )
   ADD_SUBDIRECTORY(tally)
ENDIF ( BUILD_TALLY )


# Build the dagmciface lib
IF ( STATIC_LIB )
  ADD_LIBRARY(dagmciface STATIC ${DAGMCIFACE_SRC} )
ELSEIF ( NOT STATIC_LIB )
  ADD_LIBRARY(dagmciface SHARED ${DAGMCIFACE_SRC} )
ENDIF ( STATIC_LIB )
