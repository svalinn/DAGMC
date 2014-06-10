####################################################
#
# CMakeLists.txt for top-level DAGMC repository
#
#####################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT (DAGMC)
SET ( DAGMC_VERSION_MAJOR 1 )
SET ( DAGMC_VERSION_MINOR 0 )

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

FIND_PACKAGE ( MOAB REQUIRED ) 
INCLUDE_DIRECTORIES( ${MOAB_INCLUDE_DIRS} )

###################################
#
# DAGMC PHYSICS INTERFACE LIBRARY
#
###################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MCNP5
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IF ( BUILD_MCNP5 OR BUILD_ALL )
  INCLUDE(mcnp5/CMakeLists.txt)
ENDIF ( BUILD_MCNP5 OR BUILD_ALL )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FLUKA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
<<<<<<< HEAD
IF ( BUILD_FLUKA OR BUILD_ALL )
   MESSAGE ( STATUS "Including FLUKA in dagmciface library.")
   SET ( DAG_FLUKA_SRCDIR "fluka/src" )
   SET ( DAGMCIFACE_FLUKA_SRC "${DAG_FLUKA_SRCDIR}/fluka_funcs.cpp" )
 
   SET ( DAGMCIFACE_SRC ${DAGMCIFACE_SRC} ${DAGMCIFACE_FLUKA_SRC} )
   INCLUDE_DIRECTORIES( "${DAG_FLUKA_SRCDIR}" )
ENDIF ( BUILD_FLUKA OR BUILD_ALL )
=======
#if ( WITH_FLUKA or WITH_ALL )
if ( WITH_FLUKA )
	MESSAGE("BUILDING FLUKA")
	SET ( DAG_FLUKA_SRCDIR "fluka/src" )
	SET ( DAGMCIFACE_FLUKA_SRC "${DAG_FLUKA_SRCDIR}/fluka_funcs.cpp" )
	SET ( DAGMCIFACE_SRC ${DAGMCIFACE_SRC} ${DAGMCIFACE_FLUKA_SRC} )
	INCLUDE_DIRECTORIES( "${DAG_FLUKA_SRCDIR}" )
endif()

>>>>>>> First pass at FLUKA-only build- not good

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GEANT4
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IF ( BUILD_GEANT OR BUILD_ALL )
 INCLUDE(geant4/CMakeLists.txt)
ENDIF ( BUILD_GEANT OR BUILD_ALL )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Build the dagmciface lib
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IF ( DAGMCIFACE_SRC )
  IF ( STATIC_LIB )
   ADD_LIBRARY(dagmciface STATIC ${DAGMCIFACE_SRC} )
  ELSEIF ( NOT STATIC_LIB )
    ADD_LIBRARY(dagmciface SHARED ${DAGMCIFACE_SRC} )
  ENDIF ( STATIC_LIB )
ELSE ( DAGMCIFACE_SRC )
  MESSAGE ( STATUS "Not building dagmciface physics interface library")
ENDIF ( DAGMCIFACE_SRC )



###################################
#
# DAGMC TALLY LIBRARY
#
###################################

IF ( BUILD_TALLY )
   MESSAGE ( STATUS "Building dagmctally library." )
   ADD_SUBDIRECTORY(tally)
ENDIF ( BUILD_TALLY )


