####################################################
#
# CMakeLists.txt for top-level DAGMC repository
#
#####################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT (DAGMC)
SET ( DAGMC_VERSION_MAJOR 1 )
SET ( DAGMC_VERSION_MINOR 0 )

# LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

FIND_PACKAGE ( MOAB REQUIRED ) 

MESSAGE ( STATUS "Found MOAB Version: " ${MOAB_VERSION} )

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
IF ( BUILD_FLUKA OR BUILD_ALL )
   INCLUDE(fluka/CMakeLists.txt)
ENDIF ( BUILD_FLUKA OR BUILD_ALL )

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


