####################################################
#
# CMakeLists.txt for top-level DAGMC repository
#
#####################################################

CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )
PROJECT ( DAGMC )
SET ( DAGMC_VERSION_MAJOR 1 )
SET ( DAGMC_VERSION_MINOR 0 )

###########################
#
# FIND DEPENDENCIES
#
##########################
# Needed for FindPyne.cmake
SET ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake)

FIND_PACKAGE ( MOAB REQUIRED ) 
MESSAGE ( STATUS "Found MOAB Version: " ${MOAB_VERSION} )
  MESSAGE ( STATUS "MOAB_LIBS_DIR is " ${MOAB_LIBS_DIR})
  MESSAGE ( STATUS "MOAB_INCLUDE_DIR is " ${MOAB_INCLUDE_DIR})
INCLUDE_DIRECTORIES( ${MOAB_INCLUDE_DIRS} )

SET ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} /data/opt/pyne_dev/cmake)

FIND_PACKAGE ( Pyne REQUIRED )
LINK_DIRECTORIES(${PYNE_LIBS_DIR})
INCLUDE_DIRECTORIES(${PYNE_INCLUDE_DIR})

INCLUDE_DIRECTORIES(${HDF5_DIR}/include)
###################################
#
# SETUP GTEST DEPENDENCY
# - only build during make check
#
###################################

SET ( GTEST_INC_DIR "" )
ADD_SUBDIRECTORY( gtest EXCLUDE_FROM_ALL )
INCLUDE_DIRECTORIES( ${GTEST_INC_DIR} )


###################################
#
# Install dirs 
#
###################################

SET( INSTALL_LIB_DIR     lib     CACHE PATH "Installation directory for libraries" )
SET( INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header" )
SET( INSTALL_BIN_DIR     bin     CACHE PATH "Installation directory for executables" )



###################################
#
# DAGMC TALLY LIBRARY
#
###################################

IF ( BUILD_TALLY )
  ADD_SUBDIRECTORY( tally )
ELSEIF( BUILD_TALLY )
  MESSAGE ( STATUS "Not building tagtally interface library")
ENDIF ( BUILD_TALLY )

###################################
#
# DAGMC PHYSICS INTERFACE LIBRARY
#
###################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MCNP5
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IF ( BUILD_MCNP5 OR BUILD_ALL )
  ADD_SUBDIRECTORY( uwuw )
  ADD_SUBDIRECTORY( mcnp5 )
  ADD_SUBDIRECTORY( tally )
ENDIF ( BUILD_MCNP5 OR BUILD_ALL )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FLUKA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IF ( BUILD_FLUKA OR BUILD_ALL )
  ADD_SUBDIRECTORY( uwuw )
  ADD_SUBDIRECTORY(fluka)
ENDIF ( BUILD_FLUKA OR BUILD_ALL )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GEANT4
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IF ( BUILD_GEANT4 OR BUILD_ALL )
  ADD_SUBDIRECTORY( uwuw )
  ADD_SUBDIRECTORY(geant4)
ENDIF ( BUILD_GEANT4 OR BUILD_ALL )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Build the dagmciface lib
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IF ( DI_SRC_FILES )
  ADD_SUBDIRECTORY(common)
  IF ( STATIC_LIB )
   ADD_LIBRARY(dagmciface STATIC ${DI_SRC_FILES} )
  ELSEIF ( NOT STATIC_LIB )
    ADD_LIBRARY(dagmciface SHARED ${DI_SRC_FILES} )
  ENDIF ( STATIC_LIB )

  SET_TARGET_PROPERTIES( dagmciface        PROPERTIES   
                         PUBLIC_HEADER     "${DI_PUB_HEADERS}" 
                         )
                          
  INSTALL( TARGETS       dagmciface 
           LIBRARY        DESTINATION   "${INSTALL_LIB_DIR}"
           ARCHIVE        DESTINATION   "${INSTALL_LIB_DIR}"
           PUBLIC_HEADER  DESTINATION   "${INSTALL_INCLUDE_DIR}"
         )
ELSE ( DI_SRC_FILES )
  MESSAGE ( STATUS "Not building dagmciface physics interface library")
ENDIF ( DI_SRC_FILES )
         





