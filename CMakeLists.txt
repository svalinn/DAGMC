####################################################
#
# CMakeLists.txt for top-level DAGMC repository
#
#####################################################

CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )
CMAKE_POLICY(SET CMP0004 OLD)

IF(${CMAKE_MAJOR_VERSION} GREATER 3)
  CMAKE_POLICY(SET CMP0040 OLD)
ENDIF(${CMAKE_MAJOR_VERSION} GREATER 3)

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
  MESSAGE ( STATUS "MOAB_LIBRARIES is " ${MOAB_LIBRARIES})
  MESSAGE ( STATUS "MOAB_INCLUDE_DIRS is " ${MOAB_INCLUDE_DIRS})
IF(${MOAB_FOUND})
  INCLUDE_DIRECTORIES( ${MOAB_INCLUDE_DIRS} )
ELSE(${MOAB_FOUND})
  MESSAGE(STATUS "Downloading External MOAB...")
  ADD_SUBDIRECTORY(${CMAKE_CURRENT_LIST_DIR}/cmake/moab)
ENDIF(${MOAB_FOUND})


SET ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH})

FIND_PACKAGE ( HDF5 REQUIRED )
IF(${HDF5_FOUND})
  INCLUDE_DIRECTORIES( ${HDF5_INCLUDE_DIRS} )
ELSE(${HDF5_FOUND})
  MESSGAGE(STATUS "Downloading External HDF5...")
  ADD_SUBDIRECTORY(${CMAKE_CURRENT_LIST_DIR}/cmake/hdf5)
ENDIF(${HDF5_FOUND})

# INCLUDE_DIRECTORIES(${HDF5_DIR}/include)
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
# SET THE AMALGAMTED PYNE
#
###################################
# we now have PyNE Integrated
ADD_SUBDIRECTORY(pyne)

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
  ADD_SUBDIRECTORY( mcnp5 )
  ADD_SUBDIRECTORY( tally )
ENDIF ( BUILD_MCNP5 OR BUILD_ALL )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FLUKA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IF ( BUILD_FLUKA OR BUILD_ALL )
  ADD_SUBDIRECTORY(fluka)
ENDIF ( BUILD_FLUKA OR BUILD_ALL )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GEANT4
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IF ( BUILD_GEANT4 OR BUILD_ALL )
  ADD_SUBDIRECTORY(geant4)
ENDIF ( BUILD_GEANT4 OR BUILD_ALL )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Build the dagmciface lib
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IF ( BUILD_GEANT4 OR BUILD_MCNP5 OR BUILD_FLUKA OR BUILD_ALL )
  ADD_SUBDIRECTORY( uwuw )
ENDIF ( BUILD_GEANT4 OR BUILD_MCNP5 OR BUILD_FLUKA OR BUILD_ALL )

# if we want static libs & executables
IF(STATIC_LIB)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static")
ENDIF(STATIC_LIB)


IF ( DI_SRC_FILES )
  ADD_SUBDIRECTORY(common)

  IF( STATIC_LIB )
    ADD_LIBRARY(dagmciface STATIC ${DI_SRC_FILES} )
  ELSEIF( NOT STATIC_LIB )
    ADD_LIBRARY(dagmciface SHARED ${DI_SRC_FILES} )
  ENDIF( STATIC_LIB )


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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Add the tools subdirectory
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ADD_SUBDIRECTORY(tools)


