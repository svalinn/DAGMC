# CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
MESSAGE ( STATUS "Including Amalgamated PyNE")

SET ( CODE_SRCDIR      "pyne" )
SET ( CODE_SRC_FILES   "pyne.cpp"  )
SET ( CODE_PUB_HEADERS "pyne.h" )

# make pyne library
IF( STATIC_LIB )
  add_library( pyne_dagmc STATIC ${CODE_SRC_FILES})
ELSEIF( NOT STATIC_LIB )
  add_library( pyne_dagmc SHARED ${CODE_SRC_FILES})
ENDIF( STATIC_LIB )
 
target_include_directories( pyne_dagmc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# if we want a static executable
# IF(STATIC_LIB)
#  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#  SET(BUILD_SHARED_LIBRARIES OFF)
#  SET(CMAKE_EXE_LINKER_FLAGS "-static")
# ENDIF(STATIC_LIB)

#target_link_libraries(pyne ${HDF5_hdf5_hl_LIBRARY} ${HDF5_hdf5_LIBRARY})

SET_TARGET_PROPERTIES( pyne_dagmc PROPERTIES PUBLIC_HEADER "${CODE_PUB_HEADERS}" )

INSTALL( TARGETS       pyne_dagmc
  LIBRARY        DESTINATION   "${INSTALL_LIB_DIR}"
  ARCHIVE        DESTINATION   "${INSTALL_LIB_DIR}"
  PUBLIC_HEADER  DESTINATION   "${INSTALL_INCLUDE_DIR}"
  )

target_link_libraries(pyne_dagmc ${HDF5_hdf5_hl_LIBRARY} ${HDF5_hdf5_LIBRARY})




