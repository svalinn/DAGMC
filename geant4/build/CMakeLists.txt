message(STATUS "Building DAG-Geant4")

# CXX compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Geant4_CXX_FLAGS}")
if (LIB_TYPE STREQUAL STATIC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--allow-multiple-definition")
endif (LIB_TYPE STREQUAL STATIC)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
			 if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Get Geant4 libraries, include directories, and compile definitions
include(${Geant4_USE_FILE})

# Source files
set(SRC_FILES)
file(GLOB SRC_FILES ${SRC_FILES} "*.cc")
file(GLOB SRC_FILES ${SRC_FILES} "src/*.cc")

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${Geant4_LIBRARIES})
list(APPEND LINK_LIBS ${HDF5_hdf5_hl_LIBRARY}) 
list(APPEND LINK_LIBS ${HDF5_hdf5_LIBRARY})
list(APPEND LINK_LIBS dagmciface)
list(APPEND LINK_LIBS dagsolid)
list(APPEND LINK_LIBS pyne_dagmc)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Include directories
include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/uwuw)
include_directories(${CMAKE_SOURCE_DIR}/common)

# Build DAG-Geant4
add_executable(DagGeant4 ${SRC_FILES})
target_link_libraries(DagGeant4 ${LINK_LIBS})
install(TARGETS DagGeant4 DESTINATION ${INSTALL_BIN_DIR}
                          PERMISSIONS ${PERMS})
