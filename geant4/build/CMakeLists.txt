# $Id: CMakeLists.txt 68058 2013-03-13 14:47:43Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
# cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(DagGeant4)

# set the compiler flags
SET(CMAKE_CXX_FLAGS "-std=c++11")

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

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
# file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# in case we chose static libraries executables
IF(STATIC_LIB)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static")
ENDIF(STATIC_LIB)


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(DagGeant4 exampleN01.cc ${sources} ${headers})
target_link_libraries(DagGeant4 ${Geant4_LIBRARIES}
                                ${HDF5_hdf5_hl_LIBRARY} ${HDF5_hdf5_LIBRARY}
#                                ${PYNE_LIBS_DIR}/libpyne.so
				pyne_dagmc
			        ${DAG_LINK_LIBS}
				dagsolid
				iMesh
				dagmciface)



INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/uwuw)
#INCLUDE_DIRECTORIES("/data/opt/DAGMC/geant4")

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(B1 DEPENDS DagGeant4)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS DagGeant4 DESTINATION bin)


