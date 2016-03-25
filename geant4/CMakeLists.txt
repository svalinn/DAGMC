message(STATUS "Building the DagSolid Library for Geant4")

# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Setup Geant4 include directories and compile definitions
include(${Geant4_USE_FILE})

# Source files
set(SRC_FILES DagSolid.cc)

# Public headers
set(PUB_HEADERS DagSolid.hh)

# Include directories
include_directories(.)
include_directories(./build/include)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})
include_directories(${GEANT4_DIR}/include/Geant4)

# Build DagSolid library
add_library(dagsolid ${LIB_TYPE} ${SRC_FILES})
set_target_properties(dagsolid PROPERTIES PUBLIC_HEADER ${PUB_HEADERS})
install(TARGETS dagsolid LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                         ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
                         PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}
                         LIBRARY       PERMISSIONS ${PERMS})

# Build DAG-Geant4
add_subdirectory(build)

# Build the DAG-Geant4 tests
add_subdirectory(tests)
