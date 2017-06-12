message(STATUS "Building the DagSolid Library for Geant4")

# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
set(GEANT4_REQUIRED_COMPONENTS)
if (WITH_GEANT4_UIVIS)
  list(APPEND GEANT4_REQUIRED_COMPONENTS ui_all)
  list(APPEND GEANT4_REQUIRED_COMPONENTS vis_all)
endif (WITH_GEANT4_UIVIS)
find_package(Geant4 REQUIRED COMPONENTS ${GEANT4_REQUIRED_COMPONENTS})
set(Geant4_LIBRARIES_SHARED ${Geant4_LIBRARIES})
if (BUILD_STATIC)
  list(APPEND GEANT4_REQUIRED_COMPONENTS static)
endif (BUILD_STATIC)
find_package(Geant4 REQUIRED COMPONENTS ${GEANT4_REQUIRED_COMPONENTS})
set(Geant4_LIBRARIES_STATIC ${Geant4_LIBRARIES})

# Fix bug in Geant4's CMakeConfig (v10.02 and earlier?)
if (BUILD_STATIC AND ${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
  list(FIND Geant4_LIBRARIES_STATIC -Wl,-force_load LOC)
  if (NOT ${LOC} STREQUAL "-1")
    list(REMOVE_ITEM Geant4_LIBRARIES_STATIC -Wl,-force_load)
    list(INSERT Geant4_LIBRARIES_STATIC ${LOC} -Wl,--whole-archive)
    math(EXPR LOC "${LOC}+2")
    list(INSERT Geant4_LIBRARIES_STATIC ${LOC} -Wl,--no-whole-archive)
  endif (NOT ${LOC} STREQUAL "-1")
endif (BUILD_STATIC AND ${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)

# Get Geant4 libraries, include directories, and compile definitions
include(${Geant4_USE_FILE})

# Source files
set(SRC_FILES DagSolid.cc)

# Public headers
set(PUB_HEADERS DagSolid.hh)

# Include directories
include_directories(.)
include_directories(build/include)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/dagmc)
include_directories(${GEANT4_DIR}/include/Geant4)

# Libraries
set(LINK_LIBS)
set(LINK_LIBS_STATIC)
list(APPEND LINK_LIBS        ${Geant4_LIBRARIES_SHARED})
list(APPEND LINK_LIBS        ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS        dagmc)
list(APPEND LINK_LIBS_STATIC ${Geant4_LIBRARIES_STATIC})
list(APPEND LINK_LIBS_STATIC ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS_STATIC dagmc-static)

# Build DagSolid library
add_library(dagsolid        SHARED ${SRC_FILES})
add_library(dagsolid-static STATIC ${SRC_FILES})
set_target_properties(dagsolid        PROPERTIES PUBLIC_HEADER "${PUB_HEADERS}")
set_target_properties(dagsolid-static PROPERTIES OUTPUT_NAME dagsolid)
target_link_libraries(dagsolid        ${LINK_LIBS})
target_link_libraries(dagsolid-static ${LINK_LIBS_STATIC})
install(TARGETS dagsolid        LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                                PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR})
install(TARGETS dagsolid-static ARCHIVE       DESTINATION ${INSTALL_LIB_DIR})

# Build DAG-Geant4
add_subdirectory(build)

# Build the DAG-Geant4 tests
add_subdirectory(tests)
