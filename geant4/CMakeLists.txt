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
list(APPEND GEANT4_REQUIRED_COMPONENTS static)
find_package(Geant4 REQUIRED COMPONENTS ${GEANT4_REQUIRED_COMPONENTS})
set(Geant4_LIBRARIES_STATIC ${Geant4_LIBRARIES})

# Get Geant4 libraries, include directories, and compile definitions
include(${Geant4_USE_FILE})

set(SRC_FILES DagSolid.cc)

set(PUB_HEADERS DagSolid.hh)

include_directories(${CMAKE_CURRENT_LIST_DIR})
include_directories(${CMAKE_SOURCE_DIR}/dagmc)

set(LINK_LIBS_SHARED dagmc        ${Geant4_LIBRARIES_SHARED})
set(LINK_LIBS_STATIC dagmc-static ${Geant4_LIBRARIES_STATIC})

add_library(dagsolid        SHARED ${SRC_FILES})
add_library(dagsolid-static STATIC ${SRC_FILES})
set_target_properties(dagsolid        PROPERTIES PUBLIC_HEADER "${PUB_HEADERS}")
set_target_properties(dagsolid-static PROPERTIES OUTPUT_NAME dagsolid)
target_link_libraries(dagsolid        ${LINK_LIBS_SHARED})
target_link_libraries(dagsolid-static ${LINK_LIBS_STATIC})
install(TARGETS dagsolid        LIBRARY       DESTINATION lib
                                PUBLIC_HEADER DESTINATION include)
install(TARGETS dagsolid-static ARCHIVE       DESTINATION lib)

add_subdirectory(build)
add_subdirectory(tests)
