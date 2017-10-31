message(STATUS "Building the DagSolid Library for Geant4")

find_package(Geant4 REQUIRED)

set(SRC_FILES DagSolid.cc)

set(PUB_HEADERS DagSolid.hh)

include_directories(${CMAKE_CURRENT_LIST_DIR})
include_directories(${CMAKE_SOURCE_DIR}/dagmc)

set(LINK_LIBS_SHARED dagmc        ${Geant4_LIBRARIES_SHARED})
set(LINK_LIBS_STATIC dagmc-static ${Geant4_LIBRARIES_STATIC})

dagmc_install_library(dagsolid)

add_subdirectory(build)
add_subdirectory(tests)
