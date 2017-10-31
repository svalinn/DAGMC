find_package(Geant4 REQUIRED)

set(SRC_FILES DagSolid.cc)
 
set(PUB_HEADERS DagSolid.hh)

set(LINK_LIBS_SHARED
  ${Geant4_LIBRARIES_SHARED}
  dagmc
  pyne_dagmc
  uwuw
)
set(LINK_LIBS_STATIC
  ${Geant4_LIBRARIES_STATIC}
  dagmc-static
  pyne_dagmc-static
  uwuw-static
)

include_directories(${CMAKE_CURRENT_LIST_DIR})
include_directories(${CMAKE_SOURCE_DIR}/dagmc)

dagmc_install_library(dagsolid)

add_subdirectory(build)
add_subdirectory(tests)
