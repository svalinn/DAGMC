message("")

find_package(Geant4 REQUIRED)

set(SRC_FILES DagSolid.cc)
set(PUB_HEADERS DagSolid.hh)

set(LINK_LIBS dagmc pyne_dagmc uwuw)
set(LINK_LIBS_EXTERN_NAMES Geant4_LIBRARIES)

include_directories(${CMAKE_SOURCE_DIR}/src/dagmc)
include_directories(${CMAKE_BINARY_DIR}/src/dagmc)
include_directories(${CMAKE_SOURCE_DIR}/src/pyne)

set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--disable-new-dtags")

dagmc_install_library(dagsolid)

add_subdirectory(app)
if (BUILD_TESTS)
  add_subdirectory(tests)
endif ()
