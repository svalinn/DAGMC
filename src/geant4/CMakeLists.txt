message("")

find_package(Geant4 REQUIRED)

file(GLOB SRC_FILES "*.cc")
file(GLOB PUB_HEADERS "*.hh")

set(LINK_LIBS dagmc pyne_dagmc uwuw)
set(LINK_LIBS_EXTERN_NAMES Geant4_LIBRARIES)

include_directories(${CMAKE_SOURCE_DIR}/src/dagmc)
include_directories(${CMAKE_BINARY_DIR}/src/dagmc)
include_directories(${CMAKE_SOURCE_DIR}/src/pyne)

dagmc_install_library(dagsolid)

add_subdirectory(app)
if (BUILD_TESTS)
  add_subdirectory(tests)
endif ()
