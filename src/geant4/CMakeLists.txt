message("")

find_package(Geant4 10.4.2 REQUIRED)

if ( Geant4_VERSION GREATER_EQUAL 11.0 )
  message(STATUS "Configuring DAGMC for Geant4 version >= 11.0")
  add_compile_definitions(GEANT4_GT_10_6=1)
  add_compile_definitions(GEANT4_GT_11=1)
  elseif( Geant4_VERSION GREATER_EQUAL 10.6 )
  message(STATUS "Configuring DAGMC for Geant4 version >= 10.6")
  add_compile_definitions(GEANT4_GT_10_6=1)
endif()

set(SRC_FILES DagSolid.cc)
set(PUB_HEADERS DagSolid.hh)

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
