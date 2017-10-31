message(STATUS "Building DagSolid unit tests")

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_TEST dagsolid-static)
else ()
  set(LINK_LIBS_TEST dagsolid)
endif ()

include_directories(${CMAKE_SOURCE_DIR}/geant4)
include_directories(${GTEST_INCLUDE_DIR})

set(DRIVERS dagsolid_unit_test_driver.cpp)
dagmc_setup_test(dagsolid_unit_tests cpp)

install(FILES test_geom.h5m DESTINATION tests)
