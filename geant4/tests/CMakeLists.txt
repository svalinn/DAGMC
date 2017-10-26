message(STATUS "Building DagSolid unit tests")

set(LINK_LIBS_SHARED gtest dagsolid)
set(LINK_LIBS_STATIC gtest dagsolid-static)

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_EXE ${LINK_LIBS_STATIC})
else ()
  set(LINK_LIBS_EXE ${LINK_LIBS_SHARED})
endif ()

include_directories(${CMAKE_SOURCE_DIR}/geant4)
include_directories(${GTEST_INCLUDE_DIR})

set(DRIVER dagsolid_unit_test_driver.cpp)
dagmc_setup_test(dagsolid_unit_tests cpp ${DRIVER} "${LINK_LIBS_EXE}")

install(FILES test_geom.h5m DESTINATION tests)
