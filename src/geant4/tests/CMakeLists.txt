set(DRIVERS dagsolid_unit_test_driver.cpp)

set(LINK_LIBS dagsolid)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src/geant4)

dagmc_install_test(dagsolid_unit_tests cpp)

install(FILES test_geom.h5m DESTINATION ${INSTALL_TESTS_DIR})
