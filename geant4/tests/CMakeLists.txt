set(DRIVERS dagsolid_unit_test_driver.cpp)

set(LINK_LIBS dagsolid)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${CMAKE_SOURCE_DIR}/geant4)
include_directories(${GTEST_INCLUDE_DIR})

dagmc_install_test(dagsolid_unit_tests cpp)

install(FILES test_geom.h5m DESTINATION tests)
