set(DRIVERS mcnp5_unit_test_driver.cc ${CMAKE_SOURCE_DIR}/mcnp/mcnp_funcs.cpp)

set(LINK_LIBS dagmc pyne_dagmc uwuw)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/dagmc)
include_directories(${CMAKE_SOURCE_DIR}/mcnp)

dagmc_install_test(dagmcnp_unit_tests cpp)

install(FILES test_geom_legacy.h5m      DESTINATION ${INSTALL_TESTS_DIR})
install(FILES test_geom_legacy_comp.h5m DESTINATION ${INSTALL_TESTS_DIR})
install(FILES test_reflecting.h5m       DESTINATION ${INSTALL_TESTS_DIR})
