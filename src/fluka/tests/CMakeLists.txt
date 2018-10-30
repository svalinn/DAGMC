set(DRIVERS fludag_unit_test_driver.cc)

set(LINK_LIBS fludag)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src/fluka)

# Note: this test is known to fail
dagmc_install_test(fludag_unit_tests cpp)

dagmc_install_test_file(slabs.h5m)
# We use a file from the MCNP test folder, so if not building MCNP5/6 we have to install it
if (NOT BUILD_MCNP5 AND NOT BUILD_MCNP6)
  install(FILES ${CMAKE_SOURCE_DIR}/src/mcnp/tests/test_geom_legacy.h5m DESTINATION ${INSTALL_TESTS_DIR})
  configure_file(${CMAKE_SOURCE_DIR}/src/mcnp/tests/test_geom_legacy.h5m ${CMAKE_CURRENT_BINARY_DIR}/${filename} COPYONLY)
endif (NOT BUILD_MCNP5 AND NOT BUILD_MCNP6)
