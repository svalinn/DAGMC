message(STATUS "Building FluDAG unit tests")

set(LINK_LIBS_SHARED gtest fludag)
set(LINK_LIBS_STATIC gtest fludag-static)

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_EXE ${LINK_LIBS_STATIC})
else ()
  set(LINK_LIBS_EXE ${LINK_LIBS_SHARED})
endif ()

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/fluka)

set(DRIVER fludag_unit_test_driver.cc)
dagmc_setup_test(fludag_unit_tests cpp ${DRIVER} "${LINK_LIBS_EXE}")

install(FILES slabs.h5m DESTINATION tests)
# We use a file from the MCNP test folder, so if not building MCNP5/6 we have to install it
if (NOT BUILD_MCNP5 AND NOT BUILD_MCNP6)
  install(FILES ${CMAKE_SOURCE_DIR}/mcnp/test/test_geom_legacy.h5m DESTINATION tests)
endif (NOT BUILD_MCNP5 AND NOT BUILD_MCNP6)
