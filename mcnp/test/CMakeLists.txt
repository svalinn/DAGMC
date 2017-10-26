message(STATUS "Building DAG-MCNP unit tests")

set(LINK_LIBS_SHARED gtest dagmc)
set(LINK_LIBS_STATIC gtest dagmc-static)

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_EXE ${LINK_LIBS_STATIC})
else ()
  set(LINK_LIBS_EXE ${LINK_LIBS_SHARED})
endif ()

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/dagmc)
include_directories(${CMAKE_SOURCE_DIR}/mcnp)

set(DRIVER mcnp5_unit_test_driver.cc ${CMAKE_SOURCE_DIR}/mcnp/mcnp_funcs.cpp)
dagmc_setup_test(dagmcnp_unit_tests cpp "${DRIVER}" "${LINK_LIBS_EXE}")

install(FILES test_geom_legacy.h5m      DESTINATION tests)
install(FILES test_geom_legacy_comp.h5m DESTINATION tests)
install(FILES test_reflecting.h5m       DESTINATION tests)
