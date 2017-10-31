if (BUILD_STATIC_EXE)
  set(LINK_LIBS_TEST
    dagmc-static
    pyne_dagmc-static
    uwuw-static
  )
else ()
  set(LINK_LIBS_TEST
    dagmc
    pyne_dagmc
    uwuw
  )
endif ()

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/dagmc)
include_directories(${CMAKE_SOURCE_DIR}/mcnp)

set(DRIVERS mcnp5_unit_test_driver.cc ${CMAKE_SOURCE_DIR}/mcnp/mcnp_funcs.cpp)
dagmc_setup_test(dagmcnp_unit_tests cpp)

install(FILES test_geom_legacy.h5m      DESTINATION tests)
install(FILES test_geom_legacy_comp.h5m DESTINATION tests)
install(FILES test_reflecting.h5m       DESTINATION tests)
