message(STATUS "Building DAGMC unit tests")

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_TEST dagmc-static)
else ()
  set(LINK_LIBS_TEST dagmc)
endif ()

include_directories(${GTEST_INCLUDE_DIR})

set(DRIVERS dagmc_unit_test_driver.cc)
dagmc_setup_test(dagmc_unit_tests      cpp)
dagmc_setup_test(dagmc_pointinvol_test cpp)
dagmc_setup_test(dagmc_rayfire_test    cpp)
dagmc_setup_test(dagmc_simple_test     cpp)

install(FILES test_dagmc.h5m      DESTINATION tests)
install(FILES test_dagmc_impl.h5m DESTINATION tests)
install(FILES test_geom.h5m       DESTINATION tests)
