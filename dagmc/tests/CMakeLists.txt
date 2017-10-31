set(DRIVERS dagmc_unit_test_driver.cc)

set(LINK_LIBS dagmc)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/dagmc)

dagmc_install_test(dagmc_unit_tests      cpp)
dagmc_install_test(dagmc_pointinvol_test cpp)
dagmc_install_test(dagmc_rayfire_test    cpp)
dagmc_install_test(dagmc_simple_test     cpp)

install(FILES test_dagmc.h5m      DESTINATION tests)
install(FILES test_dagmc_impl.h5m DESTINATION tests)
install(FILES test_geom.h5m       DESTINATION tests)
