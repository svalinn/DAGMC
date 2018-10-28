set(DRIVERS dagmc_unit_test_driver.cc)

set(LINK_LIBS dagmc)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src/dagmc)
include_directories(${CMAKE_BINARY_DIR}/src/dagmc)

dagmc_install_test(dagmc_unit_tests      cpp)
dagmc_install_test(dagmc_pointinvol_test cpp)
dagmc_install_test(dagmc_rayfire_test    cpp)
dagmc_install_test(dagmc_simple_test     cpp)

dagmc_install_test_file(test_dagmc.h5m)
dagmc_install_test_file(test_dagmc_impl.h5m)
dagmc_install_test_file(test_geom.h5m)
