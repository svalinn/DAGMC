set(DRIVERS mw_unit_test_driver.cc test_classes.cpp)

set(LINK_LIBS makeWatertight)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src/make_watertight)

dagmc_install_test(make_watertight_cylinder_tests        cpp)
dagmc_install_test(make_watertight_cone_tests            cpp)
dagmc_install_test(make_watertight_no_curve_sphere_tests cpp)
dagmc_install_test(make_watertight_sphere_n_box_test     cpp)
if (BUILD_CI_TESTS)
  dagmc_install_test(make_watertight_regression_tests      cpp)
endif ()

install(FILES cones.h5m           DESTINATION ${INSTALL_TESTS_DIR})
install(FILES cyl.h5m             DESTINATION ${INSTALL_TESTS_DIR})
install(FILES no_curve_sphere.h5m DESTINATION ${INSTALL_TESTS_DIR})
install(FILES sphere_n_box.h5m    DESTINATION ${INSTALL_TESTS_DIR})
