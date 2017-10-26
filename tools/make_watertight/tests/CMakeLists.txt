message(STATUS "Building make_watertight tests")

set(LINK_LIBS_SHARED gtest makeWatertight)
set(LINK_LIBS_STATIC gtest makeWatertight-static)

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_EXE ${LINK_LIBS_STATIC})
else ()
  set(LINK_LIBS_EXE ${LINK_LIBS_SHARED})
endif ()

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/tools/make_watertight)

set(DRIVER "mw_unit_test_driver.cc test_classes.cpp")
dagmc_setup_test(make_watertight_cylinder_tests        cpp ${DRIVER} "${LINK_LIBS_EXE}")
dagmc_setup_test(make_watertight_cone_tests            cpp ${DRIVER} "${LINK_LIBS_EXE}")
dagmc_setup_test(make_watertight_no_curve_sphere_tests cpp ${DRIVER} "${LINK_LIBS_EXE}")
dagmc_setup_test(make_watertight_sphere_n_box_test     cpp ${DRIVER} "${LINK_LIBS_EXE}")
if (BUILD_CI_TESTS)
  dagmc_setup_test(make_watertight_regression_tests      cpp ${DRIVER} "${LINK_LIBS_EXE}")
endif ()

install(FILES cones.h5m           DESTINATION tests)
install(FILES cyl.h5m             DESTINATION tests)
install(FILES no_curve_sphere.h5m DESTINATION tests)
install(FILES sphere_n_box.h5m    DESTINATION tests)
