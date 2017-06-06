message(STATUS "Building make_watertight tests")

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS gtest)
if (BUILD_STATIC)
  list(APPEND LINK_LIBS makeWatertight-static)
else ()
  list(APPEND LINK_LIBS makeWatertight)
endif (BUILD_STATIC)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Include directories
include_directories(..)
include_directories(${GTEST_INCLUDE_DIR})
include_directories(${MOAB_INCLUDE_DIRS})

# Setup unit tests
set(DRIVER "mw_unit_test_driver.cc test_classes.cpp")
setup_test(make_watertight_cylinder_tests        cpp ${DRIVER} "${LINK_LIBS}")
setup_test(make_watertight_cone_tests            cpp ${DRIVER} "${LINK_LIBS}")
setup_test(make_watertight_no_curve_sphere_tests cpp ${DRIVER} "${LINK_LIBS}")
setup_test(make_watertight_sphere_n_box_test     cpp ${DRIVER} "${LINK_LIBS}")
if(BUILD_CI_TESTS)
  setup_test(make_watertight_regression_tests      cpp ${DRIVER} "${LINK_LIBS}")
endif()

# Install h5m files
install(FILES cones.h5m           DESTINATION tests)
install(FILES cyl.h5m             DESTINATION tests)
install(FILES no_curve_sphere.h5m DESTINATION tests)
install(FILES sphere_n_box.h5m    DESTINATION tests)

# enable tests
enable_testing()
