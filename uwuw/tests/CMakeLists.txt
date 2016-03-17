message(STATUS "Building UWUW unit tests")

# Includes
include_directories(..)

# Libraries
set(UWUW_TEST_LIBS)
list(APPEND UWUW_TEST_LIBS gtest)
list(APPEND UWUW_TEST_LIBS dagmciface)
list(APPEND UWUW_TEST_LIBS pyne_dagmc)
list(APPEND UWUW_TEST_LIBS pthread)
list(APPEND UWUW_TEST_LIBS ${DAG_LINK_LIBS})

# Setup unit tests
set(DRIVER uwuw_unit_test_driver.cpp)
setup_test(uwuw_unit_tests       cpp ${DRIVER} "${UWUW_TEST_LIBS}")
setup_test(uwuw_unit_tests_tally cpp ${DRIVER} "${UWUW_TEST_LIBS}")

# Install h5m files
install(FILES test_uwuw.h5m DESTINATION ${INSTALL_SHARE_DIR}/tests)
install(FILES mat_lib.h5    DESTINATION ${INSTALL_SHARE_DIR}/tests)

# enable unit tests
enable_testing()
