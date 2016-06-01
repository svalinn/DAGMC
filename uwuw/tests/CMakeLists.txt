message(STATUS "Building UWUW unit tests")

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS gtest)
list(APPEND LINK_LIBS pthread)
list(APPEND LINK_LIBS dagmciface)
list(APPEND LINK_LIBS pyne_dagmc)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Include directories
include_directories(..)
include_directories(${GTEST_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})

# Setup unit tests
set(DRIVER uwuw_unit_test_driver.cpp)
setup_test(uwuw_unit_tests       cpp ${DRIVER} "${LINK_LIBS}")
setup_test(uwuw_unit_tests_tally cpp ${DRIVER} "${LINK_LIBS}")

# Install h5m files
install(FILES test_uwuw.h5m DESTINATION tests)
install(FILES mat_lib.h5    DESTINATION tests)

# Enable unit tests
enable_testing()
