message(STATUS "Building libdagmc unit tests")

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS pthread)
if (BUILD_STATIC)
  list(APPEND LINK_LIBS gtest-static)
  list(APPEND LINK_LIBS dagmc-static)
  list(APPEND LINK_LIBS pyne_dagmc-static)
else ()
  list(APPEND LINK_LIBS gtest)
  list(APPEND LINK_LIBS dagmc)
  list(APPEND LINK_LIBS pyne_dagmc)
endif (BUILD_STATIC)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Include directories
include_directories(..)
include_directories(${CMAKE_SOURCE_DIR}/dagmc)
include_directories(${GTEST_INCLUDE_DIR})

# Setup unit tests
set(DRIVER dagmciface_unit_test_driver.cc)
setup_test(dagmciface_unit_tests cpp ${DRIVER} "${LINK_LIBS}")

# Install h5m files
install(FILES test_dagmciface.h5m DESTINATION tests)
install(FILES test_dagmciface_impl.h5m DESTINATION tests)

# Enable unit tests
enable_testing()
