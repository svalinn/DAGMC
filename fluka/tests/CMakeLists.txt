message(STATUS "Building FluDAG unit tests")

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS pthread)
if (BUILD_STATIC)
  list(APPEND LINK_LIBS gtest-static)
  list(APPEND LINK_LIBS fludag-static)
else ()
  list(APPEND LINK_LIBS gtest)
  list(APPEND LINK_LIBS fludag)
endif (BUILD_STATIC)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Include directories
include_directories(..)
include_directories(${GTEST_INCLUDE_DIR})

# Setup unit tests
set(DRIVER fludag_unit_test_driver.cc)
setup_test(fludag_unit_tests cpp ${DRIVER} "${LINK_LIBS}")

# Install h5m files
install(FILES slabs.h5m DESTINATION tests)
# We use a file from the mcnp5 test folder, so if not building mcnp5 we have to
# install it
if (NOT BUILD_MCNP5)
  install(FILES ${CMAKE_SOURCE_DIR}/mcnp/test/test_geom_legacy.h5m DESTINATION tests)
endif (NOT BUILD_MCNP5)

# Enable unit tests
enable_testing()
