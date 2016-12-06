message(STATUS "Building FluDAG unit tests")

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS gtest)
list(APPEND LINK_LIBS pthread)
list(APPEND LINK_LIBS fludag)
list(APPEND LINK_LIBS dagmciface)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Include directories
include_directories(..)
include_directories(${GTEST_INCLUDE_DIR})

# Setup unit tests
set(DRIVER fludag_unit_test_driver.cc)
setup_test(fludag_unit_tests cpp ${DRIVER} "${LINK_LIBS}")

# Install h5m files
install(FILES slabs.h5m DESTINATION tests)
# we use a file from the mcnp5 test folder so if not buildiding
# mcnp5 we have to install it
if(BUILD_MCNP5)
# do nothing
else(BUILD_MCNP5)
  install(FILES ${CMAKE_SOURCE_DIR}/mcnp/test/test_geom_legacy.h5m DESTINATION tests)
endif(BUILD_MCNP5)

# Enable unit tests
enable_testing()
