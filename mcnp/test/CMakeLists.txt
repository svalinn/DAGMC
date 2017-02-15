message(STATUS "Building DAG-MCNP5 unit tests")

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS gtest)
list(APPEND LINK_LIBS pthread)
list(APPEND LINK_LIBS dagmciface)
list(APPEND LINK_LIBS dagtally)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Include directories
include_directories(..)
include_directories(${CMAKE_SOURCE_DIR}/common)
include_directories(${MOAB_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${GTEST_INCLUDE_DIR})

# Setup unit tests
set(TEST_NAME dagmcnp_unit_tests)
set(DRIVER mcnp5_unit_test_driver.cc)
add_executable(${TEST_NAME} ${TEST_NAME}.cpp ../mcnp_funcs.cpp ${DRIVER})
target_link_libraries(${TEST_NAME} "${LINK_LIBS}")
install(TARGETS ${TEST_NAME} DESTINATION tests
                             PERMISSIONS ${PERMS})
add_test(${TEST_NAME} ${TEST_NAME})

# Install h5m files
install(FILES test_geom_legacy.h5m DESTINATION tests)
install(FILES test_geom_legacy_comp.h5m DESTINATION tests)

# Enable unit tests
enable_testing()
