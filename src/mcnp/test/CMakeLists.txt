message(STATUS "Building DAG-MCNP unit tests")

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS pthread)
list(APPEND LINK_LIBS gtest)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS mcnpfuncs)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/src/mcnp)
include_directories(${CMAKE_SOURCE_DIR}/src/dagmc)
include_directories(${MOAB_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${GTEST_INCLUDE_DIR})

# Setup unit tests
set(DRIVER mcnp5_unit_test_driver.cc)
setup_test(dagmcnp_unit_tests cpp ${DRIVER} "${LINK_LIBS}")

# Install h5m files
install(FILES test_geom_legacy.h5m      DESTINATION tests)
install(FILES test_geom_legacy_comp.h5m DESTINATION tests)
install(FILES test_reflecting.h5m DESTINATION tests)

# Enable unit tests
enable_testing()
