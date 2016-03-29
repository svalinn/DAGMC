message(STATUS "Building FluDAG unit tests")

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS gtest)
list(APPEND LINK_LIBS pthread)
list(APPEND LINK_LIBS dagmciface)
list(APPEND LINK_LIBS ${HDF5_hdf5_hl_LIBRARY})
list(APPEND LINK_LIBS ${HDF5_hdf5_LIBRARY})
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS ${FLUKA_DIR}/libflukahp.a)
list(APPEND LINK_LIBS pyne_dagmc)
list(APPEND LINK_LIBS gfortran)

# Include directories
include_directories(..)

# Setup unit tests
set(DRIVER fludag_unit_test_driver.cc)
setup_test(fludag_unit_tests cpp ${DRIVER} "${LINK_LIBS}")

# Install h5m files
install(FILES slabs.h5m DESTINATION ${INSTALL_SHARE_DIR}/tests)

# Enable unit tests
enable_testing()
