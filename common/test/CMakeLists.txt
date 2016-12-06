message(STATUS "Building dagmciface unit tests")

# Libraries                                                                                                                             
set(LINK_LIBS)
list(APPEND LINK_LIBS gtest)
list(APPEND LINK_LIBS pthread)
list(APPEND LINK_LIBS dagmciface)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Include directories                                                                                                                   
include_directories(..)
include_directories(${CMAKE_SOURCE_DIR}/common)
include_directories(${GTEST_INCLUDE_DIR})

# Setup unit tests                                                                                                                      
set(DRIVER dagmciface_unit_test_driver.cc)
setup_test(dagmciface_unit_tests cpp ${DRIVER} "${LINK_LIBS}")
list(APPEND LINK_LIBS dagmciface)
list(APPEND LINK_LIBS ${HDF5_hdf5_hl_LIBRARY})
list(APPEND LINK_LIBS ${HDF5_hdf5_LIBRARY})
list(APPEND LINK_LIBS ${DAGMC_LIBRARIES})
list(APPEND LINK_LIBS pyne_dagmc)
target_link_libraries(dagmciface_unit_tests ${LINK_LIBS})
# Install h5m files                                                                                                                     
install(FILES test_dagmciface.h5m DESTINATION tests)
install(FILES test_dagmciface_impl.h5m DESTINATION tests)

# Enable unit tests                                                                                                                     
enable_testing()


