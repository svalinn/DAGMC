set(DRIVERS overlap_check_test.cpp)

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src/overlap_check)

dagmc_install_test(overlap_check_test cpp)

dagmc_install_test_file(overlap.h5m)
