set(DRIVERS overlap_test_driver.cc
            overlap_check_test.cpp
            ${CMAKE_SOURCE_DIR}/src/overlap_check/overlap.cpp
            ${CMAKE_SOURCE_DIR}/src/overlap_check/ProgressBar.cpp)

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src/overlap_check)

dagmc_install_test(overlap_check_test cpp)

dagmc_install_test_file(overlap.h5m)
dagmc_install_test_file(no_overlap.h5m)
dagmc_install_test_file(no_overlap_imp.h5m)
dagmc_install_test_file(enclosed.h5m)
dagmc_install_test_file(small_overlap.h5m)
