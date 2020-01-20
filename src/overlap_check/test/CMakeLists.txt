if(OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(LINK_LIBS dagmc)

  if(BUILD_STATIC_EXE)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
  else()
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()

else()
  set(LINK_LIBS dagmc)
endif()

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
