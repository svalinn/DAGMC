if(OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(LINK_LIBS dagmc)
else()
  set(LINK_LIBS dagmc)
endif()

set(LINK_LIBS_EXTERN_NAMES)

include_directories(${CMAKE_SOURCE_DIR}/src/overlap_check)
set(SRC_FILES overlap_check.cpp ../overlap.cpp ../ProgressBar.cpp)

dagmc_install_exe(overlap_check)
