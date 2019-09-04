message("")

find_package(Eigen3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIRS})

file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.hpp")

set(LINK_LIBS dagmc lapack blas gfortran quadmath)
set(LINK_LIBS_EXTERN_NAMES)

dagmc_install_library(dagtally)

if (BUILD_TESTS)
  add_subdirectory(tests)
endif ()
