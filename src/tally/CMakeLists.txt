message("")

file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.hpp")

set(LINK_LIBS dagmc lapack blas gfortran quadmath)
set(LINK_LIBS_EXTERN_NAMES)

dagmc_install_library(dagtally)
dagmc_export_target(dagtally)

if (BUILD_TESTS)
  add_subdirectory(tests)
endif ()
