message("")

file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.hpp")

if (BUILD_STATIC_EXE)
  set(LINK_LIBS lapack blas gfortran quadmath)
else ()
  set(LINK_LIBS)
endif ()
set(LINK_LIBS_EXTERN_NAMES MOAB_LIBRARIES HDF5_LIBRARIES)

dagmc_install_library(dagmc)

add_subdirectory(tools)

if (BUILD_TESTS)
  add_subdirectory(tests)
endif ()
