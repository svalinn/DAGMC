message("")

file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.hpp")

# Configure version header
configure_file(DagMCVersion.hpp.in DagMCVersion.hpp)
list(APPEND DAGMC_PUB_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/DagMCVersion.hpp)

if (BUILD_STATIC_EXE)
  set(LINK_LIBS lapack blas gfortran quadmath)
else ()
  set(LINK_LIBS)
endif ()
set(LINK_LIBS_EXTERN_NAMES MOAB_LIBRARIES HDF5_LIBRARIES)

include_directories(${CMAKE_BINARY_DIR}/src/dagmc)

dagmc_install_library(dagmc)

add_subdirectory(tools)

if (BUILD_TESTS)
  add_subdirectory(tests)
endif ()
