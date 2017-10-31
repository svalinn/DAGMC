file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.hpp")

set(LINK_LIBS lapack blas gfortran quadmath)
set(LINK_LIBS_EXTERN_NAMES MOAB_LIBRARIES)

dagmc_install_library(dagmc)

add_subdirectory(tests)
add_subdirectory(tools)
