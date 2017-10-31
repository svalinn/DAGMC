file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.hpp")

set(LINK_LIBS_SHARED
  ${MOAB_LIBRARIES_SHARED}
)
set(LINK_LIBS_STATIC
  ${MOAB_LIBRARIES_STATIC}
  # Needed for static exes for some reason
  lapack blas gfortran quadmath
)

dagmc_install_library(dagmc)

add_subdirectory(tests)
add_subdirectory(tools)
