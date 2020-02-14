message("")

set(SRC_FILES pyne.cpp)
set(PUB_HEADERS pyne.h)

set(LINK_LIBS)
set(LINK_LIBS_EXTERN_NAMES HDF5_LIBRARIES MOAB_LIBRARIES)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE        "-O0")
  set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0 -g")
endif ()

dagmc_install_library(pyne_dagmc)
dagmc_export_target(pyne_dagmc)
