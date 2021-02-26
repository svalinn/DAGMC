message("")

file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.hpp")

set(LINK_LIBS)
set(LINK_LIBS_EXTERN_NAMES HDF5_LIBRARIES MOAB_LIBRARIES)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE        "-O0")
  set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0 -g")
endif ()

if(NOT (${HDF5_VERSION} VERSION_LESS 1.12.0))
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DH5Oget_info_by_idx_vers=1 -DH5O_info_t_vers=1 -DH5Oget_info_by_name_vers=1")
 endif()

dagmc_install_library(pyne_dagmc)
