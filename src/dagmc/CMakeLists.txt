message("")

file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.hpp")

# Configure version header
configure_file(DagMCVersion.hpp.in DagMCVersion.hpp)
list(APPEND PUB_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/DagMCVersion.hpp)

set(LINK_LIBS)
set(LINK_LIBS_EXTERN_NAMES MOAB_LIBRARIES HDF5_LIBRARIES)

if (DOUBLE_DOWN)
  list(APPEND LINK_LIBS dd)
endif()

include_directories(${CMAKE_BINARY_DIR}/src/dagmc)

dagmc_install_library(dagmc)

add_subdirectory(tools)

if (BUILD_TESTS)
  add_subdirectory(tests)
endif ()
