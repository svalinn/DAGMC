message("")

file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.hpp")

set(LINK_LIBS dagmc pyne_dagmc)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${CMAKE_SOURCE_DIR}/src/dagmc)
include_directories(${CMAKE_BINARY_DIR}/src/dagmc)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../pyne)

dagmc_install_library(uwuw)

add_subdirectory(app)
if (BUILD_TESTS)
  add_subdirectory(tests)
endif ()
