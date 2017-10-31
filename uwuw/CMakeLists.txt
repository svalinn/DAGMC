file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.hpp")

set(LINK_LIBS_SHARED
  dagmc
  pyne_dagmc
)
set(LINK_LIBS_STATIC
  dagmc-static
  pyne_dagmc-static
)

include_directories(${CMAKE_SOURCE_DIR}/dagmc)

dagmc_install_library(uwuw)

add_subdirectory(build)
add_subdirectory(tests)
