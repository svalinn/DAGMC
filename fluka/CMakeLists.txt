find_package(Fluka REQUIRED)

set(SRC_FILES fluka_funcs.cpp)

set(PUB_HEADERS fluka_funcs.h)

set(LINK_LIBS_SHARED
  ${FLUKA_LIBRARIES}
  dagmc
  pyne_dagmc
  uwuw
)
set(LINK_LIBS_STATIC
  ${FLUKA_LIBRARIES}
  dagmc-static
  pyne_dagmc-static
  uwuw-static
)

include_directories(${CMAKE_SOURCE_DIR}/dagmc)

dagmc_install_library(fludag)

add_subdirectory(build)
add_subdirectory(tests)
