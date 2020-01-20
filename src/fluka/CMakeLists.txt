message("")

find_package(Fluka REQUIRED)

set(SRC_FILES fluka_funcs.cpp)
set(PUB_HEADERS fluka_funcs.h)

set(LINK_LIBS ${FLUKA_LIBRARIES} dagmc pyne_dagmc uwuw)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${CMAKE_SOURCE_DIR}/src/pyne)
include_directories(${CMAKE_SOURCE_DIR}/src/dagmc)
include_directories(${CMAKE_BINARY_DIR}/src/dagmc)

dagmc_install_library(fludag)

add_subdirectory(app)
if (BUILD_TESTS)
  add_subdirectory(tests)
endif ()
