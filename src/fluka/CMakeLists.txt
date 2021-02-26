message("")

find_package(Fluka REQUIRED)

file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.h")

set(LINK_LIBS ${FLUKA_LIBRARIES} dagmc pyne_dagmc uwuw)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${CMAKE_SOURCE_DIR}/src/pyne)
include_directories(${CMAKE_SOURCE_DIR}/src/dagmc)
include_directories(${CMAKE_BINARY_DIR}/src/dagmc)
include_directories(${CMAKE_SOURCE_DIR}/src/uwuw)

dagmc_install_library(fludag)

add_subdirectory(app)
if (BUILD_TESTS)
  add_subdirectory(tests)
endif ()
