message(STATUS "Building fludag library")

find_package(Fluka REQUIRED)

set(SRC_FILES fluka_funcs.cpp)

set(PUB_HEADERS fluka_funcs.h)

set(LINK_LIBS_SHARED ${FLUKA_LIBRARIES} dagmc)
set(LINK_LIBS_STATIC ${FLUKA_LIBRARIES} dagmc-static)

include_directories(${CMAKE_SOURCE_DIR}/dagmc)

dagmc_install_library(fludag)

add_subdirectory(build)
add_subdirectory(tests)
