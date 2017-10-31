message(STATUS "Building dagtally library")

set(SRC_FILES)
file(GLOB SRC_FILES ${SRC_FILES} "*.cpp")

set(PUB_HEADERS)
file(GLOB PUB_HEADERS ${PUB_HEADERS} "*.hpp")

set(LINK_LIBS_SHARED dagmc)
set(LINK_LIBS_STATIC dagmc-static)

dagmc_install_library(dagtally)

add_subdirectory(tests)
