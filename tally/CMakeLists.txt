file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.hpp")

set(LINK_LIBS_SHARED dagmc)
set(LINK_LIBS_STATIC dagmc-static)

dagmc_install_library(dagtally)

add_subdirectory(tests)
