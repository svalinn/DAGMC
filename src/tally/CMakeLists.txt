message("")

file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.hpp")

set(LINK_LIBS lapack dagmc)
set(LINK_LIBS_EXTERN_NAMES)

dagmc_install_library(dagtally)

if (BUILD_TESTS)
  add_subdirectory(tests)
endif ()
