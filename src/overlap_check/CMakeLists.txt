message("")

set(LINK_LIBS dagmc)
set(LINK_LIBS_EXTERN_NAMES)

add_subdirectory(build)

if (BUILD_TESTS)
add_subdirectory(test)
endif()
