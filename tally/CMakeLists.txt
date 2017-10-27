message(STATUS "Building dagtally library")

set(SRC_FILES)
file(GLOB SRC_FILES ${SRC_FILES} "*.cpp")

set(PUB_HEADERS)
file(GLOB PUB_HEADERS ${PUB_HEADERS} "*.hpp")

set(LINK_LIBS_SHARED dagmc)
set(LINK_LIBS_STATIC dagmc-static)

add_library(dagtally        SHARED ${SRC_FILES})
add_library(dagtally-static STATIC ${SRC_FILES})

set_target_properties(dagtally        PROPERTIES PUBLIC_HEADER "${PUB_HEADERS}")
set_target_properties(dagtally-static PROPERTIES OUTPUT_NAME dagtally)
target_link_libraries(dagtally        ${LINK_LIBS_SHARED})
target_link_libraries(dagtally-static ${LINK_LIBS_STATIC})
install(TARGETS dagtally        LIBRARY       DESTINATION lib
                                PUBLIC_HEADER DESTINATION include)
install(TARGETS dagtally-static ARCHIVE       DESTINATION lib)

add_subdirectory(tests)
