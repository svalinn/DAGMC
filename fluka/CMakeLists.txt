message(STATUS "Building fludag library")

set(FLUKA_LIBRARIES ${FLUKA_DIR}/libflukahp.a gfortran)

set(SRC_FILES fluka_funcs.cpp)

set(PUB_HEADERS fluka_funcs.h)

set(LINK_LIBS_SHARED ${FLUKA_LIBRARIES} dagmc)
set(LINK_LIBS_STATIC ${FLUKA_LIBRARIES} dagmc-static)

include_directories(${CMAKE_SOURCE_DIR}/dagmc)

add_library(fludag        SHARED ${SRC_FILES})
add_library(fludag-static STATIC ${SRC_FILES})
set_target_properties(fludag        PROPERTIES PUBLIC_HEADER ${PUB_HEADERS})
set_target_properties(fludag-static PROPERTIES OUTPUT_NAME fludag)
target_link_libraries(fludag        ${LINK_LIBS_SHARED})
target_link_libraries(fludag-static ${LINK_LIBS_STATIC})
install(TARGETS fludag        LIBRARY       DESTINATION lib
                              PUBLIC_HEADER DESTINATION include)
install(TARGETS fludag-static ARCHIVE       DESTINATION lib)

add_subdirectory(build)
add_subdirectory(tests)
