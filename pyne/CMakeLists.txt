message(STATUS "Building pyne_dagmc library")

set(CMAKE_CXX_FLAGS_RELEASE        "-O0"   )
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0 -g")

set(SRC_FILES pyne.cpp)

set(PUB_HEADERS pyne.h)

set(LINK_LIBS_SHARED ${HDF5_LIBRARIES_SHARED})
set(LINK_LIBS_STATIC ${HDF5_LIBRARIES_STATIC})

add_library(pyne_dagmc        SHARED ${SRC_FILES})
add_library(pyne_dagmc-static STATIC ${SRC_FILES})
set_target_properties(pyne_dagmc        PROPERTIES PUBLIC_HEADER "${PUB_HEADERS}")
set_target_properties(pyne_dagmc-static PROPERTIES OUTPUT_NAME pyne_dagmc)
target_link_libraries(pyne_dagmc        ${LINK_LIBS_SHARED})
target_link_libraries(pyne_dagmc-static ${LINK_LIBS_STATIC})
install(TARGETS pyne_dagmc        LIBRARY       DESTINATION lib
                                  PUBLIC_HEADER DESTINATION include)
install(TARGETS pyne_dagmc-static ARCHIVE       DESTINATION lib)
