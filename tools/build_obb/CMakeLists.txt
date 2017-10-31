set(SRC_FILES build_obb.cpp)

set(LINK_LIBS_SHARED dagmc)
set(LINK_LIBS_STATIC dagmc-static)

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_EXE ${LINK_LIBS_STATIC})
else ()
  set(LINK_LIBS_EXE ${LINK_LIBS_SHARED})
endif ()

include_directories(${CMAKE_SOURCE_DIR}/dagmc)

add_executable(build_obb ${SRC_FILES})
target_link_libraries(build_obb ${LINK_LIBS_EXE})
install(TARGETS build_obb DESTINATION bin)
