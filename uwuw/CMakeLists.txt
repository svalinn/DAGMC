message(STATUS "Building uwuw_preproc")

set(SRC_FILES uwuw_preproc.cpp)

set(LINK_LIBS_SHARED dagmc)
set(LINK_LIBS_STATIC dagmc-static)

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_EXE ${LINK_LIBS_STATIC})
else ()
  set(LINK_LIBS_EXE ${LINK_LIBS_SHARED})
endif ()

include_directories(${CMAKE_SOURCE_DIR}/dagmc)

add_executable(uwuw_preproc ${SRC_FILES})
target_link_libraries(uwuw_preproc ${LINK_LIBS_EXE})
install(TARGETS uwuw_preproc DESTINATION bin)

add_subdirectory(tests)
