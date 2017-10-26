message(STATUS "Building FluDAG")

set(SRC_FILES mainFluDAG.cpp)

set(LINK_LIBS_SHARED fludag)
set(LINK_LIBS_STATIC fludag-static)

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_EXE ${LINK_LIBS_STATIC})
else ()
  set(LINK_LIBS_EXE ${LINK_LIBS_SHARED})
endif ()

include_directories(${CMAKE_SOURCE_DIR}/fluka)

add_executable(mainfludag ${SRC_FILES})
target_link_libraries(mainfludag ${LINK_LIBS_EXE})
install(TARGETS mainfludag DESTINATION bin)
