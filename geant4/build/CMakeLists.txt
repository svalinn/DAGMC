message(STATUS "Building DAG-Geant4")

set(SRC_FILES)
file(GLOB SRC_FILES ${SRC_FILES} "*.cc")
file(GLOB SRC_FILES ${SRC_FILES} "src/*.cc")

set(LINK_LIBS_SHARED dagsolid)
set(LINK_LIBS_STATIC dagsolid-static)

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_EXE ${LINK_LIBS_STATIC})
else ()
  set(LINK_LIBS_EXE ${LINK_LIBS_SHARED})
endif ()

include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/uwuw)

add_executable(DagGeant4 ${SRC_FILES})
target_link_libraries(DagGeant4 ${LINK_LIBS_EXE})
install(TARGETS DagGeant4 DESTINATION bin)
