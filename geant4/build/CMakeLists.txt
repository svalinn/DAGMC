file(GLOB SRC_FILES "*.cc" "src/*.cc")

set(LINK_LIBS_SHARED dagsolid)
set(LINK_LIBS_STATIC dagsolid-static)

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_EXE ${LINK_LIBS_STATIC})
else ()
  set(LINK_LIBS_EXE ${LINK_LIBS_SHARED})
endif ()

include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/uwuw)

dagmc_install_exe(DagGeant4)
