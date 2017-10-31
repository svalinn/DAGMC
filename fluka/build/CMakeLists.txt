set(SRC_FILES mainFluDAG.cpp)

set(LINK_LIBS fludag)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${CMAKE_SOURCE_DIR}/fluka)

dagmc_install_exe(mainfludag)
