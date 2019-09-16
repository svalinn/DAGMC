message("")

set(LINK_LIBS dagmc)
set(LINK_LIBS_EXTERN_NAMES)

set(INSTALL_INCLUDE_DIR_SAVE ${INSTALL_INCLUDE_DIR})
set(INSTALL_INCLUDE_DIR ${INSTALL_INCLUDE_DIR}/overlap_check)

set(INSTALL_INCLUDE_DIR ${INSTALL_INCLUDE_DIR_SAVE})

include_directories(${CMAKE_SOURCE_DIR}/src/overlap_check)
set(SRC_FILES overlap_check.cpp ../overlap.cpp ../ProgressBar.cpp)

dagmc_install_exe(overlap_check)
