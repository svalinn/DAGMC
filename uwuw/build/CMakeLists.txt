set(SRC_FILES uwuw_preproc.cpp)

set(LINK_LIBS uwuw)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${CMAKE_SOURCE_DIR}/uwuw)

dagmc_install_exe(uwuw_preproc)
