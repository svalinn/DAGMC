set(SRC_FILES build_obb.cpp)

set(LINK_LIBS dagmc)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${CMAKE_SOURCE_DIR}/src/dagmc)

dagmc_install_exe(build_obb)
