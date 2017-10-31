set(LINK_LIBS makeWatertight)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${CMAKE_SOURCE_DIR}/tools/make_watertight)

set(SRC_FILES make_watertight.cpp)
dagmc_install_exe(make_watertight)
set(SRC_FILES check_watertight.cpp)
dagmc_install_exe(check_watertight)
