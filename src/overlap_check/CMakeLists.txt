message("")

file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.hpp")

set(LINK_LIBS dagmc)
set(LINK_LIBS_EXTERN_NAMES)

set(INSTALL_INCLUDE_DIR_SAVE ${INSTALL_INCLUDE_DIR})
set(INSTALL_INCLUDE_DIR ${INSTALL_INCLUDE_DIR}/overlap_check)

set(INSTALL_INCLUDE_DIR ${INSTALL_INCLUDE_DIR_SAVE})

set(SRC_FILES overlap_check.cpp overlap.cpp)
dagmc_install_exe(overlap_check)
