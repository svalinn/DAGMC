project(MCNP5)

cmake_minimum_required(VERSION 2.8)

set(PERMS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE)

set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Source)
set(UTIL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Utilities)

include(${SRC_DIR}/CMakeLists.txt)
include(${UTIL_DIR}/CMakeLists.txt)
