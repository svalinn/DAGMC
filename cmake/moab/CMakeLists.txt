cmake_minimum_required(VERSION 3.0)
project(moab)
include(ExternalProject)

set(moab_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/moab}")
set(moab_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/moab")
set(moab_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${moab_INSTALL_DIR})

message("moab_PREFIX='${moab_PREFIX}'")
message("moab_INSTALL_DIR='${moab_INSTALL_DIR}'")
message("moab_CMAKE_ARGS='${moab_CMAKE_ARGS}'")

ExternalProject_Add(moab
  PREFIX ${moab_PREFIX}
  GIT_REPOSITORY https://bitbucket.org/fathomteam/moab.git 
  INSTALL_DIR ${moab_INSTALL_DIR}
  CMAKE_ARGS ${moab_CMAKE_ARGS}
)

