cmake_minimum_required(VERSION 3.0)
project(hdf5)
include(ExternalProject)

set(hdf5_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/hdf5")
set(hdf5_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/hdf5")
set(hdf5_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${hdf5_INSTALL_DIR}
                    -DHDF5_BUILD_CPP_LIB=ON
		    -BUILD_SHARED_LIBS=ON
		    -BUILD_STATIC_EXECS=ON)

message("hdf5_PREFIX='${hdf5_PREFIX}'")
message("hdf5_INSTALL_DIR='${moab_INSTALL_DIR}'")
message("hdf5_CMAKE_ARGS='${hdf5_CMAKE_ARGS}'")

ExternalProject_Add(hdf5
  URL http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.13/src/hdf5-1.8.13.tar.gz
  PREFIX ${hdf5_PREFIX}
  INSTALL_DIR ${hdf5_INSTALL_DIR}
  CMAKE_ARGS ${hdf5_CMAKE_ARGS}
)

