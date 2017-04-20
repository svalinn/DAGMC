#!/bin/bash

# $1: compiler (gcc-4.8, gcc-5, gcc-6, clang-3.8)

set -e

source /root/etc/$1.env

mkdir -p ${build_dir}/hdf5-${hdf5_version}/bld
rm -rf ${install_dir}/hdf5-${hdf5_version}
cd ${build_dir}/hdf5-${hdf5_version}
wget https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-${hdf5_version}/src/hdf5-${hdf5_version}.tar.gz
tar -xzvf hdf5-${hdf5_version}.tar.gz
ln -snf hdf5-${hdf5_version} src
cd bld
../src/configure --enable-shared \
                 --disable-debug \
                 --prefix=${install_dir}/hdf5-${hdf5_version} \
                 CC=${CC} CXX=${CXX} FC=${FC}
make -j`grep -c processor /proc/cpuinfo`
make install
rm -rf ${build_dir}/hdf5-${hdf5_version}
