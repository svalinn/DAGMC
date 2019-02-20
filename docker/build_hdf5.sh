#!/bin/bash

# $1: compiler (gcc, clang)

set -e

source /root/etc/env.sh $1

mkdir -p ${build_dir}/hdf5-${hdf5_version}/bld
rm -rf ${install_dir}/hdf5-${hdf5_version}
cd ${build_dir}/hdf5-${hdf5_version}
if   [ "${hdf5_version:3:1}" == "." ]; then hdf5_version_major=${hdf5_version::3}
elif [ "${hdf5_version:4:1}" == "." ]; then hdf5_version_major=${hdf5_version::4}
fi
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${hdf5_version_major}/hdf5-${hdf5_version}/src/hdf5-${hdf5_version}.tar.gz
tar -xzvf hdf5-${hdf5_version}.tar.gz
ln -snf hdf5-${hdf5_version} src
cd bld
../src/configure --enable-shared \
                 --prefix=${install_dir}/hdf5-${hdf5_version} \
                 CC=${CC} CXX=${CXX} FC=${FC}
make -j${jobs}
make install
cd
rm -rf ${build_dir}/hdf5-${hdf5_version}
