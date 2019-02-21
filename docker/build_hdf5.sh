#!/bin/bash

# $1: compiler (gcc, clang)
# $2: hdf5 version (system, 1.10.4)

set -e

export compiler=$1
export hdf5_version=$2

source /root/etc/env.sh

if [ "${hdf5_version:3:1}" == "." ]; then
  hdf5_version_major=${hdf5_version::3}
elif [ "${hdf5_version:4:1}" == "." ]; then
  hdf5_version_major=${hdf5_version::4}
fi

rm -rf ${hdf5_build_dir}/bld ${hdf5_install_dir}
mkdir -p ${hdf5_build_dir}/bld
cd ${hdf5_build_dir}
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${hdf5_version_major}/hdf5-${hdf5_version}/src/hdf5-${hdf5_version}.tar.gz
tar -xzvf hdf5-${hdf5_version}.tar.gz
ln -snf hdf5-${hdf5_version} src
cd bld
../src/configure --enable-shared \
                 --prefix=${hdf5_install_dir} \
                 CC=${CC} CXX=${CXX} FC=${FC}
make -j${jobs}
make install
cd
rm -rf ${hdf5_build_dir}
