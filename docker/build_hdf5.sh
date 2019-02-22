#!/bin/bash

set -e

source /root/etc/env.sh

if [ "${HDF5_VERSION:3:1}" == "." ]; then
  HDF5_VERSION_major=${HDF5_VERSION::3}
elif [ "${HDF5_VERSION:4:1}" == "." ]; then
  HDF5_VERSION_major=${HDF5_VERSION::4}
fi

rm -rf ${hdf5_build_dir}/bld ${hdf5_install_dir}
mkdir -p ${hdf5_build_dir}/bld
cd ${hdf5_build_dir}
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5_VERSION_major}/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz
tar -xzvf hdf5-${HDF5_VERSION}.tar.gz
ln -snf hdf5-${HDF5_VERSION} src
cd bld
../src/configure --enable-shared \
                 --prefix=${hdf5_install_dir} \
                 CC=${CC} CXX=${CXX} FC=${FC}
make -j${jobs}
make install
cd
rm -rf ${hdf5_build_dir}
