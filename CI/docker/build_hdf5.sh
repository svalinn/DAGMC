#!/bin/bash

set -ex

source ${docker_env}

export hdf5_1_10_4_shasum=8f60dc4dd6ab5fcd23c750d1dc5bca3d0453bdce5c8cdaf0a4a61a9d1122adb2


if [ "${HDF5_VERSION:3:1}" == "." ]; then
  HDF5_VERSION_major=${HDF5_VERSION::3}
elif [ "${HDF5_VERSION:4:1}" == "." ]; then
  HDF5_VERSION_major=${HDF5_VERSION::4}
fi

rm -rf ${hdf5_build_dir}/bld ${hdf5_install_dir}
mkdir -p ${hdf5_build_dir}/bld
cd ${hdf5_build_dir}

wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5_VERSION_major}/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz
tar_chashum=$(sha256sum hdf5-${HDF5_VERSION}.tar.gz | cut -d' ' -f1)
if [ $hdf5_1_10_4_shasum != $tar_chashum ]; then
    echo "Bad shasum for hdf5 tar!"
    exit 1
fi

tar -xzf hdf5-${HDF5_VERSION}.tar.gz
cd bld
../hdf5-${HDF5_VERSION}/configure --enable-shared \
                 --prefix=${hdf5_install_dir} \
                 CC=${CC} CXX=${CXX}
make -j${ci_jobs}
make install
cd
rm -rf ${hdf5_build_dir}
