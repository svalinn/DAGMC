#!/bin/bash

mkdir -p ${build_dir}/hdf5-${hdf5_version}/bld && \
cd ${build_dir}/hdf5-${hdf5_version} && \
wget https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-${hdf5_version}/src/hdf5-${hdf5_version}.tar.gz && \
tar -xzvf hdf5-${hdf5_version}.tar.gz && \
ln -snf hdf5-${hdf5_version} src && \
cd bld && \
../src/configure --enable-shared \
                 --disable-debug \
                 --prefix=${install_dir}/hdf5-${hdf5_version} \
                 CC=${CC} CXX=${CXX} FC=${FC} && \
make -j8 && \
make install && \
rm -rf ${build_dir}
