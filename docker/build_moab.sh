#!/bin/bash

mkdir -p ${build_dir}/moab-${moab_version}/bld && \
cd ${build_dir}/moab-${moab_version} && \
git clone https://bitbucket.org/fathomteam/moab -b Version${moab_version} && \
ln -snf moab src && \
cd moab && \
autoreconf -fi && \
cd ../bld && \
../src/configure --enable-dagmc \
                 --disable-ahf \
                 --enable-shared \
                 --enable-optimize \
                 --disable-debug \
                 --with-hdf5=${install_dir}/hdf5-${hdf5_version} \
                 --prefix=${install_dir}/moab-${moab_version} \
                 CC=${CC} CXX=${CXX} FC=${FC} && \
make -j8 && \
make install && \
rm -rf ${build_dir}
