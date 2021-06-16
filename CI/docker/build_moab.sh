#!/bin/bash

set -ex

source ${docker_env}


rm -rf ${moab_build_dir}/bld ${moab_install_dir}
mkdir -p ${moab_build_dir}/bld
cd ${moab_build_dir}
git clone --depth 1 https://bitbucket.org/fathomteam/moab -b ${branch}
cd moab
autoreconf -fi
cd ../bld
../moab/configure --enable-pymoab \
                 --enable-shared \
                 --enable-optimize \
                 --disable-debug \
                 --disable-fortran \
                 --disable-blaslapack \
                 --with-hdf5=${hdf5_install_dir} \
                 --prefix=${moab_install_dir} \
                 CC=${CC} CXX=${CXX}
make -j${ci_jobs}
make install
cd
rm -rf ${moab_build_dir}
