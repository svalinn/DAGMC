#!/bin/bash

set -ex

source ${docker_env}

if [ ${MOAB_VERSION} == "master" ] || [ ${MOAB_VERSION} == "develop" ]; then
  branch=${MOAB_VERSION}
else
  branch=Version${MOAB_VERSION}
fi

rm -rf ${moab_build_dir}/bld ${moab_install_dir}
mkdir -p ${moab_build_dir}/bld
cd ${moab_build_dir}
git clone --depth 1 https://bitbucket.org/fathomteam/moab -b ${branch}
cd bld
cmake ../moab -DENABLE_HDF5=ON -DHDF5_ROOT=${hdf5_install_dir} \
              -DENABLE_PYMOAB=ON \
              -DENABLE_BLASLAPACK=OFF \
              -DENABLE_FORTRAN=OFF \
              -DCMAKE_INSTALL_PREFIX=${moab_install_dir} \
              -DCMAKE_C_COMPILER=${CC} \
              -DCMAKE_CXX_COMPILER=${CXX}
make -j${jobs}
make install
cd
rm -rf ${moab_build_dir}
