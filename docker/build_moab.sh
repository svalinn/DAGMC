#!/bin/bash

# $1: compiler (gcc, clang)
# $2: moab version (5.1.0, master)

set -e

export compiler=$1
export moab_version=$2

source /root/etc/env.sh

if [[ ${moab_version} == "master" ]]; then
  branch=${moab_version}
else
  branch=Version${moab_version}
fi

rm -rf ${moab_build_dir}/bld ${moab_install_dir}
mkdir -p ${moab_build_dir}/bld
cd ${moab_build_dir}
git clone https://bitbucket.org/fathomteam/moab -b ${branch}
ln -snf moab src
cd moab
autoreconf -fi
cd ../bld
../src/configure --enable-pymoab \
                 --enable-shared \
                 --enable-optimize \
                 --disable-debug \
                 --with-hdf5=${hdf5_install_dir} \
                 --prefix=${moab_install_dir} \
                 CC=${CC} CXX=${CXX} FC=${FC}
make -j${jobs}
make install
cd
rm -rf ${moab_build_dir}
