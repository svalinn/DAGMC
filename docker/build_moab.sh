#!/bin/bash

# $1: compiler (gcc, clang)
# $2: moab version (5.1.0, master)

set -e

source /root/etc/env.sh $1
moab_version=$2

mkdir -p ${build_dir}/moab-${moab_version}/bld
rm -rf ${install_dir}/moab-${moab_version}
cd ${build_dir}/moab-${moab_version}
if [[ ${moab_version} == "master" ]]; then branch=${moab_version}
else branch=Version${moab_version}
fi
git clone https://bitbucket.org/fathomteam/moab -b ${branch}
ln -snf moab src
cd moab
autoreconf -fi
cd ../bld
../src/configure --enable-pymoab \
                 --enable-shared \
                 --enable-optimize \
                 --disable-debug \
                 --with-hdf5=${install_dir}/hdf5-${hdf5_version} \
                 --prefix=${install_dir}/moab-${moab_version} \
                 CC=${CC} CXX=${CXX} FC=${FC}
make -j${jobs}
make install
cd
rm -rf ${build_dir}/moab-${moab_version}
