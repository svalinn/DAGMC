#!/bin/bash

# $1: compiler (gcc-4.8, gcc-5, gcc-6, gcc-7, clang-4.0, clang-5.0)
# $2: moab version (5.1.0, master)

set -e

source /root/etc/$1.env
moab_version=$2
hdf5_version=1.8.13

mkdir -p ${build_dir}/moab-${moab_version}/bld
rm -rf ${install_dir}/moab-${moab_version}
cd ${build_dir}/moab-${moab_version}
if [[ ${moab_version} == "master" ]]; then
  git clone https://bitbucket.org/fathomteam/moab -b ${moab_version}
else
  git clone https://bitbucket.org/fathomteam/moab -b Version${moab_version}
fi
ln -snf moab src
cd moab
autoreconf -fi
cd ../bld
../src/configure --enable-dagmc \
                 --disable-ahf \
                 --enable-shared \
                 --enable-optimize \
                 --disable-debug \
                 --with-hdf5=${install_dir}/hdf5-${hdf5_version} \
                 --prefix=${install_dir}/moab-${moab_version} \
                 CC=${CC} CXX=${CXX} FC=${FC}
make -j`grep -c processor /proc/cpuinfo`
make install
rm -rf ${build_dir}/moab-${moab_version}
