#!/bin/bash

# $1: compiler (gcc-4.8, gcc-5, gcc-6, clang-3.8)

set -e

source /root/etc/$1.env

mkdir -p ${build_dir}/geant4-${geant4_version}/bld
rm -rf ${install_dir}/geant4-${geant4_version}
cd ${build_dir}/geant4-${geant4_version}
wget http://geant4.cern.ch/support/source/geant4.${geant4_version}.tar.gz
tar -xzvf geant4.${geant4_version}.tar.gz
ln -snf geant4.${geant4_version} src
cd bld
cmake ../src -DBUILD_STATIC_LIBS=ON \
             -DCMAKE_C_COMPILER=${CC} \
             -DCMAKE_CXX_COMPILER=${CXX} \
             -DCMAKE_INSTALL_PREFIX=${install_dir}/geant4-${geant4_version}
make -j`grep -c processor /proc/cpuinfo`
make install
rm -rf ${build_dir}/geant4-${geant4_version}
