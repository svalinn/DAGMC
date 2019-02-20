#!/bin/bash

# $1: compiler (gcc, clang)

set -e

source /root/etc/env.sh $1

mkdir -p ${build_dir}/geant4-${geant4_version}/bld
rm -rf ${install_dir}/geant4-${geant4_version}
cd ${build_dir}/geant4-${geant4_version}
wget http://geant4.cern.ch/support/source/geant4.${geant4_version}.tar.gz
tar -xzvf geant4.${geant4_version}.tar.gz
ln -snf geant4.${geant4_version} src
cd bld
cmake ../src -DBUILD_STATIC_LIBS=ON \
             -DGEANT4_USE_SYSTEM_EXPAT=OFF \
             -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_C_COMPILER=${CC} \
             -DCMAKE_CXX_COMPILER=${CXX} \
             -DCMAKE_INSTALL_PREFIX=${install_dir}/geant4-${geant4_version}
make -j${jobs}
make install
cd
rm -rf ${build_dir}/geant4-${geant4_version}
