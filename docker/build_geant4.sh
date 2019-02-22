#!/bin/bash

set -e

rm -rf ${geant4_build_dir}/bld ${geant4_install_dir}
mkdir -p ${geant4_build_dir}/bld
cd ${geant4_build_dir}
wget http://geant4.cern.ch/support/source/geant4.${geant4_version}.tar.gz
tar -xzvf geant4.${geant4_version}.tar.gz
ln -snf geant4.${geant4_version} src
cd bld
cmake ../src -DBUILD_STATIC_LIBS=ON \
             -DGEANT4_USE_SYSTEM_EXPAT=OFF \
             -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_C_COMPILER=${CC} \
             -DCMAKE_CXX_COMPILER=${CXX} \
             -DCMAKE_INSTALL_PREFIX=${geant4_install_dir}
make -j${jobs}
make install
cd
rm -rf ${geant4_build_dir}
