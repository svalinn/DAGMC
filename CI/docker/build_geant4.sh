#!/bin/bash

set -ex

source ${docker_env}

rm -rf ${geant4_build_dir}/bld ${geant4_install_dir}
mkdir -p ${geant4_build_dir}/bld
cd ${geant4_build_dir}
wget https://geant4.cern.ch/support/source/geant4.${geant4_version}.tar.gz --no-check-certificate
tar_chashum=$(sha256sum ${geant4_version}.tar.gz | cut -d' ' -f1)
if [ $geant4_shasum != $tar_chashum ]; then
    echo "Bad shasum for Geant4 tar!"
    exit 1
fi

tar -xzvf geant4.${geant4_version}.tar.gz
cd bld
cmake ../geant4.${geant4_version} -DBUILD_STATIC_LIBS=ON \
             -DGEANT4_USE_SYSTEM_EXPAT=OFF \
             -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_C_COMPILER=${CC} \
             -DCMAKE_CXX_COMPILER=${CXX} \
             -DCMAKE_INSTALL_PREFIX=${geant4_install_dir}
make -j${jobs}
make install
cd
rm -rf ${geant4_build_dir}
