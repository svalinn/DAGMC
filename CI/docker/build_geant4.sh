#!/bin/bash

set -ex


source ${docker_env}

rm -rf ${geant4_build_dir}/bld ${geant4_install_dir}
mkdir -p ${geant4_build_dir}/bld
cd ${geant4_build_dir}
wget https://gitlab.cern.ch/geant4/geant4/-/archive/v${geant4_version}/${geant4_tarball} --no-check-certificate
tar_chashum=$(sha256sum ${geant4_tarball} | cut -d' ' -f1)
if [ $geant4_shasum != $tar_chashum ]; then
    echo "Bad shasum for Geant4 tar!"
    exit 1
fi

tar -xzf ${geant4_tarball}
cd bld
cmake ../${geant4_basename} -DBUILD_STATIC_LIBS=ON \
             -DGEANT4_USE_SYSTEM_EXPAT=OFF \
             -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_C_COMPILER=${CC} \
             -DCMAKE_CXX_COMPILER=${CXX} \
             -DCMAKE_INSTALL_RPATH=${geant4_install_dir}/lib \
             -DCMAKE_INSTALL_PREFIX=${geant4_install_dir}
make -j${ci_jobs}
make install
cd
rm -rf ${geant4_build_dir}
