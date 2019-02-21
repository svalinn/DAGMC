#!/bin/bash

# $1: compiler (gcc, clang)
# $2: moab version (5.1.0, master)
# $3: build static (OFF, ON)

set -e

export compiler=$1
export moab_version=$2
export build_static=$3

source /root/etc/env.sh

if [ "${TRAVIS_PULL_REQUEST}" == "false" ]; then
  build_mw_reg_tests=ON
else
  build_mw_reg_tests=OFF
fi

rm -rf ${dagmc_build_dir}/bld ${dagmc_install_dir}
mkdir -p ${dagmc_build_dir}/bld
cd ${dagmc_build_dir}
#git clone https://github.com/svalinn/DAGMC -b develop
ln -snf DAGMC src
cd bld
cmake ../src -DMOAB_DIR=${moab_install_dir} \
             -DBUILD_GEANT4=ON \
             -DGEANT4_DIR=${geant4_install_dir} \
             -DBUILD_CI_TESTS=ON \
             -DBUILD_MW_REG_TESTS=${build_mw_reg_tests} \
             -DBUILD_STATIC_EXE=${build_static} \
             -DCMAKE_C_COMPILER=${CC} \
             -DCMAKE_CXX_COMPILER=${CXX} \
             -DCMAKE_Fortran_COMPILER=${FC} \
             -DCMAKE_INSTALL_PREFIX=${dagmc_install_dir}
make -j${jobs}
make install
cd
#rm -rf ${dagmc_build_dir}
