#!/bin/bash

# $1: compiler (gcc, clang)
# $2: moab version (5.1.0, master)
# $3: build static (OFF, ON)

set -e

source /root/etc/env.sh $1
moab_version=$2
build_static=$3

if [ "${TRAVIS_PULL_REQUEST}" == "false" ]; then
  build_mw_reg_tests=ON
else
  build_mw_reg_tests=OFF
fi

rm -rf ${build_dir}/DAGMC-moab-${moab_version}/bld ${install_dir}/DAGMC-moab-${moab_version}
mkdir -p ${build_dir}/DAGMC-moab-${moab_version}/bld
cd ${build_dir}/DAGMC-moab-${moab_version}
#git clone https://github.com/svalinn/DAGMC -b develop
ln -snf DAGMC src
cd bld
cmake ../src -DMOAB_DIR=${install_dir}/moab-${moab_version} \
             -DBUILD_GEANT4=ON \
             -DGEANT4_DIR=${install_dir}/geant4-${geant4_version} \
             -DBUILD_CI_TESTS=ON \
             -DBUILD_MW_REG_TESTS=${build_mw_reg_tests} \
             -DBUILD_STATIC_EXE=${build_static} \
             -DCMAKE_C_COMPILER=${CC} \
             -DCMAKE_CXX_COMPILER=${CXX} \
             -DCMAKE_Fortran_COMPILER=${FC} \
             -DCMAKE_INSTALL_PREFIX=${install_dir}/DAGMC-moab-${moab_version}
make -j${jobs}
make install
cd
#rm -rf ${build_dir}/DAGMC-moab-${moab_version}
