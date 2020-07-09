#!/bin/bash

# DAGMC source code already exists because it was copied into
# the docker image by travis.yml

set -ex

source CI/env.sh

function build_dagmc() {

  if [ "$1" == "shared" ]; then
    local build_dir=${dagmc_build_dir_shared}
    local install_dir=${dagmc_install_dir_shared}
    local static_exe=OFF
  else
    local build_dir=${dagmc_build_dir_static}
    local install_dir=${dagmc_install_dir_static}
    local static_exe=ON
  fi

  if [ "${TRAVIS_PULL_REQUEST}" == "false" ]; then
    local build_mw_reg_tests=ON
  else
    local build_mw_reg_tests=OFF
  fi

  rm -rf ${build_dir}/bld ${install_dir}
  mkdir -p ${build_dir}/bld
  cd ${build_dir}
  cd bld
  cmake ${dagmc_build_dir} -DMOAB_DIR=${moab_install_dir} \
               -DBUILD_GEANT4=ON \
               -DGEANT4_DIR=${geant4_install_dir} \
               -DBUILD_CI_TESTS=ON \
               -DBUILD_MW_REG_TESTS=${build_mw_reg_tests} \
               -DBUILD_STATIC_EXE=${static_exe} \
               -DCMAKE_C_COMPILER=${CC} \
               -DCMAKE_CXX_COMPILER=${CXX} \
               -DCMAKE_Fortran_COMPILER=${FC} \
               -DCMAKE_INSTALL_PREFIX=${install_dir}
  make -j${jobs}
  make install
}

build_dagmc "static"

build_dagmc "shared"

cd
