#!/bin/bash

set -ex

source ${docker_env}


build_dir=${dagmc_build_dir}
install_dir=${dagmc_install_dir_shared}
static_exe=OFF


rm -rf ${build_dir}/bld ${install_dir}
mkdir -p ${build_dir}/bld
cd ${build_dir}
cd bld
cmake ${dagmc_build_dir} -DCOVERALLS=ON \
            -DCMAKE_BUILD_TYPE=Debug \
            -DMOAB_DIR=${moab_install_dir} \
            -DBUILD_GEANT4=ON \
            -DGEANT4_DIR=${geant4_install_dir} \
            -DBUILD_CI_TESTS=ON \
            -DBUILD_MW_REG_TESTS=${build_mw_reg_tests} \
            -DBUILD_STATIC_EXE=${static_exe} \
            -DCMAKE_C_COMPILER=${CC} \
            -DCMAKE_CXX_COMPILER=${CXX} \
            -DCMAKE_Fortran_COMPILER=${FC} \
            -DCMAKE_INSTALL_PREFIX=${install_dir} \
            -DDOUBLE_DOWN=${double_down} \
            -Ddd_ROOT=${double_down_install_dir}
make -j8
make coveralls