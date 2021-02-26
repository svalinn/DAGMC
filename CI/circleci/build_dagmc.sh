#!/bin/bash

# DAGMC source code already exists because it was copied into
# the docker image by circle.yml

set -ex

source ${docker_env}

function build_dagmc() {




  if [ "$1" == "shared" ]; then
    local build_dir=${dagmc_build_dir_shared}
    local install_dir=${dagmc_install_dir_shared}
    local static_exe=OFF
    local SHARED=ON
    local STATIC=OFF
    if [ "$COVERAGE" == "ON" ]; then
      local ADDITIONAL_CMAKE_FLAGS="-DCOVERAGE=ON -DCMAKE_BUILD_TYPE=Debug"
    fi
  else
    local SHARED=OFF
    local STATIC=ON
    local build_dir=${dagmc_build_dir_static}
    local install_dir=${dagmc_install_dir_static}
    local static_exe=ON
  fi

  if [ "${PULL_REQUEST}" == "false" ]; then
    local build_mw_reg_tests=ON
  else
    local build_mw_reg_tests=OFF
  fi

  if [ "${DOUBLE_DOWN}" == "ON" ]; then
    local double_down=ON
  else
    local double_down=OFF
  fi

  rm -rf ${build_dir}/bld ${install_dir}
  mkdir -p ${build_dir}/bld
  cd ${build_dir}
  cd bld
  cmake ${dagmc_build_dir} ${ADDITIONAL_CMAKE_FLAGS} \
               -DBUILD_SHARED_LIBS=${SHARED} \
               -DBUILD_STATIC_LIBS=$STATIC} \
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
  make -j${jobs}
  make install
}
if [ "$COVERAGE" != "ON" ]; then
  build_dagmc "static"
fi
build_dagmc "shared"

cd
