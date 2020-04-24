#!/bin/bash

set -ex

source ${docker_env}

# If this is not a pull request, get files needed for regression tests
if [ "${TRAVIS_PULL_REQUEST}" == "false" ]; then
  cd /tmp
  wget ${MW_REG_TEST_MODELS_URL} -O mw_reg_test_files.tar.gz -o wget.out
  tar xzvf mw_reg_test_files.tar.gz
  cd -
fi

# Build DAGMC and test (shared executables)
cd ${dagmc_build_dir}/DAGMC
build_static_exe=OFF docker/build_dagmc.sh
cd ${dagmc_build_dir}/bld
make test

# Build DAGMC and test (static executables)
cd ${dagmc_build_dir}/DAGMC
build_static_exe=ON docker/build_dagmc.sh
cd ${dagmc_build_dir}/bld
make test

# clean out config test directory for next build
cd ${dagmc_build_dir}/DAGMC
git clean -dxf .

# Test DAGMC CMake configuration file
cd ${dagmc_build_dir}/DAGMC/cmake/test_config
cmake . -DDAGMC_ROOT=${dagmc_install_dir}
make all test

cd ${dagmc_build_dir}/bld

# If this is not a pull request, get files needed for regression tests
if [ "${TRAVIS_PULL_REQUEST}" == "false" ]; then
  cd src/make_watertight/tests
  wget ${MW_REG_TEST_MODELS_URL} -O mw_reg_test_files.tar.gz
  tar xzvf mw_reg_test_files.tar.gz
  cd ../../..
fi

# Run the unit tests
make test

# Delete regression test files
if [ "${TRAVIS_PULL_REQUEST}" == "false" ]; then
  rm -f /tmp/*.h5m
  rm -f /tmp/mw_reg_test_files.tar.gz
fi
