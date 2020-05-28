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

# symlink make_watertight test files if present
# shared
cd ${dagmc_build_dir_shared}/bld/src/make_watertight/tests/
for i in /tmp/*.h5m; do ln -sf $i; done
# static
cd ${dagmc_build_dir_static}/bld/src/make_watertight/tests/
for i in /tmp/*.h5m; do ln -sf $i; done


# Test DAGMC (shared executables)
cd ${dagmc_build_dir_shared}/bld
PATH=${dagmc_install_dir_shared}/bin:$PATH make test

# Test DAGMC (static executables)
cd ${dagmc_build_dir_static}/bld
PATH=${dagmc_install_dir_static}/bin:PATH make test

# clean out config test directory for next build
cd ${dagmc_build_dir}/DAGMC
git clean -dxf .

# Test DAGMC CMake configuration file
cd ${dagmc_build_dir}/DAGMC/cmake/test_config
cmake . -DDAGMC_ROOT=${dagmc_install_dir_shared}
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
