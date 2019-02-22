#!/bin/bash

set -e

cd ${dagmc_build_dir}/DAGMC

# clean out config test directory for next build
git clean -dxf . 

# Test DAGMC CMake configuration file
cd ${dagmc_build_dir}/DAGMC/cmake/test_config
cmake . -DDAGMC_ROOT=${dagmc_install_dir}
make all test

cd ${dagmc_build_dir}/bld

# If this is not a pull request, get files needed for regression tests
if [ "${TRAVIS_PULL_REQUEST}" == "false" ]; then
  cd src/make_watertight/tests
  wget ${MW_REG_TEST_MODELS_URL} -O mw_reg_test_files.tar.gz -o wget.out
  tar xzvf mw_reg_test_files.tar.gz
  cd ../../..
fi

# Run the unit tests
make test

# Delete regression test files
if [ "${TRAVIS_PULL_REQUEST}" == "false" ]; then
  rm -f src/make_watertight/tests/*.h5m
  rm -f src/make_watertight/tests/mw_reg_test_files.tar.gz
fi
