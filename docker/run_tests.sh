#!/bin/bash

# $1: compiler (gcc, clang)
# $2: moab version (5.1.0, master)

set -e

source /root/etc/env.sh $1
moab_version=$2

# clean out config test directory for next build
git clean -dxf . 

# Test DAGMC CMake configuration file
cd ${build_dir}/DAGMC-moab-${moab_version}/src/cmake/test_config
cmake . -DDAGMC_ROOT=${install_dir}/DAGMC-moab-${moab_version}
make all test

cd ${build_dir}/DAGMC-moab-${moab_version}/bld

# If this is not a pull request, get files needed for regression tests
if [ ${TRAVIS_PULL_REQUEST} == "false" ]; then
  cd src/make_watertight/tests
  wget ${MW_REG_TEST_MODELS_URL} -O mw_reg_test_files.tar.gz -o wget.out
  tar xzvf mw_reg_test_files.tar.gz
  cd ../../..
fi

# Run the unit tests
make test

# Delete regression test files
if [ ${TRAVIS_PULL_REQUEST} == "false" ]; then
  rm -f src/make_watertight/tests/*.h5m
  rm -f src/make_watertight/tests/mw_reg_test_files.tar.gz
fi
