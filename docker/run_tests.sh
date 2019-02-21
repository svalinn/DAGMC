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

# Run the unit tests
cd ${build_dir}/DAGMC-moab-${moab_version}/bld
make test

# If this is not a pull request, run regression tests
if [ ! -z ${TRAVIS_PULL_REQUEST} ] && [ ${TRAVIS_PULL_REQUEST} == "false" ]; then
  wget ${MW_REG_TEST_MODELS_URL} -O mw_reg_test_files.tar.gz -o wget.out
  tar xzvf mw_reg_test_files.tar.gz
  ./make_watertight_regression_tests
  rm -rf *.h5m mw_reg_test_files.tar.gz
  if [ $? != 0 ]; then
    exit 1
  fi
fi
