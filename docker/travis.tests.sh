#!/bin/bash

set -e

source ${docker_env}

# If this is not a pull request, get files needed for regression tests
if [ "${TRAVIS_PULL_REQUEST}" == "false" ]; then
  cd /tmp
  wget ${MW_REG_TEST_MODELS_URL} -O mw_reg_test_files.tar.gz -o wget.out
  tar xzvf mw_reg_test_files.tar.gz
  cd -
fi

# Build DAGMC and test (shared executables)
# GEANT4's internal RPATH's aren't quite right,
# so we need to set the LD_LIBRARY_PATH for the
# test to run
cd ${dagmc_build_dir}/DAGMC
LD_LIBRARY_PATH=${geant4_install_dir}/lib:$LD_LIBRARY_PATH \
build_static_exe=OFF docker/build_dagmc.sh

cd ${dagmc_build_dir}/bld/src/make_watertight/tests/
for i in /tmp/*.h5m; do ln -sf $i; done
cd ${dagmc_build_dir}/bld
CTEST_OUTPUT_ON_FAILURE=1 make test

# Build DAGMC and test (static executables)
cd ${dagmc_build_dir}/DAGMC
build_static_exe=ON docker/build_dagmc.sh

cd ${dagmc_build_dir}/bld/src/make_watertight/tests/
for i in /tmp/*.h5m; do ln -sf $i; done
cd ${dagmc_build_dir}/bld
CTEST_OUTPUT_ON_FAILURE=1 make test

# clean out config test directory for next build
git clean -dxf .

# Test DAGMC CMake configuration file
cd ${dagmc_build_dir}/DAGMC/cmake/test_config
cmake . -DDAGMC_ROOT=${dagmc_install_dir}
for i in /tmp/*.h5m; do ln -sf src/make_watertight/tests/$i; done
make all test

cd ${dagmc_build_dir}/bld

# Delete regression test files
if [ "${TRAVIS_PULL_REQUEST}" == "false" ]; then
  rm -f /tmp/*.h5m
  rm -f /tmp/mw_reg_test_files.tar.gz
fi
