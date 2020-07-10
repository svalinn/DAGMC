#!/bin/bash

set -ex

source ${docker_env}

# symlink make_watertight regression test files if present
# shared
cd ${dagmc_build_dir_shared}/bld/src/make_watertight/tests/
for i in /tmp/*.h5m; do ln -sf $i; done
ls .
# static
cd ${dagmc_build_dir_static}/bld/src/make_watertight/tests/
for i in /tmp/*.h5m; do ln -sf $i; done
ls .

# Test DAGMC (shared executables)
cd ${dagmc_build_dir_shared}/bld
PATH=${dagmc_install_dir_shared}/bin:$PATH make test

# Test DAGMC (static executables)
cd ${dagmc_build_dir_static}/bld
PATH=${dagmc_install_dir_static}/bin:$PATH make test

# clean out config test directory for next build
cd ${dagmc_build_dir}/DAGMC
git clean -dxf .

# Test DAGMC CMake configuration file
cd ${dagmc_build_dir}/DAGMC/cmake/test_config
cmake . -DDAGMC_ROOT=${dagmc_install_dir_shared}
make all test

# Delete regression test files
if [ "${TRAVIS_PULL_REQUEST}" == "false" ]; then
  rm -f /tmp/*.h5m
  rm -f /tmp/mw_reg_test_files.tar.gz
fi
