#!/bin/bash

# Manual testing that mimics CI
# --------------------------------
# 
# Run `CI/circleci/install.sh` in one of the CI Docker images with the root dir 
# of the DAGMC repo as the working directory.  Note that CI only runs on 
# MOAB_VERSION=5.1.0 and HDF5_VERSION=1.10.4
#
# If you are in the DAGMC root directory:
#
#    $ docker run -v ${PWD}:/root/build_dir/DAGMC -w /root/build_dir/DAGMC \
#              -it svalinn/dagmc-ci-ubuntu-16.04-clang-ext-hdf5_1.10.4-moab_5.1.0
#
# Inside the docker container:
#
#    $ CI/circleci/install.sh
#
#    $ CI/circleci/tests.sh

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
cd ${dagmc_build_dir}
git clean -dxf .

# Test DAGMC CMake configuration file
cd ${dagmc_build_dir}/cmake/test_config
cmake . -DDAGMC_ROOT=${dagmc_install_dir_shared}
make all test

# Delete regression test files
if [ "${PULL_REQUEST}" == "false" ]; then
  rm -f /tmp/*.h5m
  rm -f /tmp/mw_reg_test_files.tar.gz
fi
