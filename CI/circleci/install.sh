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
#              -it svalinn/dagmc-ci-ubuntu-16.04-clang-ext-hdf5_1.10.4-moab_5.1.0 \
#              CI/circleci/install.sh
#


set -ex

source ${docker_env}

cd ${dagmc_build_dir}

# Only build MOAB master and develop; v5.1.0 is already in the docker image
if [ "${MOAB_VERSION}" == "master" ] || [ "${MOAB_VERSION}" == "develop" ]; then
  CI/docker/build_moab.sh
fi

# If this is not a pull request, get files needed for regression tests
if [ "${PULL_REQUEST}" == "false" ]; then
  cd /tmp
  wget ${MW_REG_TEST_MODELS_URL} -O mw_reg_test_files.tar.gz -o wget.out
  tar xzvf mw_reg_test_files.tar.gz
  cd -
fi

# Build the double_down project if needed
if [ "${DOUBLE_DOWN}" == "ON" ]; then
  CI/circleci/build_double_down.sh
fi

# Build and install DAGMC
CI/circleci/build_dagmc.sh
