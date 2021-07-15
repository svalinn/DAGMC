#!/bin/bash

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
