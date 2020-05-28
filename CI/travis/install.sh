#!/bin/bash

set -ex

source ${docker_env}

cd ${dagmc_build_dir}/DAGMC

# Only build MOAB master and develop; v5.1.0 is already in the docker image
if [ "${MOAB_VERSION}" == "master" ] || [ "${MOAB_VERSION}" == "develop" ]; then
  CI/docker/build_moab.sh
fi
