#!/bin/bash

set -e

source ${docker_env}

cd ${dagmc_build_dir}/DAGMC

# Only build MOAB master and develop; v5.1.0 is already in the docker image
if [ "${MOAB_VERSION}" == "master" ] || [ "${MOAB_VERSION}" == "develop" ]; then
  docker/build_moab.sh
fi

# Build the double_down project if
if [ "${DOUBLE_DOWN}" == "yes" ]; then
  docker/build_double_down.sh
fi

# Build DAGMC (shared executables)
build_static_exe=OFF docker/build_dagmc.sh

# Build DAGMC (static executables)
if [ "${DOUBLE_DOWN}" != "yes"]; then
  build_static_exe=ON docker/build_dagmc.sh
fi
