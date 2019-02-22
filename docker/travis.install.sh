#!/bin/bash

set -e

cd ${dagmc_build_dir}/DAGMC

# Only build MOAB master; v5.1.0 is already in the docker image
if [ "${MOAB_VERSION}" == "master" ]; then
  docker/build_moab.sh
fi

# Build DAGMC (shared executables)
BUILD_STATIC=OFF docker/build_dagmc.sh

# Build DAGMC (static executables)
BUILD_STATIC=ON docker/build_dagmc.sh
