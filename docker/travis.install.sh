#!/bin/bash

# $1: compiler (gcc, clang)
# $2: moab version (5.1.0, master)

set -e

export compiler=$1
export moab_version=$2

source /root/etc/env.sh

cd ${dagmc_build_dir}/DAGMC

# Only build MOAB master; v5.1.0 is already in the docker image
if [ "${moab_version}" == "master" ]; then
  bash docker/build_moab.sh ${compiler} ${moab_version}
fi

# Build DAGMC (shared executables)
bash docker/build_dagmc.sh ${compiler} ${moab_version} OFF

# Build DAGMC (static executables)
bash docker/build_dagmc.sh ${compiler} ${moab_version} ON
