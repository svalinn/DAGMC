#!/bin/bash

# $1: astyle only (OFF, ON)
# $2: compiler (gcc, clang)
# $3: moab version (5.1.0, master)

set -e

astyle_only=$1
compiler=$2
moab_version=$3

cd /root/build/${compiler}/DAGMC-moab-${moab_version}/DAGMC

# If only running astyle, don't need to build anything
if [ "${astyle_only}" == "ON" ]; then
  exit 0
fi

# Only build MOAB master; v5.1.0 is already in the docker image
if [ "${moab_version}" == "master" ]; then
  bash docker/build_moab.sh ${compiler} ${moab_version}
fi

# Build DAGMC (shared executables)
bash docker/build_dagmc.sh ${compiler} ${moab_version} OFF

# Build DAGMC (static executables)
bash docker/build_dagmc.sh ${compiler} ${moab_version} ON
