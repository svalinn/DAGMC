#!/bin/bash

# $1: astyle only (OFF, ON)
# $2: documentation only (OFF, ON)
# $3: compiler (gcc, clang)
# $4: moab version (5.1.0, master)

set -e

astyle_only=$1
doc_only=$2
compiler=$3
moab_version=$4

source /root/etc/env.sh ${compiler}

# If only running astyle or only building documentation, don't need to build
# anything
if [ "${astyle_only}" == "ON" ] || [ "${doc_only}" == "ON" ]; then
  exit 0
fi

cd ${build_dir}/DAGMC-moab-${moab_version}/DAGMC

# Only build MOAB master; v5.1.0 is already in the docker image
if [ "${moab_version}" == "master" ]; then
  bash docker/build_moab.sh ${compiler} ${moab_version}
fi

# Build DAGMC (shared executables)
bash docker/build_dagmc.sh ${compiler} ${moab_version} OFF

# Build DAGMC (static executables)
bash docker/build_dagmc.sh ${compiler} ${moab_version} ON
