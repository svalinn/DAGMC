#!/bin/bash

set -e

cd /root/build/${COMPILER}/DAGMC-moab-${MOAB_VERSION}/DAGMC

# If only running astyle, don't need to build anything
if [ "${ASTYLE_ONLY}" == "ON" ]; then
  exit 0
fi

# Only build MOAB master; v5.1.0 is already in the docker image
if [ "${MOAB_VERSION}" == "master" ]; then
  bash docker/build_moab.sh ${COMPILER} ${MOAB_VERSION}
fi

# Build DAGMC (shared executables)
bash docker/build_dagmc.sh ${COMPILER} ${MOAB_VERSION} OFF

# Build DAGMC (static executables)
bash docker/build_dagmc.sh ${COMPILER} ${MOAB_VERSION} ON
