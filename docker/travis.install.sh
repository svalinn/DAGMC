#!/bin/bash

set -e

cd /root/build/${COMPILER}/DAGMC-moab-${MOAB_VERSION}/DAGMC

# If only running astyle, don't need to build anything
if [ "${ASTYLE_ONLY}" == "ON" ]; then
  exit 0
fi

# Build MOAB
bash docker/build_moab.sh  ${COMPILER} ${MOAB_VERSION}

# Build DAGMC (shared executables)
bash docker/build_dagmc.sh ${COMPILER} ${MOAB_VERSION} OFF

# Build DAGMC (static executables)
bash docker/build_dagmc.sh ${COMPILER} ${MOAB_VERSION} ON
