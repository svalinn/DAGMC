#!/bin/bash

ASTYLE_ONLY=$1
COMPILER=$2
MOAB_VERSION=$3
BUILD_GEANT4=$4

if [ "${ASTYLE_ONLY}" == "ON" ]; then
  exit 0
fi

cd /root/build/${COMPILER}/DAGMC-moab-${MOAB_VERSION}/DAGMC

# Build MOAB
bash docker/build_moab.sh ${COMPILER} ${MOAB_VERSION}

# Build DAGMC (shared executables)
bash docker/build_dagmc.sh ${COMPILER} OFF ${MOAB_VERSION} ${BUILD_GEANT4}

# Build DAGMC (static executables)
bash docker/build_dagmc.sh ${COMPILER} ON ${MOAB_VERSION} ${BUILD_GEANT4}
