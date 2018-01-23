#!/bin/bash

ASTYLE_ONLY=$1
COMPILER=$2
MOAB_VERSION=$3
BUILD_GEANT4=$4

cd /root/build/${COMPILER}/DAGMC-moab-${MOAB_VERSION}/DAGMC

if [ "${ASTYLE_ONLY}" == "ON" ]; then
  # Run astyle
  bash docker/run_astyle.sh
else
  # Run the tests
  export TRAVIS_PULL_REQUEST=${TRAVIS_PULL_REQUEST}
  export MW_REG_TEST_MODELS_URL=${MW_REG_TEST_MODELS_URL}
  bash docker/run_tests.sh ${COMPILER} ${MOAB_VERSION} ${BUILD_GEANT4}
fi
