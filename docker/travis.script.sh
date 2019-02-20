#!/bin/bash

set -e

cd /root/build/${COMPILER}/DAGMC-moab-${MOAB_VERSION}/DAGMC

# Only run astyle
if [ "${ASTYLE_ONLY}" == "ON" ]; then
  bash docker/run_astyle.sh
  exit 0
fi

# Run the tests
export TRAVIS_PULL_REQUEST=${TRAVIS_PULL_REQUEST}
export MW_REG_TEST_MODELS_URL=${MW_REG_TEST_MODELS_URL}
bash docker/run_tests.sh ${COMPILER} ${MOAB_VERSION}
