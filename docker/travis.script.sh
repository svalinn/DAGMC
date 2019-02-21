#!/bin/bash

# $1: astyle only (OFF, ON)
# $2: compiler (gcc, clang)
# $3: moab version (5.1.0, master)

set -e

astyle_only=$1
compiler=$2
moab_version=$3

cd /root/build/${compiler}/DAGMC-moab-${moab_version}/DAGMC

# Only run astyle
if [ "${astyle_only}" == "ON" ]; then
  bash docker/run_astyle.sh
  exit 0
fi

# Run the tests
export TRAVIS_PULL_REQUEST=${TRAVIS_PULL_REQUEST}
export MW_REG_TEST_MODELS_URL=${MW_REG_TEST_MODELS_URL}
bash docker/run_tests.sh ${compiler} ${moab_version}
