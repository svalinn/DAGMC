#!/bin/bash

# $1: housekeeping only (OFF, ON)
# $2: compiler (gcc, clang)
# $3: moab version (5.1.0, master)

set -e

hk_only=$1
compiler=$2
moab_version=$3

source /root/etc/env.sh ${compiler}

cd ${build_dir}/DAGMC-moab-${moab_version}/DAGMC

# Only do housekeeping tasks
if [ "${hk_only}" == "ON" ]; then
  # Run astyle
  bash docker/run_astyle.sh
  # Build documentation
  make html
  exit 0
fi

# Run the tests
export TRAVIS_PULL_REQUEST=${TRAVIS_PULL_REQUEST}
export MW_REG_TEST_MODELS_URL=${MW_REG_TEST_MODELS_URL}
bash docker/run_tests.sh ${compiler} ${moab_version}
