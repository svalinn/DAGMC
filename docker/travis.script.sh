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
  # Check for news file if this is a PR
  if [ "${TRAVIS_PULL_REQUEST}" != "false" ]; then
    news_file=$(printf 'news/PR-%04u' $TRAVIS_PULL_REQUEST)
    if [ ! -f "${news_file}" ]; then
      echo "Error: news file ${news_file} not found."
      exit 1
    fi
  fi

  # Run astyle
  bash docker/run_astyle.sh

  # Build documentation
  make html

  exit 0
fi

# Run the tests
bash docker/run_tests.sh ${compiler} ${moab_version}
