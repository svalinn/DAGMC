#!/bin/bash

set -e

source ${env_file}

cd ${dagmc_build_dir}/DAGMC

# Check for news file if this is a PR into svalinn/DAGMC
if [ "${TRAVIS_REPO_SLUG}" == "svalinn/DAGMC" ] && \
   [ "${TRAVIS_PULL_REQUEST}" != "false" ]; then
  news_file=$(printf 'news/PR-%04u.rst' ${TRAVIS_PULL_REQUEST})
  if [ -f "${news_file}" ]; then
    echo "News file ${news_file} found!"
  else
    echo "ERROR: News file ${news_file} not found. Please create a news file."
    exit 1
  fi
fi

# Run astyle check
astyle --options=astyle_google.ini \
       --exclude=src/gtest \
       --exclude=src/mcnp/mcnp5/Source \
       --exclude=src/mcnp/mcnp6/Source \
       --ignore-exclude-errors \
       --suffix=none \
       --recursive \
       --verbose \
       --formatted \
       "*.cc" "*.cpp" "*.h" "*.hh" "*.hpp"
astyle_diffs=`git status --porcelain`
if [ -z "${astyle_diffs}" ]; then
  echo "Style guide checker passed!"
else
  echo "ERROR: Style guide checker failed. Please run astyle."
  git diff
  exit 1
fi

# Build documentation
make html

exit 0
