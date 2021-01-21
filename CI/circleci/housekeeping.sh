#!/bin/bash

set -ex

source ${docker_env}

cd ${dagmc_build_dir}

# Check for news file if this is a PR into svalinn/DAGMC
if [ "${REPO_SLUG}" == "svalinn/DAGMC" ] && \
   [ "${PULL_REQUEST}" != "false" ]; then
  news_file=$(printf 'news/PR-%04u.rst' ${PULL_REQUEST})
  if [ -f "${news_file}" ]; then
    echo "News file ${news_file} found!"
  else
    echo "ERROR: News file ${news_file} not found. Please create a news file."
    exit 1
  fi
fi

# Run clang-format check
find src/ \( -name "*.hpp" -o -name "*.cpp" -o -name "*.hh" -o -name "*.cc" -o -name "*.h" \) \
          \( -not -path "src/gtest*" -not -path "src/mcnp/mcnp?/Source/*" -not -path "src/pyne*" \)  \
          -exec clang-format -style=file -i {} \;
clang_diffs=`git status --porcelain`
if [ -z "${clang_diffs}" ]; then
  echo "Style guide checker passed!"
else
  echo "ERROR: Style guide checker failed. Please run clang-format."
  echo "clang_diffs: ${clang_diffs}"
  exit 1
fi

# Build documentation
make html
