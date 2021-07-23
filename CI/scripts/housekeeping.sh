#!/bin/bash

set -ex

source ${docker_env}

cd ${dagmc_build_dir}

# Ensure clang-format is there
if ! command -v clang-format &> /dev/null
then
    echo "clang-format could not be found"
    exit 1
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
