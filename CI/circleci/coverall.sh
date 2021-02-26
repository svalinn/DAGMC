#!/bin/bash

set -ex

source ${docker_env}

cd ${dagmc_build_dir_shared}/bld
/root/.local/bin/coveralls --root . -E ".*gtest.*" -E ".*CMakeFiles.*"