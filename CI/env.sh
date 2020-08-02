#!/bin/bash

export jobs=`grep -c processor /proc/cpuinfo`

export CTEST_OUTPUT_ON_FAILURE=1

export dagmc_build_dir=${build_dir}/DAGMC
export dagmc_install_dir=${install_dir}/DAGMC
export dagmc_build_dir_shared=${dagmc_build_dir}-shared
export dagmc_install_dir_shared=${dagmc_install_dir}-shared

export dagmc_build_dir_static=${dagmc_build_dir}-static
export dagmc_install_dir_static=${dagmc_install_dir}-static


# Compiler related variables

export  FC=`which gfortran`

if [ "$COMPILER" == "gcc" ]; then
  export  CC=`which gcc`
  export CXX=`which g++`
  export jobs="4"
elif [ "$COMPILER" == "clang" ]; then
  export  CC=`which clang`
  export CXX=`which clang++`
fi
