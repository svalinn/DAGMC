#!/bin/bash

export jobs=`grep -c processor /proc/cpuinfo`

export CTEST_OUTPUT_ON_FAILURE=1

if [ "${HDF5_VERSION}" == "system" ]; then
  export hdf5_build_dir=
  export hdf5_install_dir=/usr/lib/x86_64-linux-gnu/hdf5/serial
else
  export hdf5_build_dir=${build_dir}/hdf5-${HDF5_VERSION}
  export hdf5_install_dir=${install_dir}/hdf5-${HDF5_VERSION}
fi

export hdf5_1_10_4_shasum=8f60dc4dd6ab5fcd23c750d1dc5bca3d0453bdce5c8cdaf0a4a61a9d1122adb2


export geant4_build_dir=${build_dir}/geant4-${geant4_version}
export geant4_install_dir=${install_dir}/geant4-${geant4_version}

export moab_build_dir=${build_dir}/moab-${MOAB_VERSION}-hdf5-${HDF5_VERSION}
export moab_install_dir=${install_dir}/moab-${MOAB_VERSION}-hdf5-${HDF5_VERSION}

export double_down_build_dir=${build_dir}/double-down/
export double_down_install_dir=${install_dir}/double-down/

export dagmc_build_dir=${build_dir}/DAGMC-moab-${MOAB_VERSION}-hdf5-${HDF5_VERSION}
export dagmc_install_dir=${install_dir}/DAGMC-moab-${MOAB_VERSION}-hdf5-${HDF5_VERSION}

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
