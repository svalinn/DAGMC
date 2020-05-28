#!/bin/bash

export jobs=`grep -c processor /proc/cpuinfo`

export geant4_version=10.05
export geant4_shasum=4b05b4f7d50945459f8dbe287333b9efb772bd23d50920630d5454ec570b782d

export build_dir=/root/build/${COMPILER}
export install_dir=/root/opt/${COMPILER}

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

export dagmc_build_dir=${build_dir}/DAGMC-moab-${MOAB_VERSION}-hdf5-${HDF5_VERSION}
export dagmc_install_dir=${install_dir}/DAGMC-moab-${MOAB_VERSION}-hdf5-${HDF5_VERSION}

if [ "$COMPILER" == "gcc" ]; then
  export  CC=`which gcc`
  export CXX=`which g++`
  export  FC=`which gfortran`
elif [ "$COMPILER" == "clang" ]; then
  export  CC=`which clang`
  export CXX=`which clang++`
  export  FC=`which gfortran`
fi
