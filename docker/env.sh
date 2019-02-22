#!/bin/bash

export jobs=`grep -c processor /proc/cpuinfo`

export geant4_version=10.05

export build_dir=/root/build/${COMPILER}
export install_dir=/root/opt/${COMPILER}

if [ "${HDF5_VERSION}" == "system" ]; then
  export hdf5_build_dir=
  export hdf5_install_dir=/usr/lib/x86_64-linux-gnu/hdf5/serial
else
  export hdf5_build_dir=${build_dir}/hdf5-${HDF5_VERSION}
  export hdf5_install_dir=${install_dir}/hdf5-${HDF5_VERSION}
fi

export geant4_build_dir=${build_dir}/geant4-${geant4_version}
export geant4_install_dir=${install_dir}/geant4-${geant4_version}

export moab_build_dir=${build_dir}/moab-${MOAB_VERSION}-hdf5-${HDF5_VERSION}
export moab_install_dir=${install_dir}/moab-${MOAB_VERSION}-hdf5-${HDF5_VERSION}

export dagmc_build_dir=${build_dir}/DAGMC-moab-${MOAB_VERSION}-hdf5-${HDF5_VERSION}
export dagmc_install_dir=${install_dir}/DAGMC-moab-${MOAB_VERSION}-hdf5-${HDF5_VERSION}

if [ "$1" == "gcc" ]; then
  export CC=gcc
  export CXX=g++
  export FC=gfortran
elif [ "$1" == "clang" ]; then
  export CC=clang
  export CXX=clang++
  export FC=gfortran
fi
