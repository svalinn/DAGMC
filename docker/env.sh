#!/bin/bash

export jobs=`grep -c processor /proc/cpuinfo`

export hdf5_version=1.10.4
export geant4_version=10.05

export build_dir=/root/build/${compiler}
export install_dir=/root/opt/${compiler}

export hdf5_build_dir=${build_dir}/hdf5-${hdf5_version}
export hdf5_install_dir=${install_dir}/hdf5-${hdf5_version}

export geant4_build_dir=${build_dir}/geant4-${geant4_version}
export geant4_install_dir=${install_dir}/geant4-${geant4_version}

export moab_build_dir=${build_dir}/moab-${moab_version}
export moab_install_dir=${install_dir}/moab-${moab_version}

export dagmc_build_dir=${build_dir}/DAGMC-moab-${moab_version}
export dagmc_install_dir=${install_dir}/DAGMC-moab-${moab_version}

if [ "$1" == "gcc" ]; then
  export CC=gcc
  export CXX=g++
  export FC=gfortran
elif [ "$1" == "clang" ]; then
  export CC=clang
  export CXX=clang++
  export FC=gfortran
fi
