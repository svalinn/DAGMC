#!/bin/bash

export jobs=`grep -c processor /proc/cpuinfo`

export geant4_version=10.05

export build_dir=/root/build/${compiler}
export install_dir=/root/opt/${compiler}

if [ "${hdf5_version}" == "system" ]; then
  export hdf5_build_dir=
  export hdf5_install_dir=/usr/lib/x86_64-linux-gnu/hdf5/serial
else
  export hdf5_build_dir=${build_dir}/hdf5-${hdf5_version}
  export hdf5_install_dir=${install_dir}/hdf5-${hdf5_version}
fi

export geant4_build_dir=${build_dir}/geant4-${geant4_version}
export geant4_install_dir=${install_dir}/geant4-${geant4_version}

export moab_build_dir=${build_dir}/moab-${moab_version}-hdf5-${hdf5_version}
export moab_install_dir=${install_dir}/moab-${moab_version}-hdf5-${hdf5_version}

export dagmc_build_dir=${build_dir}/DAGMC-moab-${moab_version}-hdf5-${hdf5_version}
export dagmc_install_dir=${install_dir}/DAGMC-moab-${moab_version}-hdf5-${hdf5_version}

if [ "$1" == "gcc" ]; then
  export CC=gcc
  export CXX=g++
  export FC=gfortran
elif [ "$1" == "clang" ]; then
  export CC=clang
  export CXX=clang++
  export FC=gfortran
fi
