#!/bin/bash

export jobs=`grep -c processor /proc/cpuinfo`

export hdf5_version=1.10.4
export geant4_version=10.05

export build_dir=/root/build/$1
export install_dir=/root/opt/$1

if [ "$1" == "gcc" ]; then
  export CC=gcc
  export CXX=g++
  export FC=gfortran
elif [ "$1" == "clang" ]; then
  export CC=clang
  export CXX=clang++
  export FC=gfortran
fi
