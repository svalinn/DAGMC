#!/bin/bash

# HDF5 VERSIONS VARIABLE
export HDF5_VERSION=1.10.4
export HDF5_VERSION_major=1.10
export hdf5_1_10_4_shasum=8f60dc4dd6ab5fcd23c750d1dc5bca3d0453bdce5c8cdaf0a4a61a9d1122adb2


# Geant4 version and corresponding SHASUM
export geant4_version=10.5.1
export geant4_basename=geant4-v${geant4_version}
export geant4_tarball=${geant4_basename}.tar.gz
export geant4_shasum=2397eb859dc4de095ff66059d8bda9f060fdc42e10469dd7890946293eeb0e39

# Embree Variables
export EMBREE_TAG='v3.6.1'

# MOAB Version variable
# Supported Branches are: 
# - develop
# - master
# - Version5.1.0 used for CI
export  MOAB_BRANCH=${MOAB_VERSION}
