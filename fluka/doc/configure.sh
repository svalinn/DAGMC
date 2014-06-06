#! /bin/bash
# Any changes to this script may require changes to
# batlab scripts in repo svalinn/DAGMC-CI

EXTRA_ARGS=$@

rm -rf CMakeCache.txt
rm tests/slabs.h5m
# Protect from removing the real src and tests dirs
rm -rf ../bld/src
rm -rf ../bld/tests

cmake \
-D MOAB_HOME=$HOME/data/opt/dagmc_bld/moab \
-D HDF5_HOME=$HOME/dagmc_bld/HDF5 \
-D PYNE_HOME=$HOME/.local/lib/python2.7/site-packages/pyne \
$EXTRA_ARGS \
..
