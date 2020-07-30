#!/bin/bash
set -ex

# Embree Variables
EMBREE_TAG='v3.6.1'
EMBREE_REPO='https://github.com/embree/embree'
EMBREE_INSTALL_DIR=$HOME/EMBREE/

CURRENT_DIR=$(pwd)

# Embree Install
cd $HOME
mkdir EMBREE && cd EMBREE
git clone -b $EMBREE_TAG $EMBREE_REPO
mkdir build && cd build
cmake ../embree -DCMAKE_INSTALL_PREFIX=$EMBREE_INSTALL_DIR \
      -DEMBREE_ISPC_SUPPORT=OFF \
      -DEMBREE_TASKING_SYSTEM=INTERNAL \
      -DEMBREE_TUTORIALS=OFF \
      -DEMBREE_TBB_ROOT=/usr
make -j2 && make -j2 install
rm -rf $HOME/EMBREE/embree
