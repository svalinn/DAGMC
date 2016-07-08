#!/bin/bash
MOAB_VERSION=$1
start_dir=$PWD
cd /root
git clone https://bitbucket.org/fathomteam/moab
cd moab
git checkout $MOAB_VERSION
autoreconf -fi
mkdir bld
cd bld
/./configure --enable-dagmc --enable-shared --disable-debug --enable-optimize --with-hdf5 --prefix=/root/moab
make -j2
make install
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$root/moab/lib/"
export PATH="/root/moab/bin:$PATH"
cd $start_dir
