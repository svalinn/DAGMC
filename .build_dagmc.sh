#!/bin/bash

export LD_LIBRARY_PATH="/root/moab/lib"
mkdir bld
cd bld
cmake ../. -DBUILD_TALLY=ON \
           -DBUILD_TESTS=ON \
           -DBUILD_GEANT4=ON \
           -DGEANT4_DIR=/root/geant4.10.00.p02 \
           -DCMAKE_INSTALL_PREFIX=/root/dagmc
make 
make install
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/root/dagmc/lib"
cd ..

