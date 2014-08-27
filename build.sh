#!/bin/bash
geant4_dir="/mnt/data/opt/geant4/"
g++ g4mat.cpp -I/home/davisa/.local/include/pyne -I/mnt/data/opt/dagmc/hdf5/include -I${geant4_dir}/include -I${geant4_dir}/include/Geant4 -L/mnt/data/opt/moab_dev/lib -L/mnt/data/opt/dagmc/hdf5/lib -L/mnt/data/opt/geant4/lib -L/home/davisa/.local/lib -lhdf5 -lpyne -lG4materials -lG4global -lG4clhep -lG4intercoms -o g4mat
