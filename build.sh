#!/bin/bash
if [ $1 = "laptop" ] ; then
    geant4_dir="/data/opt/geant4/"
    g++ DagSolidMaterial.cpp -I/home/davisa/.local/include/pyne -I/data/opt/dagmc/hdf5/include -I${geant4_dir}/include -I${geant4_dir}/include/Geant4 -L/data/opt/moab_dev/lib -L/data/opt/dagmc/hdf5/lib -L/data/opt/geant4/lib -L/home/davisa/.local/lib -lhdf5 -lpyne -lG4materials -lG4global -lG4clhep -lG4intercoms -o g4mat

    g++ src/DagSolidMaterial.cpp -c -I/data/opt/DagGeant4/include -I/home/davisa/.local/include/pyne -I/data/opt/dagmc/hdf5/include -I${geant4_dir}/include -I${geant4_dir}/include/Geant4 -L/data/opt/moab_dev/lib -L/data/opt/dagmc/hdf5/lib -L/data/opt/geant4/lib -L/home/davisa/.local/lib -lhdf5 -lpyne -lG4materials -lG4global -lG4clhep -lG4intercoms -o g4mat.o

    g++ test.cpp g4mat.o -I/data/opt/DagGeant4/include -I/home/davisa/.local/include/pyne -I/data/opt/dagmc/hdf5/include -I${geant4_dir}/include -I${geant4_dir}/include/Geant4 -L/data/opt/moab_dev/lib -L/data/opt/dagmc/hdf5/lib -L/data/opt/geant4/lib -L/home/davisa/.local/lib -lhdf5 -lpyne -lG4materials -lG4global -lG4clhep -lG4intercoms -o g4mat2


elif [ $1 = "work" ] ; then
    geant4_dir="/mnt/data/opt/geant4/"
    g++ DagSolidMaterial.cpp -I/home/davisa/.local/include/pyne -I/mnt/data/opt/dagmc/hdf5/include -I${geant4_dir}/include -I${geant4_dir}/include/Geant4 -L/mnt/data/opt/moab_dev/lib -L/mnt/data/opt/dagmc/hdf5/lib -L/mnt/data/opt/geant4/lib -L/home/davisa/.local/lib -lhdf5 -lpyne -lG4materials -lG4global -lG4clhep -lG4intercoms -o g4mat
fi
