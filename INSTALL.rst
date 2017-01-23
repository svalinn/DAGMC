Install Instructions
====================

The new DAGMC toolkit has a full CMake install and build method, for all codes used downstream, even
going so far to replace the MCNP build method with a CMake file.

Prerequisites
===========
In order to install you must have done the following absolute minimum

1) cloned the DAGMC repository, (http://github.com/svalinn/DAGMC)
2) Installed MOAB
3) Installed HDF5
4) Installed Fluka

How To
========
1) Copy the "Source" directory for MCNP5v16 from the LANL/RSICC CD to the mcnp5/ directory in the DAGMC source tree
   cp -r <path to cdrom/MCNP5/Source mcnp5/.
2) Apply the patch from the patch folder
   patch -p1 < patch/dagmc.patch.5.1.60
3) Assuming your patch was succesfully applied, i.e. there were no warnings or errors then we can now start building,
   assuming that you are in the base level of the DAGMC repo, create a new directory and navigate to it.
   mkdir bld
   cd bld
4) We can now configure DAGMC for building, you must include the CMAKE_INSTALL_PREFIX option with a folder where
   you would like the toolkit to be installed, this directory need not exist.

   For example we could just build the DAGMC interfaces and DAG-MCNP5
   cmake ../. -DMOAB_DIR=$MOAB_PATH/lib -DBUILD_MCNP5=ON -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
   or to building MCNP5 in parallel
   cmake ../. -DMOAB_DIR=$MOAB_PATH/lib -DBUILD_MCNP5=ON -DMPI_BUILD=ON \
   -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
   We could build Fluka as well
   cmake ../. -DMOAB_DIR=$MOAB_PATH/lib -DBUILD_MCNP5=ON -DMPI_BUILD=ON \
     -DBUILD_FLUKA=ON -DFLUKA_DIR=$FLUPRO -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
   or only Fluka
   cmake ../. -DMOAB_DIR=$MOAB_PATH/lib -DBUILD_FLUKA=ON -DFLUKA_DIR=$FLUPRO \
             -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
   and or Geant4
   cmake ../. -DMOAB_DIR=$MOAB_PATH/lib -DBUILD_MCNP5=ON -DMPI_BUILD=ON \
          -DBUILD_FLUKA=ON -DFLUKA_DIR=$FLUPRO -DBUILD_GEANT4=ON -DGEANT4_DIR=/mnt/data/opt/geant4  \
          -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

5) Assuming that the cmake step was succesful, i.e. no errors were reported then we can make by issuing the make command
   make
6) If there were no errors, then we can install the DAGMC suite of libraries and tools by issing the install command
   make install
7) If everything was successful, you may have the mcnp5 and mainfludag executables in the bin folder, the libraries in lib
   and the header files in the include folder


  


