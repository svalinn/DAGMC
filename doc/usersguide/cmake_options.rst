CMAKE Install options
----------------------------------------
There are now large number of install options available to the person
building DAGMC; we use the CMake system, which builds custom makefiles
for your particular system.

The options follow.  Note that, unless set, boolean options default to ```OFF```:

 * ```-DBUILD_MCNP5``` controls whether to build the DAG-MCNP5 executable:  allowed options are ```ON``` and ```OFF```.
 * ```-DBUILD_FLUKA``` controls whether to build the FluDAG exectuable:  allowed
   options are ```ON``` and ```OFF```.  If you would like to build the FluDAG interface, 
   you must also set the ```-DFLUKA_DIR``` variable.
 * ```-DFLUKA_DIR``` sets the path to the Fluka libraries and scripts.  This variable must be set 
   in order to build the FluDAG executable.  It typically should be set to the $FLUPRO environment variable
   but if you have not got it set, then you should set it to the path to the Fluka directory.
 * ```-DBUILD_GEANT4``` controls whether to build the DagSolid libraries and the DagGeant4 
   exectuable.  Options are ```ON``` and ```OFF```; you may also need to set the ```-DGEANT4_DIR``` variable if
   the file ```DAGMC/cmake/FindGeant4.cmake``` fails to find Geant4.
 * ```-DGEANT4_DIR``` sets a hint to the FindGeant4.cmake file to where your Geant4 installation is.  If the FindGeant4.cmake
   file fails to locate Geant4 you should set this to variable to the path to the Geant4 installation directory.
 * ```-DBUILD_TALLY``` sets whether build the DagTally interface:  allowed options are ```ON``` and ```OFF```.  
   Note that if you have ```-DBUILD_MCNP5=ON``` then ```-DBUILD_TALLY=ON``` is automatically set.
 * ```-DBUILD_ALL``` sets whether to build all interfaces:  allowed options are ```ON``` and ```OFF```.
 * ```-DSTATIC_LIB``` sets whether to make a shared or static build of the Dag libraries.  If you would like them 
 static then set it to ```ON``; this variable defaults to ```OFF```.
