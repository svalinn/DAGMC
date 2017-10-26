DAGMC CMake variables
=====================

This page describes the specific CMake variables that affect the DAGMC build.
Note that unless set, all boolean options default to ``OFF``.

    * ``-DBUILD_MCNP5=ON`` Build DAG-MCNP5.

    * ``-DMCNP5_DATAPATH=/path/to/mcnp/data`` Set the path where the DAG-MCNP5
      executable will look for nuclear data. If this variable is not set, it
      will look for the data in the location specified by the ``$DATAPATH``
      environment variable instead.

    * ``-DBUILD_MCNP6=ON`` Build DAG-MCNP6.

    * ``-DMCNP6_DATAPATH=/path/to/mcnp/data`` Set the path where the DAG-MCNP5
      executable will look for nuclear data. If this variable is not set, it
      will look for the data in the location specified by the ``$DATAPATH``
      environment variable instead.

    * ``-DMPI_BUILD=ON`` If building DAG-MCNP5 and/or DAG-MCNP6, build MPI
      versions.

    * ``-DBUILD_FLUKA=ON`` Build FluDAG. If this option is turned on,
      ``-DFLUKA_DIR`` must also be specified.

    * ``-DFLUKA_DIR=/path/to/fluka`` Set the path to Fluka. It should typically
      be set to the ``$FLUPRO`` environment variable.

    * ``-DBUILD_GEANT4=ON`` Build DAG-Geant4 and the DagSolid library. If this
      option is turned on, ``-DGEANT4_DIR`` may need to be specified if the
      CMake files can't find Geant4 by themselves.

    * ``-DGEANT4_DIR=/path/to/geant4`` Set the path to Geant4.

    * ``-DBUILD_TALLY=ON`` Build the DagTally interface. Note that if you are
      using ``-DBUILD_MCNP5=ON`` or ``-DBUILD_MCNP6=ON`` then this is
      automatically set to ON.

    * ``-DBUILD_ALL=ON`` Build DAG-MCNP5, DAG-MCNP6, FluDAG, and DAG-Geant4.

    * ``-DBUILD_STATIC_EXE=ON`` Build static executables.
