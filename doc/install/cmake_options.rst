DAGMC CMake variables
=====================

This page describes the specific CMake variables that affect the DAGMC build.
Note that unless set, all boolean options default to ``OFF``.

    * ``-DMOAB_ROOT=/path/to/moab`` Path to MOAB.

    * ``-DBUILD_MCNP5=ON`` Build DAG-MCNP5. (Default: OFF)

    * ``-DBUILD_MCNP6=ON`` Build DAG-MCNP6. (Default: OFF)

    * ``-DMCNP5_DATAPATH=/path/to/mcnp/data`` Set the path where the DAG-MCNP5
      executable will look for nuclear data. If this variable is not set, it
      will look for the data in the location specified by the ``$DATAPATH``
      environment variable instead.

    * ``-DMCNP6_DATAPATH=/path/to/mcnp/data`` Set the path where the DAG-MCNP5
      executable will look for nuclear data. If this variable is not set, it
      will look for the data in the location specified by the ``$DATAPATH``
      environment variable instead.

    * ``-DBUILD_MCNP_PLOT=ON`` Enable DAG-MCNP5/6 plotting capability. (Default:
      OFF)

    * ``-DBUILD_MCNP_MPI=ON`` If building DAG-MCNP5/6, build MPI versions.
      (Default: OFF)

    * ``-DBUILD_MCNP_PYNE_SOURCE=ON`` If building DAG-MCNP5/6, build with PyNE
      mesh source routine support. (Default: OFF)

    * ``-DBUILD_FLUKA=ON`` Build FluDAG and the fludag library. If this option
      is turned on, ``-DFLUKA_DIR`` must also be specified. (Default: OFF)

    * ``-DFLUKA_DIR=/path/to/fluka`` Set the path to Fluka. It should typically
      be set to the ``$FLUPRO`` environment variable.

    * ``-DBUILD_GEANT4=ON`` Build DAG-Geant4 and the DagSolid library. If this
      option is turned on, ``-DGEANT4_DIR`` must be specified. (Default: OFF)

    * ``-DGEANT4_DIR=/path/to/geant4`` Set the path to Geant4.

    * ``-DBUILD_TALLY=ON`` Build the DagTally interface. (Default: ON)

    * ``-DBUILD_BUILD_OBB=ON`` Build the build_obb tool. (Default: ON)

    * ``-DBUILD_MAKE_WATERTIGHT=ON`` Build the make_watertight tool. (Default:
      ON)

    * ``-DBUILD_TESTS=ON`` Build unit tests where appropriate. (Default: ON)

    * ``-DBUILD_CI_TESTS=ON`` Build everything needed to run the continuous
      integration tests. (Default: OFF)

    * ``-DBUILD_SHARED_LIBS=ON`` Build shared libraries. (Default: ON)

    * ``-DBUILD_STATIC_LIBS=ON`` Build static libraries. (Default: ON)

    * ``-DBUILD_STATIC_EXE=ON`` Build static executables. (Default: OFF)

    * ``-DBUILD_PIC=ON`` Build with position-independent code. (Default: OFF)

    * ``-DBUILD_RPATH=ON`` Build with RPATH functionality. (Default: ON)
