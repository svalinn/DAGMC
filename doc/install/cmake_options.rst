DAGMC CMake variables
=====================

This page describes the specific CMake variables that affect the DAGMC build.
Note that unless set, all boolean options default to ``OFF``.

    * ``-DBUILD_MCNP5=ON`` Build DAG-MCNP5. (Default: OFF)

    * ``-DBUILD_MCNP5_PLOT=ON`` Enable DAG-MCNP5 plotting capability. (Default:
      OFF)

    * ``-DMCNP5_DATAPATH=/path/to/mcnp/data`` Set the path where the DAG-MCNP5
      executable will look for nuclear data. If this variable is not set, it
      will look for the data in the location specified by the ``$DATAPATH``
      environment variable instead.

    * ``-DBUILD_MCNP6=ON`` Build DAG-MCNP6. (Default: OFF)

    * ``-DBUILD_MCNP6_PLOT=ON`` Enable DAG-MCNP6 plotting capability. (Default:
      OFF)

    * ``-DMCNP6_DATAPATH=/path/to/mcnp/data`` Set the path where the DAG-MCNP5
      executable will look for nuclear data. If this variable is not set, it
      will look for the data in the location specified by the ``$DATAPATH``
      environment variable instead.

    * ``-DMPI_BUILD=ON`` If building DAG-MCNP5 and/or DAG-MCNP6, build MPI
      versions. (Default: OFF)

    * ``-DBUILD_FLUKA=ON`` Build FluDAG and the fludag library. If this option
      is turned on, ``-DFLUKA_DIR`` must also be specified. (Default: OFF)

    * ``-DFLUKA_DIR=/path/to/fluka`` Set the path to Fluka. It should typically
      be set to the ``$FLUPRO`` environment variable.

    * ``-DBUILD_GEANT4=ON`` Build DAG-Geant4 and the DagSolid library. If this
      option is turned on, ``-DGEANT4_DIR`` must be specified. (Default: OFF)

    * ``-DGEANT4_DIR=/path/to/geant4`` Set the path to Geant4.

    * ``-DBUILD_TALLY=ON`` Build the DagTally interface. (Default: ON)

    * ``-DBUILD_ASTYLE=ON`` Build the Astyle code formatter. (Default: ON)

    * ``-DBUILD_BUILD_OBB=ON`` Build the build_obb tool. (Default: ON)

    * ``-DBUILD_MAKE_WATERTIGHT=ON`` Build the make_watertight tool. (Default:
      ON)

    * ``-DBUILD_TESTS=ON`` Build unit tests where appropriate. (Default: ON)

    * ``-DBUILD_STATIC_EXE=ON`` Build static executables. (Default: OFF)

    * ``-DBUILD_PIC=ON`` Build with position-independent code. (Default: OFF)
