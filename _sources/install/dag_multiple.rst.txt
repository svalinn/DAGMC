..  |DAG-Code| replace:: DAGMC with multiple physics codes

..  include:: header.txt

..  include:: get_dagmc.txt

Before configuring DAGMC
~~~~~~~~~~~~~~~~~~~~~~~~

Before you configure DAGMC, you will need to perform the steps listed in the
guides for the individual codes with which you want to install DAGMC.

For DAG-MCNP5/6, you will need to
:ref:`apply one of the MCNP source code patches <mcnp_patch>`.

For FluDAG, you will need to :ref:`install FLUKA <install_fluka>`.

For DAG-Geant4, you will need to :ref:`install Geant4 <install_geant4>`.

..  include:: configure_dag-code_header.txt

The following CMake command will build DAG-Geant4 and FluDAG as well as MPI
versions of DAG-MCNP5 and DAG-MCNP6.
::

    $ cmake .. -DBUILD_MCNP5=ON \
               -DBUILD_MCNP6=ON \
               -DMPI_BUILD=ON \
               -DBUILD_GEANT4=ON \
               -DGEANT4_DIR=$HOME/dagmc_bld/Geant4 \
               -DBUILD_FLUKA=ON \
               -DFLUKA_DIR=$FLUPRO \
               -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

If the CMake configuration proceeded successfully, you are now ready to install
DAGMC.

..  include:: build_dagmc.txt

..  include:: test_dagmc.txt
