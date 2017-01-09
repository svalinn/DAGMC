..  |DAG-Code| replace:: DAG-MCNP5

..  include:: header.txt

..  include:: get_dagmc.txt

..  _mcnp5_patch:

Apply the MCNP5 source code patch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DAG-MCNP5 uses `MCNP5 <https://mcnp.lanl.gov>`_ from LANL. It is
export-controlled software so you will need to request it from
`RSICC <https://rsicc.ornl.gov>`_.

If you are building DAG-MCNP5, you need to copy the MCNP5 source code from the
DVD into the DAGMC repository and patch it so it can be used with DAGMC.
::

    $ cd mcnp/mcnp5
    $ cp -r <path_to_dvd>/MCNP5/Source .
    $ patch -p0 < patch/dagmc.patch.5.1.60

Assuming the patch was succesfully applied, i.e. there were no warnings or
errors, you are now ready to configure DAGMC to produce the desired build.

..  include:: configure_dag-code_header.txt

**Example 1:** Build the DAGMC interfaces and DAG-MCNP5, using the
``$DATAPATH`` environment variable to specify the location of the MCNP data.
::

    $ cmake .. -DBUILD_MCNP5=ON \
               -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 2:** Build the DAGMC interfaces and DAG-MCNP5, assuming that the
``$DATAPATH`` environment variable is undefined.
::

    $ cmake .. -DBUILD_MCNP5=ON \
               -DMCNP5_DATAPATH=<path to MCNP data> \
               -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 3:** Build an MPI version of DAG-MCNP5.
::

    $ cmake .. -DBUILD_MCNP5=ON \
               -DMPI_BUILD=ON \
               -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

If the CMake configuration proceeded successfully, you are now ready to install
DAGMC.

..  include:: build_dagmc.txt

..  include:: test_dagmc.txt
