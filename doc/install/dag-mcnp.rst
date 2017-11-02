..  |DAG-Code| replace:: DAG-MCNP5/6

..  include:: header.txt

..  include:: get_dagmc.txt

..  _mcnp_patch:

Apply the MCNP source code patch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DAGMC can be built with MCNP5_ or MCNP6_ or both, from Los Alamos National
Laboratory. It is export-controlled software so you will need to request it from
RSICC_.

If you are building |DAG-Code|, you will need to copy the MCNP source code from
the DVD into the DAGMC repository and apply a patch so it can be used with
DAGMC. The patch you use must correspond to your version of MCNP. Currently
supported versions of MCNP5 are 5.1.40, 5.1.51, 5.1.60
::

    $ cd src/mcnp/mcnp5
    $ cp -r <path_to_dvd>/MCNP5/Source .
    $ chmod -R u+rw Source
    $ patch -p0 < patch/dagmc.5.1.60.patch

Currently supported versions of MCNP6 are 6_beta2, 6.1, and 6.1.1beta.
::

    $ cd src/mcnp/mcnp6
    $ cp -r <path_to_dvd>/MCNP6/Source .
    $ chmod -R u+rw Source
    $ patch -p0 < patch/dagmc.6.1.1beta.patch

Assuming the patch or patches were succesfully applied, i.e. there were no
warnings or errors, you are now ready to configure DAGMC to produce the desired
build.

..  include:: configure_dag-code_header.txt

Note that all of these examples examples assume that the ``$DATAPATH``
environment variable is set. If it is not set, then the ``-DMCNP5_DATAPATH``
and/or ``-DMCNP6_DATAPATH`` cmake options must be included instead.

**Example 1:** Build the DAGMC interfaces and DAG-MCNP5.
::

    $ cmake .. -DBUILD_MCNP5=ON \
               -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 2:** Build an MPI version of DAG-MCNP5.
::

    $ cmake .. -DBUILD_MCNP5=ON \
               -DMPI_BUILD=ON \
               -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 3:** Build the DAGMC interfaces and DAG-MCNP6.
::

    $ cmake .. -DBUILD_MCNP6=ON \
               -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 4:** Build an MPI version of DAG-MCNP6.
::

    $ cmake .. -DBUILD_MCNP6=ON \
               -DMPI_BUILD=ON \
               -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 5:** Build both DAG-MCNP5 and DAG-MCNP6.
::

    $ cmake .. -DBUILD_MCNP5=ON \
               -DBUILD_MCNP6=ON \
               -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 6:** Build MPI versions of both DAG-MCNP5 and DAG-MCNP6.
::

    $ cmake .. -DBUILD_MCNP5=ON \
               -DBUILD_MCNP6=ON \
               -DMPI_BUILD=ON \
               -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

If the CMake configuration proceeded successfully, you are now ready to install
DAGMC.

..  include:: build_dagmc.txt

..  include:: test_dagmc.txt

..  _MCNP5: https://mcnp.lanl.gov
..  _MCNP6: https://mcnp.lanl.gov
..  _RSICC: https://rsicc.ornl.gov
