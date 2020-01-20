.. |DAG-Code| replace:: DAGMC
..  _OpenMC: https://mit-crpg.github.io/openmc

Installing for use with OpenMC
==============================

**Note: DagMC can simultaneously be built with support for other physics codes
while being installed as a dependency for OpenMC.**

This document explains how to install DAGMC for use with OpenMC, assuming you have
already installed the required `dependencies <dependencies.html>`_.

DAGMC is an optional dependency of OpenMC_. Therefore, the install process for
this code only generates and installs the DAGMC libraries necessary for linkage
with OpenMC.

.. include:: get_dagmc.txt

.. include:: configure_dag-code_header.txt


Installing DAGMC as a dependency of OpenMC
------------------------------------------

From the build directory, run::

    $ cmake ../src -DMOAB_DIR=$HOME/dagmc_bld/MOAB \
                   -DBUILD_TALLY=ON \
                   -DCMAKE_INSTALL_PATH=$INSTALL_PATH

If the CMake configuration proceeded successfully, you are now ready to install
DAGMC.

Use Make to install the DAGMC libraries.
::

    $ make
    $ make install

If the build was successful, the binaries, libraries, header files, and tests
will be installed to the ``bin``, ``lib``, ``include``, and ``tests``
subdirectories of ``$INSTALL_PATH`` respectively.

Finally, ensure that the DagMC library is in your ``$LD_LIBRARY_PATH``
::

   $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/dagmc_bld/dagmc/lib

Installing OpenMC with DagMC geometry enabled
---------------------------------------------

To install OpenMC with support for DagMC please refer to OpenMC's installation
from source instructions `here <OpenMCInstallFromSource_>`_.

..  _OpenMCInstallFromSource: https://openmc.readthedocs.io/en/latest/usersguide/install.html?#installing-from-source

..  include:: test_dagmc.txt
