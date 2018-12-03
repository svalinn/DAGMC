.. |DAG-Code| replace:: OpenMC

.. include:: header.txt

DAGMC is a dependency of OpenMC. Therefore, the install process for this code
only generates and installs the DAGMC libraries to be used during compilation of
OpenMC.

.. include:: get_dagmc.txt

.. include:: configure_dag-code-header.txt

Installing DAGMC as a dependency of OpenMC
------------------------------------------

From the build directory, run::

    $ cmake .. -DMOAB_DIR=$HOME/dagmc_bld/MOAB \
               -DBUILD_TALLY=ON \
               -DCMAKE_INSTALL_PATH=$INSTALL_PATH

If the CMake configuration proceeded successfully, you are now ready to install
DAGMC.

..  include:: build_dagmc.txt

..  include:: test_dagmc.txt
