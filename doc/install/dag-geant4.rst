..  |DAG-Code| replace:: DAG-Geant4

..  include:: header.txt

..  include:: get_dagmc.txt

..  _install_geant4:

Install Geant4
~~~~~~~~~~~~~~

DAG-Geant4 uses Geant4_ from CERN. It is open-source software so you do not need
to register an account.

Refer to the `getting started <Geant4_getting_started_>`_ page for information
about downloading and installing Geant4. The following commands can be used to
download the Geant4 source code and set it up for building:
::

    $ cd $HOME/dagmc_bld
    $ mkdir -p Geant4/bld
    $ cd Geant4
    $ wget http://geant4.cern.ch/support/source/geant4.10.04.tar.gz
    $ tar -xzvf geant4.10.04.tar.gz
    $ ln -s geant4.10.04 src

Geant4 uses a CMake build, and we recommend using the following flags when
installing it with the purpose of coupling with DAGMC:
::

    $ cd bld
    $ cmake ../src -DGEANT4_INSTALL_DATA=ON \
                   -DGEANT4_USE_QT=ON \  # or -DGEANT4_USE_OPENGL_X11=ON
                   -DGEANT4_USE_SYSTEM_EXPAT=OFF
    $ make
    $ make install

The following CMake command will build DAG-Geant4, assuming you built Geant4 as
specified in the Geant4 build instructions above.
::

    $ cmake .. -DBUILD_GEANT4=ON \
               -DGEANT4_DIR=$HOME/dagmc_bld/Geant4 \
               -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

If the CMake configuration proceeded successfully, you are now ready to install
DAGMC.

..  include:: build_dagmc.txt

To use DagGeant4, you first must make add the Geant4 library directory to your
``$LD_LIBRARY_PATH``.
::

    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/dagmc_bld/Geant4/lib

Test |DAG-Code|
~~~~~~~~~~~~~~~

..  include:: configure_dag-code_header.txt

To run the DagSolid unit tests, use the following command. Make sure that the
Geant4 library directory is in your ``$LD_LIBRARY_PATH`` as specified above.
::

    $ cd $INSTALL_PATH/tests
    $ ./dagsolid_unit_tests

If the tests ran successfully, the last few lines of the screen output will look
like this:
::

    [       OK ] DagSolidTest.surface_area_test (5 ms)
    [----------] 16 tests from DagSolidTest (228 ms total)

    [----------] Global test environment tear-down
    [==========] 16 tests from 1 test case ran. (228 ms total)
    [  PASSED  ] 16 tests.

..  include:: test_dagmc.txt

..  _Geant4: http://geant4.cern.ch
..  _Geant4_getting_started: http://geant4.cern.ch/support/gettingstarted.shtml
