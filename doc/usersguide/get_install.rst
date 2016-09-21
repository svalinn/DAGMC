Install Guide
=============

These instructions explain how to install DAGMC and its dependencies.

Before you begin
~~~~~~~~~~~~~~~~

This guide assumes that you are building DAGMC in the subdirectory ``dagmc_bld``
of your home directory; i.e.
::
    $ cd $HOME
    $ mkdir dagmc_bld
    $ cd dagmc_bld

Prerequisites
~~~~~~~~~~~~~

 * `LAPACK <http://www.netlib.org/lapack>`_
 * `HDF5 <http://www.hdfgroup.org/HDF5>`_
 * `MOAB <http://sigma.mcs.anl.gov/moab-library>`_

LAPACK
------

Debian linux users can install LAPACK with:
::
    $ sudo apt-get install libblas-dev liblapack-dev

Redhat linux users can do likewise with:
::
    $ sudo yum install libblas-dev liblapack-dev

Refer to the `LAPACK website <http://www.netlib.org/lapack>`_ if you need to
build from source.

HDF5
------

Debian linux users can install the latest HDF5 release with:
::
    $ sudo apt-get install hdf5-dev

Redhat linux users can do likewise with:
::
    $ sudo yum install hdf5-dev

You can also elect to install HDF5 from source. The tarball containing the HDF5
source code can also be downloaded from the
`HDF5 website <https://support.hdfgroup.org/HDF5/>`_.
Note that if you choose this option, we recommend you obtain HDF5 version 1.8.13
instead of the newest version. The following commands can be used to install
HDF5 from source.
::
    $ cd $HOME/dagmc_bld
    $ mkdir -p HDF5/bld
    $ cd HDF5
    $ wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.13/src/hdf5-1.8.13.tar.gz
    $ tar -xzvf hdf5-1.8.13.tar.gz
    $ ln -s hdf5-1.8.13 src
    $ cd bld
    $ ../src/configure --enable-shared \
                       --prefix=$HOME/dagmc_bld/HDF5
    $ make
    $ make install

MOAB
------

As of DAGMC version 2.0, MOAB version 4.9.2 or higher is required. The following
commands can be used to download MOAB from its `source repository
<https://bitbucket.org/fathomteam/moab>`_ and set it up for building.
::
    $ cd $HOME/dagmc_bld
    $ mkdir -p MOAB/bld
    $ cd MOAB
    $ git clone https://bitbucket.org/fathomteam/moab
    $ cd moab
    $ git checkout master
    $ autoreconf -fi
    $ cd ..
    $ ln -s moab src

The command ``git checkout master`` is redundant but is included here for
completeness. If you would prefer to use a specific version instead of the
master branch, e.g. version 4.9.2, use ``git checkout Version4.9.2`` instead.

The following commands should be used to build MOAB.
Note that the ``--enable-dagmc`` configure option is required.
::
    $ cd bld
    $ ../src/configure --enable-dagmc \
                       --enable-optimize \
                       --enable-shared \
                       --disable-debug \
                       --with-hdf5=$HOME/dagmc_bld/HDF5 \
                       --prefix=$HOME/dagmc_bld/MOAB
    $ make
    $ make install

Physics codes
~~~~~~~~~~~~~

MCNP5
------

DAG-MCNP5 uses `MCNP5 <https://mcnp.lanl.gov>`_ from LANL. It is
export-controlled software so you will need to request it from
`RSICC <https://rsicc.ornl.gov>`_.

Nothing needs to be done with MCNP5 prior to installing DAGMC.

Geant4
------

DAG-Geant4 uses `Geant4 <http://geant4.cern.ch>`_ from CERN. It is open-source
software so you do not need to register an account.

Refer to the
`getting started <http://geant4.cern.ch/support/gettingstarted.shtml>`_ page for
information about downloading and installing Geant4. The following commands can
be used to download the Geant4 source code and set it up for building:
::
    $ cd $HOME/dagmc_bld
    $ mkdir -p Geant4/bld
    $ cd Geant4
    $ wget http://geant4.cern.ch/support/source/geant4.10.02.p02.tar.gz
    $ tar -xzvf geant4.10.02.p02.tar.gz
    $ ln -s geant4.10.02.p02 src

Geant4 uses a CMake build, and we recommend using the following flags when
installing it with the purpose of coupling with DAGMC:
::
    $ cd bld
    $ cmake ../src -DGEANT4_INSTALL_DATA=ON \
                   -DGEANT4_USE_QT=ON \  # or -DGEANT4_USE_OPENGL_X11=ON
                   -DGEANT4_USE_SYSTEM_EXPAT=OFF
    $ make
    $ make install

FLUKA
------

FluDAG uses `FLUKA <http://www.fluka.org/fluka.php>`_ from CERN/INFN. In order
to download FLUKA you need to become a registered user, which you can do at the
`FLUKA register <https://www.fluka.org/fluka.php?id=secured_intro>`_ page.

Save your user ID and password for future FLUKA updates. We recommend an x64
worfklow and thus you should download the 64-bit executable. The name of the
downloaded tarball is of the form ``fluka20xx.xx-linux-gfor64bitAA.tar.gz``.
Refer to the
`installation instructions <http://www.fluka.org/fluka.php?id=ins_run&mm2=3>`_
when building FLUKA.

Take care to follow the FLUKA site instructions when setting the
``$FLUPRO`` and ``$FLUFOR`` environment variables.

Tripoli4
--------

DAG-Tripoli4 uses Tripoli4, which is is distributed by CEA/Saclay as a binary
executable. For access to DAG-Tripoli4, please contact `Jean-Christophe Trama
<mailto:jean-christophe.trama@cea.fr>`_.

Environment variables
~~~~~~~~~~~~~~~~~~~~~

After installing HDF5 and MOAB, you need to make sure the system can find them
when it comes time to build DAGMC. This is done by adding some directories to
your ``$PATH`` and ``$LD_LIBRARY_PATH``.
::
    $ export PATH=$PATH:$HOME/.local/bin: \
                        $HOME/dagmc_bld/HDF5/bin: \
                        $HOME/dagmc_bld/MOAB/bin
    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/.local/lib: \
                                              $HOME/dagmc_bld/HDF5/lib: \
                                              $HOME/dagmc_bld/MOAB/lib

You can use the following commands to test whether HDF5 and MOAB were built
successfully.
::
    $ which h5ls
    $ which mbconvert

If you installed Geant4, you will also need to add the Geant4 directories to
your ``$PATH`` and ``$LD_LIBRARY_PATH``.
::
    $ export PATH=$PATH:$HOME/dagmc_bld/Geant4/bin
    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/dagmc_bld/Geant4/lib

DAGMC build procedure
~~~~~~~~~~~~~~~~~~~~~

Get DAGMC
---------

The first step is to clone the `DAGMC repository
<https://github.com/svalinn/DAGMC>`_.
::
    $ cd $HOME/dagmc_bld
    $ mkdir DAGMC
    $ cd DAGMC
    $ git clone https://github.com/svalinn/DAGMC
    $ cd DAGMC
    $ git checkout develop

If you are building DAG-MCNP5, you need to copy the MCNP5 source code from the
DVD into the DAGMC repository and patch it so it can be used with DAGMC.
::
    $ cd mcnp/mcnp5
    $ cp -r <path_to_dvd>/MCNP5/Source .
    $ patch -p0 < patch/dagmc.patch.5.1.60

If you are building FluDAG, you will need to patch FLUKA's run script ``rfluka``
in order to allow for some DAGMC-specific options.
::
    $ cd $FLUPRO/flutil
    $ patch -Nb rfluka $HOME/dagmc_bld/DAGMC/fluka/rfluka.patch

Assuming the patch was succesfully applied, i.e. there were no warnings or
errors, you can now configure DAGMC to produce the desired build.

Configure DAGMC
-------------

CMake variables are used to configure DAGMC with your desired build options. A
few examples will be shown here, but you can see a list of all possible options
`here <cmake_options.html>`_.

First, create and enter the build directory.
::
    $ cd $HOME/dagmc_bld/DAGMC
    $ mkdir bld
    $ cd bld

Then, choose where you want to install DAGMC. This is where the binaries,
libraries, header files, etc. will be placed. This guide uses ``$INSTALL_PATH``
to represent this location.
::
    $ INSTALL_PATH=$HOME/dagmc_bld/DAGMC

**Example 1:** Build the DAGMC interfaces and DAG-MCNP5, using the
``$DATAPATH`` environment variable to specify the location of the MCNP data.
::
    $ cmake ../src -DBUILD_MCNP5=ON \
                   -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 2:** Build the DAGMC interfaces and DAG-MCNP5, assuming that the
``$DATAPATH`` environment variable is undefined.
::
    $ cmake ../src -DBUILD_MCNP5=ON \
                   -DMCNP5_DATAPATH=<path to MCNP data> \
                   -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 3:** Build an MPI version of DAG-MCNP5.
::
    $ cmake ../src -DBUILD_MCNP5=ON \
                   -DMPI_BUILD=ON \
                   -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 4:** Build DAG-Geant4 (assuming you built Geant4 as specified in the
Geant4 build instructions above).
::
    $ cmake ../src -DBUILD_GEANT4=ON \
                   -DGEANT4_DIR=$HOME/dagmc_bld/Geant4 \
                   -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 5:** Build FluDAG. Note that ``$FLUPRO`` should have previously been
defined as part of the FLUKA install.
::
    $ cmake ../src -DBUILD_FLUKA=ON \
                   -DFLUKA_DIR=$FLUPRO \
                   -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 6:** Build an MPI version of DAG-MCNP5 as well as DAG-Geant4 and
FluDAG.
::
    $ cmake ../src -DBUILD_MCNP5=ON \
                   -DMPI_BUILD=ON \
                   -DBUILD_GEANT4=ON \
                   -DGEANT4_DIR=$HOME/dagmc_bld/Geant4 \
                   -DBUILD_FLUKA=ON \
                   -DFLUKA_DIR=$FLUPRO \
                   -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

If the CMake configuration proceeded successfully, you can now install DAGMC.

Install DAGMC
------------

Use Make to install DAGMC.
::
    $ make
    $ make install

If the build was successful, the binaries, libraries, header files, and tests
will be installed to the ``bin``, ``lib``, ``include``, and ``tests``
subdirectories of ``$INSTALL_PATH`` respectively.

In order to use DAGMC, make these additions to your paths:
::
    $ export PATH=$PATH:$INSTALL_PATH/bin
    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INSTALL_PATH/lib

Test DAGMC
----------

We regularly run the DAGMC test suite on
`Travis <https://travis-ci.org/svalinn/DAGMC>`_ as part of our continuous
integration system. You may however wish to run the tests yourself in order to
verify you have installed DAGMC correctly.

To run the FluDAG unit tests, use
::
    $ cd $INSTALL_PATH/tests
    $ ./fludag_unit_tests

If the tests ran successfully, the last few lines of the screen output will look
like this:
::
    [       OK ] FluDAGTest.GFireGoodPropStep (5 ms)
    [----------] 3 tests from FluDAGTest (108 ms total)

    [----------] Global test environment tear-down
    [==========] 3 tests from 1 test case ran. (108 ms total)
    [  PASSED  ] 3 tests.

To run the DagSolid unit tests use the following command. Make sure that the
Geant4 directories are in your ``$PATH`` and ``$LD_LIBRARY_PATH`` as specified
in the `Environment variables`_ section.
::
    $ cd $INSTALL_PATH/tests
    $ ./dagsolid_unit_tests

Again, with successful execution the last few lines of the screen output will
look like this:
::
    [       OK ] DagSolidTest.surface_area_test (5 ms)
    [----------] 16 tests from DagSolidTest (228 ms total)

    [----------] Global test environment tear-down
    [==========] 16 tests from 1 test case ran. (228 ms total)
    [  PASSED  ] 16 tests.

If the tests have completed successfully, you are now ready to run your first
DAGMC problem. See the `DAGMC workflow guides <workflow.html>`_ for more
information.

.. toctree::
   :maxdepth: 1

   cmake_options
