DAGMC Toolkit
----------------------------------------
Getting and Installing the Toolkit

DAGMC is a toolkit that provides direct geometry support to any Monte
Carlo radiation transport code.  The primary code used for development
and testing of DAGMC is `MCNP5 <http://laws.lanl.gov/vhosts/mcnp.lanl.gov/mcnp5.shtml>`_,
developed by `Los Alamos National Laboratory <http://www.lanl.gov>`_
and distributed by the `Radiation Safety Information Computing Center
<http://rsicc.ornl.gov>`_.  There has also been experience with MCNPX
(LANL), Tripoli4 (CEA/Saclay), `FLUKA <http://www.fluka.org>`_ (CERN/INFN), 
and `Geant4 <http://www.geant4.cern.ch/>`_ (CERN).

These instructions describe the basic steps for downloading and
installing the software libraries that provide the DAGMC toolkit for
integration with Monte Carlo codes.  After this, code-specific
instructions will be give for each code. 

Toolkit Installation
++++++++++++++++++++++++++++

This section details the installation and build steps for the prerequisite packages for the the DAGMC
toolkit with specific physics codes.

Requirements
~~~~~~~~~~~~~

In order to install you must have done the following:

 * Cloned the `DAGMC <http://github.com/svalinn/DAGMC>`_ repository
 * Installed `HDF5 <http://www.hdfgroup.org/HDF5/>`_
 * Install `CGM <http://trac.mcs.anl.gov/projects/ITAPS/wiki/CGM>`_, using the --with-cubit option
 * Installed `MOAB <http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_,
   using options --with-cgm --with-hdf5 --with-dagmc --without-netcdf 
   If you need to prepare meshed geometries the following are also required
   a) Install `CUBIT <http://cubit.sandia.gov>`_ v12.2 or v13.1
 * Installed Lapack
   Note that the MCNP build automatically builds the dagtally library, which uses Lapack 
 * Installed PyNE
 * Installed `FLUKA <http://www.fluka.org>`_ - and/or - 
 * Installed Geant4


Building
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With installation of the DAGMC Toolkit, the dependency stack will look like this:

* Some physics package, e.g. MCNP5
   * `MOAB/DAGMC <http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_
   * `PyNE <http://pyne.io/install.html>`_
   * `HDF5 <http://www.hdfgroup.org/HDF5/release/obtain5.html>`_
   * `CGM <http://trac.mcs.anl.gov/projects/ITAPS/wiki/CGM>`_ 
       * ACIS v19, or `CUBIT <http://cubit.sandia.gov>`_ v12.2 or v13.1 (with CGM `trunk <http://ftp.mcs.anl.gov/pub/fathom/cgm-nightly-trunk.tar.gz>`_ only)


Assumptions and conventions that are used in these instructions:

* all operations are carried out in the a subdirectory ``dagmc_bld`` of a user's home directory
* the path to CUBIT files is known, e.g. ``/path/to/cubit``.  This is the directory that contains the Python script file ``cubit`` and a ``bin`` subdirectory.  
* all tarballs reside in user's home directory

If these do not apply to you, please modify your steps accordingly.

     *(For a shortcut to installing DAG-MCNP5.1.60 you may be able to use the build_dagmc_stack.bash script .)*

CGM
=====

Create a CGM directory to build in:
::
    prompt%> mkdir -p $HOME/dagmc_bld/CGM/bld
    prompt%> cd $HOME/dagmc_bld/CGM

If you are building CGM with Cubit v13.1 or v12.2 you *must* clone the cgm repository using, respectively, one of the following methods:
::
    prompt%> git clone https://bitbucket.org/fathomteam/cgm/
    prompt%> cd cgm
    prompt%> git checkout 13.1.1
    prompt%> autoreconf -fi
    prompt%> cd ..
    prompt%> ln -s cgm src

-or-
::
    prompt%> git clone https://bitbucket.org/fathomteam/cgm/
    prompt%> cd cgm
    prompt%> git checkout 12.2.2
    prompt%> autoreconf -fi
    prompt%> cd ..
    prompt%> ln -s cgm src


Configure and build CGM:
::
    prompt%> cd bld
    prompt%> ../src/configure --enable-optimize \
              --enable-shared --disable-debug \
              --with-cubit=/path/to/cubit  \
              --prefix=$HOME/dagmc_bld/CGM
    prompt%> make
    prompt%> make install


HDF5
======

Debian users may conveniently install the latest HDF5 release with the command:
::
    prompt%> sudo apt-get install hdf5-dev

Fedora users can do likewise with this command:
::
    prompt%> sudo yum install hdf5-dev

The HDF5 tarball can also be downloaded from the HDF5 `website <http://www.hdfgroup.org/HDF5/release/obtain5.html>`_.  On a Linux machine the wget command may be used to get the most recent release, which is currently hdf5-1.8.13:
::
    prompt%> wget \
    http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.13/src/hdf5-1.8.13.tar.gz

See the `HDF5 ftp site <http://www.hdfgroup.org/ftp/HDF5/releases>`_ for available versions.

In the case of a tarball, create a directory and install HDF5:
::
    prompt%> mkdir -p $HOME/dagmc_bld/HDF5/bld
    prompt%> cd $HOME/dagmc_bld/HDF5
    prompt%> tar xzf ~/hdf5-1.8.11.tar.gz
    prompt%> ln -s hdf5-1.8.11 src
    prompt%> cd bld
    prompt%> ../src/configure --enable-shared --prefix=$HOME/dagmc_bld/HDF5
    prompt%> make
    prompt%> make install


MOAB
======

Note:  MOAB version 4.7.0 is the earliest version that may be used.

Create a MOAB directory to install in
::
    prompt%> mkdir -p $HOME/dagmc_bld/MOAB/bld
    prompt%> cd $HOME/dagmc_bld/MOAB


If installing MOAB from the git repository:
::
    prompt%> git clone https://bitbucket.org/fathomteam/moab/
    prompt%> cd moab
    prompt%> git checkout Version4.7.0
    prompt%> autoreconf -fi
    prompt%> cd ..
    prompt%> ln -s trunk src

In all MOAB cases:
::
    prompt%> cd bld
    prompt%> ../src/configure --enable-optimize \
              --enable-shared --disable-debug \
              --with-cgm=$HOME/dagmc_bld/CGM  \
              --with-hdf5=$HOME/dagmc_bld/HDF5 \
	      --enable-dagmc \
              --prefix=$HOME/dagmc_bld/MOAB
    prompt%> make
    prompt%> make install


PyNE
=====
`PyNE <http://pyne.io>`_ is a Python-based nuclear materials data handling package.  
Integration of the DAGMC Toolkit with any physics package, e.g.  FLUKA (FluDAG) or 
Geant4 (DAGSolid), now requires that this library be installed.  Directions for 
installing PyNE are `here <http://pyne.io/install.html>`_.  We recommend building the 
dependencies individually rather than using the Conda Build method.

Post Install
~~~~~~~~~~~~~~

Having installed all the prerequisite tools, namely Cubit, CGM, HDF5, MOAB and PyNE, the user
must ensure that the system has access to the libraries and programs that have been built.
Therefore modify the $PATH and $LD_LIBRARY_PATH environments accordingly:
:: 

    prompt%> export PATH=$PATH:$HOME/dagmc_bld/path/to/cubit/bin: \
                               $HOME/.local/bin: \
                               $HOME/dagmc_bld/HDF5/bin: \
                               $HOME/dagmc_bld/MOAB/bin
    prompt%> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH: \
                               $HOME/dagmc_bld/path/to/cubit/bin: \
                               $HOME/.local/lib: \
                               $HOME/dagmc_bld/HDF5/lib: \
                               $HOME/dagmc_bld/MOAB/lib:/$HOME/dagmc_bld/CGM/lib
 

Toolkit Applications
+++++++++++++++++++++++++++++++++++++++++++++

Install DAGMC
~~~~~~~~~~~~~

Clone the DAGMC repository
::
    prompt%> cd $HOME/dagmc_bld
    prompt%> git clone https://github.com/svalinn/DAGMC
    prompt%> cd DAGMC
    prompt%> git checkout develop


*(A version of the build instructions, INSTALL.rst, is in the DAGMC directory)*.

Install FLUKA
~~~~~~~~~~~~~~
FluDAG uses `FLUKA <http://www.fluka.org>`_ from CERN/INFN with the DAGMC Toolkit.

In order to download FLUKA you need to become a registered user, which you can do at 
the `FLUKA register <https://www.fluka.org/fluka.php?id=secured_intro>`_ page from a link on the main FLUKA page.
Save the user id and password for future FLUKA updates.  We recommend an x64 worfklow and as such you should download
the 64 bit executable.  The download name is of the form *fluka20xx.xx-linux-gfor64bitAA.tar.gz*.
See the `site <http://www.fluka.org/fluka.php?id=ins_run&mm2=3>`_ for instructions.

Once the FLUPRO environment variables have been set, confirm that you have a working install of Fluka and proceed to
the next steps.  

Install Geant4
~~~~~~~~~~~~~~~~
`Geant4 <http://geant4.cern.ch>`_, a toolkit for the simulation of the passage of particles through matter, can be found 
`here <http://geant4.cern.ch/support/gettingstarted.shtml>`_,  including a link to instructions for installation.


Build DAGMC Interfaces
~~~~~~~~~~~~~~~~~~~~~~

The DAGMC toolkit now has a full CMake install and build method for all codes used downstream.  It even
replaces the MCNP build method with a CMake file.

Note that in addition to the detailed instructions above for building the MOAB stack, you must also install
Lapack using your favorite method, for example, "sudo apt-get install liblapack-dev libblas-dev".

Populate and Patch 
============================================
Populate the mcnp5 subdirectory of DAGMC and apply the dagmc patch.

Copy the "Source" directory for MCNP5v16 from the LANL/RSICC CD to the mcnp5/ directory in the DAGMC source tree
::
    prompt%> cd $HOME/dagmc_bld
    prompt%> mkdir -p $HOME/damc_bld/mcnp5
    prompt%> cp -r <path to cdrom/MCNP5/Source mcnp5/

Apply the patch from the patch folder
::
    prompt%> patch -p1 < patch/dagmc.patch.5.1.60


Configure 
===================

Assuming the patch was succesfully applied, i.e. there were no warnings or errors, then we can now 
configure the DAGMC cmake system for the desired build.  From the base level of the DAGMC repo 
create a build directory and navigate to it.
::
    prompt%> cd $HOME/dagmc_bld/DAGMC
    prompt%> mkdir bld
    prompt%> cd bld


We can now configure DAGMC for building.  The CMake system can be configured to build any 
or all of the following
   * MCNP5 with or without MPI
   * GEANT4 (DAGSolid)
   * FLUKA  (FluDAG) 
   
You will need to include the CMAKE_INSTALL_PREFIX=install_dir option as part of the configuration.  When the 
build command 'make install' is invoked, libraries, executables, tests, and include files are installed in 
subdirectories under install_dir.  It is common to use -DCMAKE_INSTALL_PREFIX=..', which creates and populates 
these directories one level above the build directory, that is, in the DAGMC directory.  
Note the '-D' in front of CMAKE_INSTALL_PREFIX, and all of the configuration variables, defines the variable
for the cmake system.

For the examples that follow, it is assumed you are in the bld directory of DAGMC:
::
    prompt%> cd $HOME/dagmc_bld/DAGMC/bld

In the examples, the environment variable, "INSTALL_PATH", can point to any location
where you want the libraries ($INSTALL_PATH/lib), executables ($INSTALL_PATH/bin), and
other build products to be installed.  It is typically set to the DAGMC directory, i.e.
::
    prompt%> export INSTALL_PATH=$HOME/dagmc_bld/DAGMC


**Example 1:**  Build the DAGMC interfaces and DAG-MCNP5
::
    prompt%> cmake ../. -DBUILD_MCNP5=ON -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH


DAG-MCNP5 with MCNP5 in parallel:

**Example 2:**  Build MCNP5 in parallel
::
    prompt%> cmake ../. -DBUILD_MCNP5=ON -DMPI_BUILD=ON \
                        -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH



**Example 3:**  Build MCNP5 in parallel and build the dagmc-enabled FLUKA 
Note that $FLUPRO should have been previously defined as part of the FLUKA install.
::
    prompt%> cmake ../. -DBUILD_MCNP5=ON -DMPI_BUILD=ON \
                        -DBUILD_FLUKA=ON -DFLUKA_DIR=$FLUPRO \
			-DCMAKE_INSTALL_PREFIX=$INSTALL_PATH


**Example 4:** Build only FluDAG
::
    prompt%> cmake ../. -DBUILD_FLUKA=ON -DFLUKA_DIR=$FLUPRO \
                        -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH



**Example 5:**  Build MCNP, FluDAG and Geant4-enabled DAGMC
::
    prompt%> cmake ../. -DBUILD_MCNP5=ON  -DMPI_BUILD=ON \
                        -DBUILD_FLUKA=ON  -DFLUKA_DIR=$FLUPRO \
			-DBUILD_GEANT4=ON -DGEANT4_DIR=path/to/geant4 \
                        -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

Compile and Install
~~~~~~~~~~~~~~~~~~~~~

Assuming that the cmake step was succesful, i.e. no errors were reported, compile by issuing the make command
::
    prompt%> make

If there were no errors, install the DAGMC suite of libraries and tools by issuing the install command
::
    prompt%> make install

If everything was successful, you may have the mcnp5 and mainfludag executables in the $INSTALL_PATH/bin folder, 
the libraries in $INSTALL_PATH/lib and the header files in the $INSTALL_PATH/include folder

Test
~~~~

You may wish to run the tests in the $INSTALL_PATH/tests directory to verify correct installation.  To do this requires
that $INSTALL_PATH/bin be in your PATH and $INSTALL_PATH/lib be in your LD_LIBRARY_PATH:
::
    prompt%> export PATH=$PATH:$INSTALL_PATH/bin
    prompt%> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INSTALL_PATH/lib

Note that this assumes you have previously set the environment variables per the `Post Install`_ section.

With these environment variables you can run fludag_unit_tests:
::
    prompt%> cd $INSTALL_PATH/tests
    prompt%> ./fludag_unit_tests
 
With successful execution the last few lines of the screen output will look similar to:
::
    [       OK ] FluDAGTest.GFireGoodPropStep (5 ms)
    [----------] 3 tests from FluDAGTest (108 ms total)

    [----------] Global test environment tear-down
    [==========] 3 tests from 1 test case ran. (108 ms total)
    [  PASSED  ] 3 tests.

To run dagsolid_unit_test, in addition to the settings just mentioned, you must also execute
a script that was created at the time geant4 was built:
::
    prompt%> source path/to/geant4/bld/geant4make.sh
    prompt%> cd $INSTALL_PATH/tests
    prompt%> ./dagsolid_unit_tests

Again, with successful execution the last few lines of screen output are:
::
    [       OK ] DagSolidTest.surface_area_test (5 ms)
    [----------] 16 tests from DagSolidTest (228 ms total)

    [----------] Global test environment tear-down
    [==========] 16 tests from 1 test case ran. (228 ms total)
    [  PASSED  ] 16 tests.

With testing successfully completed you are now ready to run a `problem <index.html>`_.

DAG-Tripoli4 Access
~~~~~~~~~~~~~~~~~~~

Tripoli4 is distributed by CEA/Saclay as a binary executable.  For
access to DAG-Tripoli4, please contact `Jean-Christophe Trama
<mailto:jean-christophe.trama@cea.fr>`_.

