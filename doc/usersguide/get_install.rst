Installing the DAGMC Toolkit
----------------------------------------
DAGMC is a toolkit that provides direct geometry support to any Monte
Carlo radiation transport code. These instructions describe the basic 
steps for downloading and installing the software libraries that provide 
the DAGMC toolkit for integration with Monte Carlo codes.  After this, 
code-specific instructions will be give for each code. 

Toolkit Installation
++++++++++++++++++++++++++++
This section details the installation and build steps for the prerequisite 
packages for the the DAGMC toolkit with specific physics codes.

Requirements
~~~~~~~~~~~~~

In order to install you must have done the following:

 * Cloned the `DAGMC <http://github.com/svalinn/DAGMC>`_ repository
 * Installed `HDF5 <http://www.hdfgroup.org/HDF5/>`_
 * Installed `MOAB <http://sigma.mcs.anl.gov/moab-library/>`_,
   using options --with-hdf5 --enable-dagmc 
 * Installed Lapack.  __Note:__ the MCNP build automatically builds the dagtally library, which requires Lapack 
 * Installed `FLUKA <http://www.fluka.org/fluka.php>`_ - and/or - 
 * Installed `Geant4 <http://geant4.cern.ch/>`_


Building
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
With installation of the DAGMC Toolkit, the dependency stack will look like this:

* Some physics package, e.g. MCNP5, Geant4 etc
   * `MOAB/DAGMC <http://bitbucket.org/fathomteam/moab>`_
   * `HDF5 <http://www.hdfgroup.org/HDF5/release/obtain5.html>`_

Assumptions and conventions that are used in these instructions:

* all operations are carried out in the a subdirectory ``dagmc_bld`` of a user's home directory

If these do not apply to you, please modify your steps accordingly.

     *(For a shortcut to installing DAG-MCNP5.1.60 you may be able to use the build_dagmc_stack.bash script .)*

HDF5
======
Debian linux users may conveniently install the latest HDF5 release with the command:
::
    prompt%> sudo apt-get install hdf5-dev

Redhat linux users can do likewise with this command:
::
    prompt%> sudo yum install hdf5-dev

Otherwise, the HDF5 tarball can also be downloaded from the HDF5 `website <http://www.hdfgroup.org/HDF5/release/obtain5.html>`_.  
::
    prompt%> wget \
    http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.13/src/hdf5-1.8.13.tar.gz

See the `HDF5 ftp site <http://www.hdfgroup.org/ftp/HDF5/releases>`_ for available versions.

In the case of a tarball, create a directory and install HDF5:
::
    prompt%> mkdir -p $HOME/dagmc_bld/HDF5/bld
    prompt%> cd $HOME/dagmc_bld/HDF5
    prompt%> tar xzf ~/hdf5-1.8.13.tar.gz
    prompt%> ln -s hdf5-1.8.13 src
    prompt%> cd bld
    prompt%> ../src/configure --enable-shared --prefix=$HOME/dagmc_bld/HDF5
    prompt%> make
    prompt%> make install


MOAB
======

The master branch of MOAB is currently at version 4.9.2, which is the earliest version that may be used. We have
recently updated the DAGMC interface which is now on Version 2.0

Create a MOAB directory to install in
::
    prompt%> mkdir -p $HOME/dagmc_bld/MOAB/bld
    prompt%> cd $HOME/dagmc_bld/MOAB

If installing MOAB from the git repository:
::
    prompt%> git clone https://bitbucket.org/fathomteam/moab/
    prompt%> cd moab
    prompt%> git checkout master
    prompt%> autoreconf -fi
    prompt%> cd ..
    prompt%> ln -s moab src

The command to "git checkout master" is, in general, redundant but is included here for completeness.

In all MOAB cases:
::
    prompt%> cd bld
    prompt%> ../src/configure --enable-optimize \
              --enable-shared --disable-debug \
              --with-hdf5=$HOME/dagmc_bld/HDF5 \
              --enable-dagmc \
              --prefix=$HOME/dagmc_bld/MOAB
    prompt%> make
    prompt%> make install


Post Install
~~~~~~~~~~~~~~

Having installed all the prerequisite tools, HDF5, and MOAB, the user
must ensure that the system has access to the libraries and programs that have been built.
Therefore modify the $PATH and $LD_LIBRARY_PATH environments accordingly:
:: 

    prompt%> export PATH=$PATH:$HOME/.local/bin: \
                               $HOME/dagmc_bld/HDF5/bin: \
                               $HOME/dagmc_bld/MOAB/bin
    prompt%> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH: \
                               $HOME/.local/lib: \
                               $HOME/dagmc_bld/HDF5/lib: \
                               $HOME/dagmc_bld/MOAB/lib

One should be able to sucessfully run the commands
::
   prompt%> which mbconvert
   prompt%> which h5ls

This is indicative of a succesful depdendency build.

Toolkit Applications
+++++++++++++++++++++++++++++++++++++++++++++

Install FLUKA
~~~~~~~~~~~~~~
FluDAG uses `FLUKA <http://www.fluka.org/fluka.php>`_ from CERN/INFN with the DAGMC Toolkit.

In order to download FLUKA you need to become a registered user, which you can do at 
the `FLUKA register <https://www.fluka.org/fluka.php?id=secured_intro>`_ page 
from a link on the main FLUKA page.

Save the user id and password for future FLUKA updates.  We recommend an x64 
worfklow and as such you should download the 64 bit executable.  The download 
name is of the form *fluka20xx.xx-linux-gfor64bitAA.tar.gz*.  See the 
`site <http://www.fluka.org/fluka.php?id=ins_run&mm2=3>`_ for instructions.

Follow the FLUKA site instructions to set the FLUPRO and FLUFOR environment 
variables.  Currently, you must patch FLUKA's run script, rfluka, to allow for some DAGMC
specific options.
::
    prompt%> cd $FLUPRO/flutil
    prompt%> cp rfluka rfluka.orig
    prompt%> patch rfluka $HOME/dagmc_bld/DAGMC/fluka/rfluka.patch

Confirm that you have a working install of Fluka and proceed to the next steps.

Install Geant4
~~~~~~~~~~~~~~~~
`Geant4 <http://geant4.cern.ch>`_, a toolkit for the simulation of the passage of particles through matter, can be found 
`here <http://geant4.cern.ch/support/gettingstarted.shtml>`_,  including a link to instructions for installation. We recommend the following flags
::
   -DCMAKE_INSTALL_PREFIX=<path to install location>
   -DGEANT4_INSTALL_DATA=ON
   -DGEANT4_USE_QT=ON or -DGEANT4_USE_OPENGL_X11=ON
   -DGEANT4_USE_SYSTEM_EXPAT=OFF


Build DAGMC Interfaces
~~~~~~~~~~~~~~~~~~~~~~
The DAGMC toolkit now has a full CMake install and build method for all codes used downstream.  It even
replaces the MCNP build method with a CMake file. Note that in addition to the detailed instructions above 
for building the MOAB stack, you may also need to install Lapack using, for example, "sudo apt-get install 
liblapack-dev libblas-dev".

DAGMC Build Procedure
~~~~~~~~~~~~~~~~~~~~~
Clone the DAGMC repository
::
    prompt%> cd $HOME/dagmc_bld
    prompt%> git clone https://github.com/svalinn/DAGMC
    prompt%> cd DAGMC
    prompt%> git checkout develop

If building MCNP5 one must populate and patch the MCNP5 source in the DAGMC subdirectory 
first. Copy the "Source" directory for MCNP5v16 from the LANL/RSICC CD to the 
mcnp/mcnp5 directory in the DAGMC source tree
::
    prompt%> cd $HOME/dagmc_bld/DAGMC/mcnp/mcnp5
    prompt%> cp -r <path to cdrom>/MCNP5/Source .

Apply the patch from the mcnp5 folder of the DAGMC source tree, i.e. dagmc/mcnp/mcnp5
::
    prompt%> patch -p0 < patch/dagmc.patch.5.1.60

Assuming the patch was succesfully applied, i.e. there were no warnings or 
errors, we can now configure the DAGMC cmake system for the desired build.  
Configuration
~~~~~~~~~~~~~~~~
The CMake system can be used to configure a build of any or all of the 
following, see `cmake options <cmake_options.html>`_ for a list of all possible options, 
which include

   * MCNP5 with or without MPI
   * GEANT4 (DagSolid)
   * FLUKA  (FluDAG)
   * TALLY (Tally interface)
   
You will need to include the CMAKE_INSTALL_PREFIX=install_dir option as part of the configuration.  When the 
build command 'make install' is invoked, libraries, executables, tests, and include files are installed in 
subdirectories under install_dir.  It is common to use -DCMAKE_INSTALL_PREFIX=..', which creates and populates 
these directories one level above the build directory, that is, in the DAGMC directory.  
Note that the '-D' in front of CMAKE_INSTALL_PREFIX, and all of the configuration variables, defines the variable
for the cmake system.

From the base level of the DAGMC repository create a build directory and navigate to it.
::
    prompt%> cd $HOME/dagmc_bld/DAGMC
    prompt%> mkdir bld
    prompt%> cd bld

In the examples, the environment variable, "INSTALL_PATH", can point to any location
where you want the libraries ($INSTALL_PATH/lib), executables ($INSTALL_PATH/bin), and
other build products to be installed.  It is typically set to the DAGMC directory, i.e.
::
    prompt%> export INSTALL_PATH=$HOME/dagmc_bld/DAGMC

**Example 1:**  Build the DAGMC interfaces and DAG-MCNP5, assuming that 
the DATAPATH environment variable is undefined.
::
    prompt%> cmake ../. -DBUILD_MCNP5=ON -DMCNP5_DATAPATH=<path to MCNP data> \
                        -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
**Example 2:**  Build MCNP5 in parallel.
::
    prompt%> cmake ../. -DBUILD_MCNP5=ON -DMPI_BUILD=ON \
                        -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 3:**  Build MCNP5 in parallel and build the dagmc-enabled FLUKA.
Note that $FLUPRO should have been previously defined as part of the FLUKA install.
::
    prompt%> cmake ../. -DBUILD_MCNP5=ON -DMPI_BUILD=ON \
                        -DBUILD_FLUKA=ON -DFLUKA_DIR=$FLUPRO \
			-DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 4:** Build only FluDAG.
::
    prompt%> cmake ../. -DBUILD_FLUKA=ON -DFLUKA_DIR=$FLUPRO \
                        -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 5:**  Build MCNP, FluDAG and Geant4-enabled DAGMC.
::
    prompt%> cmake ../. -DBUILD_MCNP5=ON  -DMPI_BUILD=ON \
                        -DBUILD_FLUKA=ON  -DFLUKA_DIR=$FLUPRO \
			-DBUILD_GEANT4=ON -DGEANT4_DIR=path/to/geant4 \
                        -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

**Example 6:**  Build MCNP, FluDAG, Geant4-enabled DAGMC and the Tally library and tests.
::
    prompt%> cmake ../. -DBUILD_MCNP5=ON  -DMPI_BUILD=ON \
                        -DBUILD_FLUKA=ON  -DFLUKA_DIR=$FLUPRO \
			-DBUILD_GEANT4=ON -DGEANT4_DIR=/path/to/geant4 \
			-DBUILD_TALLY=ON \
                        -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

Compile and Install
~~~~~~~~~~~~~~~~~~~~~

Assuming that the CMake step was succesful, i.e. no errors were reported, compile by issuing the make command:
::
    prompt%> make

If there were no errors, install the DAGMC suite of libraries and tools by issuing the install command:
::
    prompt%> make install

If everything was successful, you may have the mcnp5 and mainfludag executables in the $INSTALL_PATH/bin folder, 
the libraries in $INSTALL_PATH/lib and the header files in the $INSTALL_PATH/include folder.

Post Install
~~~~~~~~~~~~
If your build was successful, you must add the dagmc_bld/lib folder to your LD_LIBRARY_PATH
::
    prompt%> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INSTALL_PATH/lib

Testing
~~~~~~~~~~~

We regularly run the DAGMC test suite as part of our continuous integration system, for which we use 
`Travis <https://travis-ci.org/svalinn/DAGMC>`_. You may however, wish to run the tests in the 
$INSTALL_PATH/tests directory to verify correct installation.  To do this requires
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

Note that the path to geant4make.sh is different from the path to the geant4 install 
directory, defined with -DGEANT4_DIR=path/to/geant4, in the DAGMC compilation examples.

Again, with successful execution the last few lines of screen output are:
::
    [       OK ] DagSolidTest.surface_area_test (5 ms)
    [----------] 16 tests from DagSolidTest (228 ms total)

    [----------] Global test environment tear-down
    [==========] 16 tests from 1 test case ran. (228 ms total)
    [  PASSED  ] 16 tests.

With testing successfully completed you are now ready to run your first DAGMC `problem <workflow/uw2.html>`_.

DAG-Tripoli4 Access
~~~~~~~~~~~~~~~~~~~

Tripoli4 is distributed by CEA/Saclay as a binary executable.  For
access to DAG-Tripoli4, please contact `Jean-Christophe Trama
<mailto:jean-christophe.trama@cea.fr>`_.

