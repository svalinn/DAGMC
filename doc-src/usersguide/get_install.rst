Getting and Installing the DAGMC Toolkit
----------------------------------------

DAGMC is a toolkit that provides direct geometry support to any Monte
Carlo radiation transport code.  The primary code used for development
and testing of DAGMC is `MCNP5 <http://laws.lanl.gov/vhosts/mcnp.lanl.gov/mcnp5.shtml>`_,
developed by `Los Alamos National Laboratory <http://www.lanl.gov>`_
and distributed by the `Radiation Safety Information Computing Center
<http://rsicc.ornl.gov>`_.  There has also been experience with MCNPX
(LANL), Tripoli4 (CEA/Saclay),`FLUKA <http://www.fluka.org>`_ (CERN/INFN), 
and `Geant4 <http://www.geant4.cern.ch/>`_ (CERN).

These instructions describe the basic steps for downloading and
installing the software libraries that provide the DAGMC toolkit for
integration with Monte Carlo codes.  After this, code-specific
instructions will be give for each code. 

Toolkit Installation
++++++++++++++++++++++++++++

This section details the installation and build steps for the prerequisite packages for the the DAGMC
toolkit with specific physics codes.

Prerequisites
~~~~~~~~~~~~~

In order to install you must have done the following:

1) Cloned the `DAGMC <http://github.com/svalinn/DAGMC>`_ repository
2) Installed `HDF5 <http://www.hdfgroup.org/HDF5/>`_
3) Install `CGM <http://trac.mcs.anl.gov/projects/ITAPS/wiki/CGM>`_, using the --with-cubit option
4) Installed `MOAB <http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_,
   using options --with-cgm --with-hdf5 --with-dagmc --without-netcdf 
   If you need to prepare meshed geometries the following are also required
   a) Install `CUBIT <http://cubit.sandia.gov>`_ v12.2 or v13.1
5) Installed Lapack
   Note that the MCNP build automatically builds the dagtally library, which uses Lapack 
6) Installed `FLUKA <http://www.fluka.org>`_


Installation of Prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With installation of the DAGMC Toolkit, the dependency stack will look like this:

* Some physics package, e.g. MCNP5
   * `MOAB/DAGMC <http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_
   * `HDF5 <http://www.hdfgroup.org/HDF5/release/obtain5.html>`_
   * `CGM <http://trac.mcs.anl.gov/projects/ITAPS/wiki/CGM>`_ 
       * ACIS v19, or `CUBIT <http://cubit.sandia.gov>`_ v12.2 or v13.1 (with CGM `trunk <http://ftp.mcs.anl.gov/pub/fathom/cgm-nightly-trunk.tar.gz>`_ only)


Here are some assumptions/conventions that are used in these instructions:

* all operations are carried out in the a subdirectory ``dagmc_bld`` of a user's home directory
* the path to CUBIT files is known, e.g. ``/path/to/cubit``.  This is the directory that contains the Python script file ``cubit`` and a ``bin`` subdirectory.  
* all tarballs reside in user's home directory

If these do not apply to you, please modify your steps accordingly.

     *(For a shortcut to installing DAG-MCNP5.1.60 you may be able to use the build_dagmc_stack.bash script .)*

CGM
````

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


If installing CGM version 12.2 from a tarball, ``CGM-12.2.0.tar.gz``:
::
    prompt%> tar xzf ~/CGM-12.2.0.tar.gz
    prompt%> ln -s CGM-12.2.0 src

In all CGM cases:
::
    prompt%> cd bld
    prompt%> ../src/configure --enable-optimize \
              --enable-shared --disable-debug \
              --with-cubit=/path/to/cubit  \
              --prefix=$HOME/dagmc_bld/CGM
    prompt%> make
    prompt%> make install


HDF5
````

The HDF5 tarball can be downloaded from the `website <http://www.hdfgroup.org/HDF5/release/obtain5.html>`_ or, on a Linux machine, using the wget command, e.g.
::
    prompt%> wget http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.11/src/hdf5-1.8.11.tar.gz

See `ftp versions <https:// http://www.hdfgroup.org/ftp/HDF5/releases>`_ for available versions.

Create a directory and install HDF5:
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
````

Note:  MOAB version 4.7.0 is the earliest version that may be used.

Create a MOAB directory to install in
::
    prompt%> mkdir -p $HOME/dagmc_bld/MOAB/bld
    prompt%> cd $HOME/dagmc_bld/MOAB


If installing MOAB from the git repository:
::
    prompt%> git clone https://bitbucket.org/fathomteam/moab/
    prompt%> cd moab
    prompt%> git checkout 4.7.0
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
              --prefix=$HOME/dagmc_bld/MOAB
    prompt%> make
    prompt%> make install


Post Install
``````````````

Having installed all the prerequisite tools, namely Cubit, CGM, HDF5 and MOAB, the user
must ensure that the system has access to the libraries and programs that have been built.
Therefore modify the $PATH and $LD_LIBRARY_PATH environments accordingly:
:: 

    prompt%> export PATH=$PATH:/$HOME/dagmc_bld/path/to/cubit/bin: \
                               /$HOME/dagmc_bld/HDF5/bin: \
                               /$HOME/dagmc_bld/MOAB/bin
    prompt%> export PATH=$PATH:/$HOME/dagmc_bld/path/to/cubit/bin:  \
                               /$HOME/dagmc_bld/HDF5/lib: \
                               /$HOME/dagmc_bld/MOAB/lib:/$HOME/dagmc_bld/CGM/lib
 

Applying DAGMC to Specific Monte Carlo Codes
+++++++++++++++++++++++++++++++++++++++++++++

Install DAGMC
~~~~~~~~~~~~~

Clone the DAGMC repository
::
    prompt%> cd $HOME/dagmc_bld
    prompt%> git clone https://github.com/svalinn/DAGMC
    prompt%> cd DAGMC
    prompt%> git checkout develop

Note:  a version of the build instructions, INSTALL.rst, is in the DAGMC directory.

Install Instructions
~~~~~~~~~~~~~~~~~~~~

The DAGMC toolkit now has a full CMake install and build method, for all codes used downstream, even
going so far to replace the MCNP build method with a CMake file.

Note that in addition to the detailed instructions above for building the MOAB stack, you must also install

1) Lapack, using your favorite method
2) `Fluka <http://www.fluka.org>`_

How To
========
1) Copy the "Source" directory for MCNP5v16 from the LANL/RSICC CD to the mcnp5/ directory in the DAGMC source tree
   cp -r <path to cdrom/MCNP5/Source mcnp5/.
2) Apply the patch from the patch folder
   patch -p1 < patch/dagmc.patch.5.1.60
3) Assuming your patch was succesfully applied, i.e. there were no warnings or errors then we can now start building,
   assuming that you are in the base level of the DAGMC repo, create a new directory and navigate to it.
   mkdir bld
   cd bld
4) We can now configure DAGMC for building, you must include the CMAKE_INSTALL_PREFIX option with a folder where
   you would like the toolkit to be installed, this directory need not exist.

   For example we could just build the DAGMC interfaces and DAG-MCNP5
   cmake ../. -DMOAB_DIR=$MOAB_PATH/lib -DBUILD_MCNP5=ON -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
   or to building MCNP5 in parallel
   cmake ../. -DMOAB_DIR=$MOAB_PATH/lib -DBUILD_MCNP5=ON -DMPI_BUILD=ON \
   -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
   We could build Fluka as well
   cmake ../. -DMOAB_DIR=$MOAB_PATH/lib -DBUILD_MCNP5=ON -DMPI_BUILD=ON \
     -DBUILD_FLUKA=ON -DFLUKA_DIR=$FLUPRO -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
   or only Fluka
   cmake ../. -DMOAB_DIR=$MOAB_PATH/lib -DBUILD_FLUKA=ON -DFLUKA_DIR=$FLUPRO \
             -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
   and or Geant4
   cmake ../. -DMOAB_DIR=$MOAB_PATH/lib -DBUILD_MCNP5=ON -DMPI_BUILD=ON \
          -DBUILD_FLUKA=ON -DFLUKA_DIR=$FLUPRO -DBUILD_GEANT4=ON -DGEANT4_DIR=/mnt/data/opt/geant4  \
          -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

5) Assuming that the cmake step was succesful, i.e. no errors were reported then we can make by issuing the make command
   make
6) If there were no errors, then we can install the DAGMC suite of libraries and tools by issing the install command
   make install
7) If everything was successful, you may have the mcnp5 and mainfludag executables in the bin folder, the libraries in lib
   and the header files in the include folder


  


   
MOSTLY CUT AFTER HERE
DAG-MCNP5 Build/Install
~~~~~~~~~~~~~~~~~~~~~~~~

If you would like to use DAGMC with MCNP5, known as DAG-MCNP5, you will also need:

* MCNP5.1.60 source code from `RSICC <http://rsicc.ornl.gov>`_
* a local copy of UW-Madison's `DAGMC git repo <https://github.com/svalinn/DAGMC>`_ 

Automatic Installation
=======================

A package has been prepared that includes many of the requires
software libraries and an automated build script.  Because the DAGMC
team is not authorized to distribute `CUBIT
<http://cubit.sandia.gov>`_ nor `MCNP5.1.60 source code
<http://mcnp.lanl.gov>`_, you must acquire those through the
appropriate channels on your own.

Once you have both of those things, you should be able to use the
DagmcBuildPackage to create a working install of DAG-MCNP5.1.60.

Manual Installation
=====================

The following steps are required to install DAG-MCNP5.  Most of these steps are described in more detail below.

1. Install the DAGMC Toolkit as described above
2. Clone a copy of the DAGMC git repo.
3. Apply the appropriate patch from DAGMC/MCNP5/patch/ to your copy of the MCNP5 source code
4. Build & install the patched version of MCNP5

Some assumptions/conventions:

* all operations are carried out in the a subdirectory ``dagmc_bld`` of a user's home directory
* path to CUBIT files is known, e.g. ``/path/to/cubit``
* all tarballs reside in user's home directory
* MCNP5 source code is available in location ``$HOME/dagmc_bld/MCNP5``
* A cloned DAGMC git repo can be found at ``$HOME/dagmc_bld/DAGMC``

Apply Patch
############

*Apply DAGMC Patch to MCNP5 v1.60*

Perform the following steps:
::
    prompt%> cd $HOME/dagmc_bld/MCNP5
    prompt%> patch -p1 < /path/to/DAGMC/MCNP5/patch/dagmc.patch.5.1.60


Build DAG-MCNP5
################

*Build DAG-MCNP5 from modified code*

One of the easiest ways to build DAG-MCNP5 is directly using the
``makefile`` from the command-line.  To do this, you must know the
``makefile`` options to build a version of MCNP5 without DAGMC,
usually in the form:
::
    prompt%> make build CONFIG="seq plot gfortran" FC=gfortran MARCH=M64``

or similar.  Starting from these options, you can build DAG-MCNP5 from
a patched source code with:
::
    prompt%> make build CONFIG="seq plot gfortran dagmc" FC=gfortran MARCH=M64 \
                 MOAB_DIR=$HOME/dagmc_bld/MOAB CUBIT_DIR=/path/to/cubit/bin \
		 DAGMC_DIR=$HOME/dagmc_bld/DAGMC/MCNP5/dagmc


If you are less familiar with building MCNP5 from the ``makefile`` you
may want to use the interactive ``install`` script provided by LANL:
::
    prompt%> ./install

Within the ``install`` program you will need to set the DAGMC build options:

* turn on DAGMC mode
* provide the path to MOAB: ``$HOME/dagmc_bld/MOAB``
* provide the path to CUBIT: ``/path/to/cubit``

Your executable should be available as ``$HOME/dagmc_bld/MCNP5/Source/src/mcnp5``.

DAG-Tripoli4 Access
~~~~~~~~~~~~~~~~~~~

Tripoli4 is distributed by CEA/Saclay as a binary executable.  For
access to DAG-Tripoli4, please contact `Jean-Christophe Trama
<mailto:jean-christophe.trama@cea.fr>`_.

FluDAG Build
~~~~~~~~~~~~

FluDAG uses `FLUKA <http://www.fluka.org>`_ from CERN/INFN with the DAGMC Toolkit.
The steps to build and install FluDAG follow.

*Install  FLUKA*
==================

In order to download FLUKA you need to become a registered user, which you can do at 
the `FLUKA register <https://www.fluka.org/fluka.php?id=secured_intro>`_ page from a link on the main FLUKA page.
Save the user id and password for future FLUKA updates. We recommend a x64 worfklow and as such you should download
the 64 bit exectubale the download name is of the form *fluka20xx.xx-linux-gfor64bitAA.tar.gz*, `See <http://www.fluka.org/fluka.php?id=ins_run&mm2=3>`_.

Once the FLUPRO environment variables have been set and you have a confirmed working install of Fluka, please proceed to
the FluDAG install below.

*FluDAG Installation*
=====================
 
Get the FluDAG Development release repository by cloning :ref:`DAGMC` (if you haven't done so)
::

    prompt%> cd $HOME/dagmc_bld/DAGMC
    prompt%> git checkout develop

In order to install we run CMake and provide the path to the MOAB installation, the $FLUPRO
path is picked up implcitly
::

    prompt%> cd $HOME/DAGMC/FluDAG
    prompt%> mkdir bld
    prompt%> cd bld
    prompt%> cmake ../src/. -DMOAB_HOME=$HOME/dagmc_bld/MOAB
    prompt%> make
    prompt%> mv src/mainfludag .

Upon successful compilation the directory *bld* will have the *mainfludag* executable in it. 

Fluka is typically run by users with the rfluka script, we patch this script below, it will still
be compatable with standard Fluka inputs
::

    prompt%> cp $FLUPRO/flutil/rfluka $FLUPRO/flutil/rfluka.orig 
    prompt%> patch -p1 < ../src/rfluka.patch $FLUPRO/flutil/rfluka

In order to test FluDAG, an environment variable, named 'FLUDAG', with the path to the *bld* 
directory must be set:

Add this statement to your login script:
::

    export FLUDAG=${HOME}/dagmc_bld/DAGMC/FluDAG/bld/

The FluDAG tests are in a separate repository, which can be cloned from github
::

    prompt%> cd $HOME/dagmc_bld
    prompt%> git clone https://github.com/svalinn/fludag_testing.git

To run all the tests type:
::

    prompt%> cd $HOME/fludag_testing
    prompt%> ./run_test test_input

Some of the tests are slow, so the above command will take some time.
If you want to run just the fast tests, or just the magnetic tests:
:: 

    prompt%> ./run_test test_fast
    prompt%> ./run_test test_magnetic

The slow tests can be run separately:
::

    prompt%> ./run_test test_slow

Some of the tests check the installation and can be run separately:
::

    prompt%> ./run_test test_install


