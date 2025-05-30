Installing dependencies
=======================

This document describes how to install the dependencies of DAGMC. You can
either follow a package manager based route or perform a source install.
Depending upon your preferences and circumstances one may suit you rather
than the other, however, soon it will be possible to install all dependencies
using a package manager.

The following dependencies are required:

    * Eigen3_
    * HDF5_
    * MOAB_

For users following the source install route, we assume that you are building
DAGMC in the subdirectory ``dagmc_bld`` of your home directory; i.e.
::

    $ cd $HOME
    $ mkdir dagmc_bld
    $ cd dagmc_bld

Eigen3
~~~~~~

The best way to install the Eigen3 package is with your package manager.

Package manager installation
----------------------------

Debian linux users can install the latest HDF5 release with:
::

    $ sudo apt-get install libeigen3-dev

Redhat linux users can do likewise with:
::

    $ sudo yum install eigen3-devel

HDF5
~~~~

Source installation
-------------------

The tarball containing the HDF5 source code can also be downloaded from the
`HDF5 website <HDF5_>`_. Note that if you choose this option, we recommend you
obtain HDF5 version 1.14.3 instead of the newest version. The following commands
can be used to install HDF5 from source.
::

    $ mkdir -p $HOME/dagmc_bld/HDF5/build
    $ cd $HOME/dagmc_bld/HDF5
    $ git clone -b 1.14.3 --depth 1 https://github.com/HDFGroup/hdf5.git && \
    $ cd build
    $ cmake ../hdf5 \
        -DCMAKE_INSTALL_PREFIX=$HOME/dagmc_bld/HDF5 \
        -DBUILD_SHARED_LIBS=ON
    $ make
    $ make check
    $ make install

Package manager installation
----------------------------

Debian linux users can install the latest HDF5 release with:
::

    $ sudo apt-get install libhdf5-dev

Redhat linux users can do likewise with:
::

    $ sudo yum install -y epel-release
    $ sudo yum config-manager --enable epel
    $ sudo yum install hdf5-devel

MOAB installation
~~~~~~~~~~~~~~~~~

As of DAGMC version 3.2.3, MOAB version 5.5.1 is preferred. The following
commands can be used to download MOAB from its `source repository <MOAB_>`_ and
set it up for building.
::

    $ cd $HOME/dagmc_bld
    $ mkdir -p MOAB/bld
    $ cd MOAB
    $ git clone --depth 1 -b 5.5.1 https://bitbucket.org/fathomteam/moab

To build moab using the default packages:
::

    $ cd bld
    $ cmake ../moab -DENABLE_HDF5=ON -DHDF5_ROOT=$HOME/dagmc_bld/HDF5 \
              -DCMAKE_BUILD_TYPE=Release \
              -DENABLE_PYMOAB=ON \
              -DENABLE_BLASLAPACK=OFF \
              -DENABLE_FORTRAN=OFF \
              -DCMAKE_INSTALL_PREFIX=$HOME/dagmc_bld/MOAB \
              -DBUILD_SHARED_LIBS=ON \
              -DCMAKE_INSTALL_RPATH=$HOME/dagmc_bld/HDF5/lib:$HOME/dagmc_bld/MOAB/lib
    $ make
    $ make check
    $ make install

If you whish to avoid building MOAB from source, it is possible to build it at
the same time as DAGMC by setting the ``PULL_INSTALL_MOAB`` option to the 
desired MOAB version. This will automatically download and build MOAB as a 
dependency of DAGMC. When using this option it is assumed that HDF5 is properly
installed on your computer and you will need to provide its location for the 
MOAB installation. For example, to build MOAB version 5.5.1, you would add the 
following flag to the cmake command when building DAGMC:
::

    $ cmake ../dagmc -DPULL_INSTALL_MOAB=5.5.1 -DHDF5_ROOT=PATH_TO_HDF5

Making sure the dependencies were installed correctly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to test HDF5 and MOAB, some of their directories must be added to your
``$PATH`` and ``$LD_LIBRARY_PATH``.
::

    $ export PATH=$PATH:$HOME/dagmc_bld/HDF5/bin
    $ export PATH=$PATH:$HOME/dagmc_bld/MOAB/bin
    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/dagmc_bld/HDF5/lib
    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/dagmc_bld/MOAB/lib

The following commands can be used to test whether HDF5 and MOAB were built
successfully.
::

    $ which h5ls
    $ which mbconvert

If you have installed the dependencies correctly, you are now ready to
`install DAGMC <dagmc.html>`_.

..  _HDF5: http://www.hdfgroup.org/HDF5
..  _MOAB: https://sigma.mcs.anl.gov/moab-library/
..  _Eigen3: http://eigen.tuxfamily.org/index.php