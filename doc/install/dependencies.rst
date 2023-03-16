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
obtain HDF5 version 1.8.13 instead of the newest version. The following commands
can be used to install HDF5 from source.
::

    $ mkdir -p $HOME/dagmc_bld/HDF5/bld
    $ cd $HOME/dagmc_bld/HDF5
    $ wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.13/src/hdf5-1.8.13.tar.gz
    $ tar -xzvf hdf5-1.8.13.tar.gz
    $ ln -s hdf5-1.8.13 src
    $ cd bld
    $ ../src/configure --enable-shared \
                       --prefix=$HOME/dagmc_bld/HDF5
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

    $ sudo yum install libhdf5-dev

MOAB installation
~~~~~~~~~~~~~~~~~

As of DAGMC version 3.2, MOAB version 5.3.0 or higher is required. The following
commands can be used to download MOAB from its `source repository <MOAB_>`_ and
set it up for building.
::

    $ cd $HOME/dagmc_bld
    $ mkdir -p MOAB/bld
    $ cd MOAB
    $ git clone --depth 1 -b 5.3.0 https://bitbucket.org/fathomteam/moab

To build moab using the default packages:
::

    $ cd bld
    $ cmake ../moab -DENABLE_HDF5=ON -DHDF5_ROOT=${hdf5_install_dir} \
              -DENABLE_PYMOAB=ON \
              -DENABLE_BLASLAPACK=OFF \
              -DENABLE_FORTRAN=OFF \
              -DCMAKE_INSTALL_PREFIX=$HOME/dagmc_bld/MOAB \
              -DCMAKE_C_COMPILER=${CC} \
              -DCMAKE_CXX_COMPILER=${CXX} \
              -DBUILD_SHARED_LIBS=ON \
              -DCMAKE_INSTALL_RPATH=${hdf5_install_dir}/lib:${moab_install_dir}/lib
    $ make
    $ make check
    $ make install

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

If you have installed the dependencies corretly, you are now ready to
`install DAGMC <dagmc.html>`_.

..  _HDF5: http://www.hdfgroup.org/HDF5
..  _MOAB: http://press3.mcs.anl.gov/sigma/moab-library
..  _Eigen3: http://eigen.tuxfamily.org/index.php