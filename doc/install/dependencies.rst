Installing dependencies
=======================

This document describes how to install the dependencies of DAGMC. You can
either follow a package manager based route or perform a source install.
Depending upon your preferences and circumstances one may suit you rather
than the other, however, soon it will be possible to install all dependencies
using a package manager.

The following dependencies are required:

    * LAPACK_
    * HDF5_
    * MOAB_

For users following the source install route, we assume that you are building
DAGMC in the subdirectory ``dagmc_bld`` of your home directory; i.e.
::

    $ cd $HOME
    $ mkdir dagmc_bld
    $ cd dagmc_bld

LAPACK
~~~~~~

Source installation
-------------------

If you don't have administrator privileges, refer to the
`LAPACK website <LAPACK_>`_ for information on how to build
LAPACK from source.

Package manager installation
----------------------------

Ubuntu/Debian linux users can install LAPACK with:
::

    $ sudo apt-get install libblas-dev liblapack-dev

Redhat linux users can do likewise with:
::

    $ sudo yum install libblas-dev liblapack-dev

HDF5
~~~~~~

Source installation
-------------------

If electing to install HDF5 from source. The tarball containing the HDF5
source code can also be downloaded from the `HDF5 website <HDF5_>`_.
Note that if you choose this option, we recommend you obtain HDF5 version 1.8.13
instead of the newest version. The following commands can be used to install
HDF5 from source.
::

    $ mkdir -p $HOME/dagmc_bld/HDF5/bld
    $ cd $HOME/dagmc_bld/HDF5
    $ wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.13/src/hdf5-1.8.13.tar.gz
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

    $ sudo apt-get install hdf5-dev

Redhat linux users can do likewise with:
::

    $ sudo yum install hdf5-dev

MOAB installation
~~~~~~~~~~~~~~~~~

As of DAGMC version 2.0, MOAB version 4.9.2 or higher is required. The following
commands can be used to download MOAB from its `source repository
<MOAB_>`_ and set it up for building.
::

    $ cd $HOME/dagmc_bld
    $ mkdir -p MOAB/bld
    $ cd MOAB
    $ git clone https://bitbucket.org/fathomteam/moab
    $ cd moab
    $ git checkout Version4.9.2
    $ autoreconf -fi
    $ cd ..
    $ ln -s moab src

If you have followed the source install route, then the following commands
should be used to build MOAB. Note that the ``--enable-dagmc`` configure option
is required.
::

    $ cd bld
    $ ../src/configure --enable-dagmc \
                       --enable-optimize \
                       --enable-shared \
                       --disable-debug \
                       --with-hdf5=$HOME/dagmc_bld/HDF5 \
                       --prefix=$HOME/dagmc_bld/MOAB
    $ make
    $ make check
    $ make install

If you have followed the package manager install route, then the following
commands should be used to build MOAB. Note that the ``--enable-dagmc``
configure option is required.
::

    $ cd bld
    $ ../src/configure --enable-dagmc \
                       --enable-optimize \
                       --enable-shared \
                       --disable-debug \
                       --with-hdf5 \
                       --prefix=$HOME/dagmc_bld/MOAB
    $ make
    $ make check
    $ make install


Making sure the dependencies were installed correctly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you installed HDF5 from source, you will need to make sure the system can
find it when it comes time to build DAGMC. This is done by adding some
directories to your ``$PATH`` and ``$LD_LIBRARY_PATH``. (This is not required if
you used a package manager to install HDF5.)
::

    $ export PATH=$PATH:$HOME/dagmc_bld/HDF5/bin
    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/dagmc_bld/HDF5/lib

You will also need to make sure the system can find MOAB.
::

    $ export PATH=$PATH:$HOME/dagmc_bld/MOAB/bin
    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/dagmc_bld/MOAB/lib

After including HDF5 and MOAB in your paths as described above, the following
commands can be used to test whether HDF5 and MOAB were built successfully.
::

    $ which h5ls
    $ which mbconvert

Note that your ``$PATH`` and ``$LD_LIBRARY_PATH`` will revert to their original
state when you open a new terminal, so it may be a good idea to add these
``export`` commands to your ``.bashrc`` file.

If you have installed the dependencies corretly, you are now ready to
`install DAGMC <dagmc.html>`_.

..  _LAPACK: http://www.netlib.org/lapack
..  _HDF5: http://www.hdfgroup.org/HDF5
..  _MOAB: http://press3.mcs.anl.gov/sigma/moab-library
