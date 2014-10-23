PyNE
-----

This documents describes a procedure to get, build, and install PyNE along with its 
most basic prerequisites, in order to be used by the DAGMC Toolkit.  

Requirements
+++++++++++++
PyNE has many features that require different packages to be installed.  For 
the DAGMC Toolkit, only those packages that are needed are installed with this
procedure.

These instructions assume you are in the dagmc_bld directory, i.e.
::
   prompt%> cd $HOME/dagmc_bld

PyNE
~~~~
Clone pyne from the pyne repository:
::
    prompt%> git clone git://github.com/makeclean/pyne
    prompt%> cd pyne
    prompt%> git checkout add_particles_tally
    prompt%> cd ..

Miniconda
~~~~~~~~~~~
Download the build and install script for Miniconda:
::
    prompt%> wget http://repo.continuum.io/miniconda/Miniconda-3.0.5-Linux-x86_64.sh


PyTAPS
~~~~~~
Download the PyTAPS tarball:
::
    prompt%> wget https://pypi.python.org/packages/source/P/PyTAPS/PyTAPS-1.4.tar.gz
    prompt%> tar zxvf PyTAPS-1.4.tar.gz


Build PyNE
++++++++++

Prepare and run the Miniconda script.  This creates a new subdirectory 'anaconda'
::
    prompt%> chmod +x Miniconda-3.0.5-Linux-x86_64.sh
    prompt%> ./Miniconda-3.0.5-Linux-x86_64.sh -b -p `pwd`/anaconda

Set up paths to know about anaconda:
:: 
    prompt%> export LD_LIBRARY_PATH=`pwd`/anaconda/lib:$LD_LIBRARY_PATH
    prompt%> export C_INCLUDE_PATH=`pwd`/anaconda/include:$C_INCLUDE_PATH
    prompt%> export CPLUS_INCLUDE_PATH=`pwd`/anaconda/include:$CPLUS_INCLUDE_PATH
    prompt%> export LIBRARY_PATH=`pwd`/anaconda/lib:$LIBRARY_PATH

    prompt%> export PATH=`pwd`/anaconda/bin:`pwd`/anaconda/usr/local/bin:$PATH

Install named prerequisites with conda:
::
    prompt%> conda install conda-build jinja2 nose setuptools pytables scipy cython cmake
    prompt%> conda install patchelf

Install PyTAPS:
::
    prompt%> cd PyTAPS-1.4
    prompt%> python setup.py build
    prompt%> python setup.py install --skip-build --prefix=`pwd`/../anaconda
    prompt%> cd ..

Install pyne:
::
    prompt%> cd pyne
    prompt%> python setup.py install -- --prefix=`pwd`/../anaconda --hdf5=$HOME/dagmc_bld/HDF5 
    prompt%> cd ..

Verify the PyNE library has just been built:
::
    prompt%> ls -l anaconda/lib/libpyne.so
