PyNE Installation
-----------------

This documents describes a procedure to get, build, and install PyNE, along 
with its most basic prerequisites, in order to be used by the DAGMC Toolkit.  

Requirements
+++++++++++++
PyNE has many features that require different packages to be installed.  For 
the DAGMC Toolkit, only those packages that are needed are installed with this
procedure.

These instructions assume 
 - you have built an HDF5-dependent MOAB, per the instructions `here <get_install.html>`_.
 - You have set your LD_LIBRARY_PATH environment variable per the instructions in 
   the Post Install `section <get_install.html#post-install>`_.
   

If you are installing on a Debian-type machine (e.g. Ubuntu), you may perform
a system install of many of the PyNE prerequisites with the following commands: 
::
    prompt%> sudo apt-get install python-numpy python-scipy
    prompt%> sudo apt-get install python-tables python-nose cython

If you are on a system that does not permit you to do system installs, you will
need to install PyNE dependencies per the instructions at the 
`PyNE site <http://pyne.io/install.html>`_.  
   
PyTAPS
++++++++++
Add the LD_LIBRARY_PATH to $LIBRARY_PATH and set up $CPLUS_INCLUDE_PATH and $C_INCLUDE_PATH
::
    prompt%> export LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBRARY_PATH
    prompt%> export CPLUS_INCLUDE_PATH=$HOME/dagmc_bld/MOAB/include:$CPLUS_INCLUDE_PATH
    prompt%> export C_INCLUDE_PATH=$HOME/dagmc_bld/MOAB/include:$C_INCLUDE_PATH

Download the PyTAPS tarball:
::
    prompt%> cd $HOME/dagmc_bld
    prompt%> wget https://pypi.python.org/packages/source/P/PyTAPS/PyTAPS-1.4.tar.gz
    prompt%> tar zxvf PyTAPS-1.4.tar.gz
    prompt%> cd PyTAPS
    prompt%> python setup.py --iMesh-path=$HOME/dagmc_bld/MOAB install --user
    
(Note: you may get error messages)

PyNE
++++

The following instructions assume you are in the dagmc_bld directory, i.e.
::
    prompt%> cd $HOME/dagmc_bld

Clone pyne from the pyne repository:
::
    prompt%> git clone git://github.com/pyne/pyne.git

Move to the new pyne repository and ensure you are on the 'develop' branch:
::
    prompt%> cd pyne
    prompt%> git checkout develop

Run the python command to build pyne:
::
    prompt%> python setup.py install --hdf5=path/to/hdf5 --user

Create the database of nuclear materials:
::
    prompt%> cd scripts
    prompt%> ./nuc_data_make

