Moab Tools
========================================

The tools/functions that prove useful here will be folded into Meshkit.

Among these tools is the make_watertight algorithm. This algorithm is used to 'seal' 
unwatertight models facted using DAGMC (http://svalinn.github.io/DAGMC). The accepted input
format is an .h5m file. After sealing a new file will be prouced with '_zip' appended 
to the original .h5m filename.



Dependencies
------------

This code relies on a fully installed version of either MOAB 4.6.0 or MOAB 4.6.2 which can be found at https://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB.


Building Code 
-------------

1) clone the code from this github repot using the mw-moab4.6.2 tag:

   ```git clone https://github.com/pshriwise/moab_tools```

   followed by:

   ```cd ./moab_tools```

   ```git checkout tags/mw-moab4.6.2```

2) go to `/moab_tools/make_watertight/` and update the Makefile so that the first line reads:

   ```include /path/to/MOAB/install/lib/moab.make```

   or create an environment variable MOAB_MAKE with this path

3) run `make` in this directory

The make wateright algorithm should now be ready to run!

Running Algorithm
-----------------

In order to run the algorithm simply input:

```make_watertight /path/to/file/filename.h5m```

There is also a checking utility, `check_watertight', that can be used to 
evaluate the success of the sealing algorithm. It is run using:

```check_watertight /path/to/file/filename_zip.h5m```

Testing
-------

Tests exist for make_watertight in `/moab_tools/make_watertight/test/`.
These tests can be installed by updating the Makefile in this particular directory
and then running `make`. The tests can then be run (from within this directory) using the command:

```test_cyl cyl.h5m```
