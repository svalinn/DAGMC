Make Watertight
========================================

This algorithm is used to 'seal' unwatertight models facted using DAGMC (http://svalinn.github.io/DAGMC). The accepted input
format is an .h5m file. After sealing a new file will be prouced with '_zip' appended 
to the original .h5m filename.

Running Algorithm
-----------------

In order to run the algorithm simply input:

```make_watertight /path/to/file/filename.h5m```

There is also a checking utility, `check_watertight', that can be used to 
evaluate the success of the sealing algorithm. It is run using:

```check_watertight /path/to/file/filename_zip.h5m```

Testing
-------

Tests for make_watertight will be installed in /your/dagmc/install/location/tests.

The tests can be run (from within this directory) using the command:

./make_watertight_cone_tests && ./make_watertight_cylinder_tests


