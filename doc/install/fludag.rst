..  |DAG-Code| replace:: FluDAG

..  include:: header.txt

..  include:: get_dagmc.txt

..  _install_fluka:

Install FLUKA
~~~~~~~~~~~~~

FluDAG uses `FLUKA <http://www.fluka.org/fluka.php>`_ from CERN/INFN. In order
to download FLUKA you need to become a registered user, which you can do at the
`FLUKA register <https://www.fluka.org/fluka.php?id=secured_intro>`_ page.

Save your user ID and password for future FLUKA updates. We recommend an x64
worfklow and thus you should download the 64-bit executable. The name of the
downloaded tarball is of the form ``fluka20xx.xx-linux-gfor64bitAA.tar.gz``.
Refer to the
`installation instructions <http://www.fluka.org/fluka.php?id=ins_run&mm2=3>`_
when building FLUKA.

Take care to follow the FLUKA site instructions when setting the
``$FLUPRO`` and ``$FLUFOR`` environment variables.

..  include:: configure_dag-code_header.txt

The following CMake command will build FluDAG. Note that ``$FLUPRO`` should have
previously been defined as part of the FLUKA install.
::

    $ cmake .. -DBUILD_FLUKA=ON \
               -DFLUKA_DIR=$FLUPRO \
               -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

If the CMake configuration proceeded successfully, you are now ready to install
DAGMC.

..  include:: build_dagmc.txt

Test |DAG-Code|
~~~~~~~~~~~~~~~

To run the FluDAG unit tests, use
::

    $ cd $INSTALL_PATH/tests
    $ ./fludag_unit_tests

If the tests ran successfully, the last few lines of the screen output will look
like this:
::

    [       OK ] FluDAGTest.GFireGoodPropStep (5 ms)
    [----------] 3 tests from FluDAGTest (108 ms total)

    [----------] Global test environment tear-down
    [==========] 3 tests from 1 test case ran. (108 ms total)
    [  PASSED  ] 3 tests.

..  include:: test_dagmc.txt
