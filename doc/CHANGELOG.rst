================
DAGMC Changelog
================

.. current developments

v3.2.0
====================

**Added:**

* `DagMC::git_sha` method (PR#707)
* Added `libeigen3-dev` package to be installed by `apt` (PR#683)
* Tool for checking DagMC models for overlaps. (PR#641)
* turn off fortran support in MOAB build (PR#684)
* use '&&' to join successive build steps to fail build on first failure (PR#684)
* remove specification of unused FORTRAN compiler for HDF5 & MOAB builds (PR#684)
* Note on adding DagMC libraries to LD_LIBRARY_PATH for OpenMC usage. (PR#625)
* Link to installation instructions for OpenMC with DagMC. (PR#625)
* Documentation section on material assignment by name for OpenMC model
  prep. (PR#616)
* MCNP6 version of pyne mesh source.F90. (PR#604)
* Documentation on model prep for OpenMC simulations with a DAGMC
  geometry. (PR#599)
* Added a script which updates amalgamated pyne. (PR#595)
* OpenMC material writing capability for UWUW workflow via updates to the
  amalgamated PyNE source. (PR#594)
* Community best practices including issue templates, pull request templates,
  contributing guidelines, and code of conduct. (PR#589)
* Added PyNE mesh source routine functionality. This can be enabled by setting
  the ``BUILD_MCNP_PYNE_SOURCE`` compiler definition to ``ON``. (PR#585)
* Added ability for users to disable building static or shared libraries (PR#572)
* Patch file for DAG-MCNP6.2 (PR#569)
* Default to a Release build. This results in optimization flags being used
  everywhere as appropriate. (PR#555)

  * Note that MCNP is still configured to use no more than ``-O1``
* Add macros to ``cmake/DAGMC_macros.cmake``. This results in much less
  duplicated cmake elsewhere. The following macros were added: (PR#555)

  * ``dagmc_setup_build``: Sets core variables used throughout the rest of the
    project.
  * ``dagmc_setup_options``: Defines cmake build options.
  * ``dagmc_setup_flags``: Defines compiler flags.
  * ``dagmc_get_link_libs``: Used by the ``dagmc_install_X`` macros to
    determine the names of the libraries that need to be linked.
  * ``dagmc_make_configure_file``: Setup the ``DAGMCConfig.cmake`` file.
  * ``dagmc_install_library``: Install a library.
  * ``dagmc_install_exe``: Install an executable.
  * ``dagmc_install_test``: Install a unit test.
  * ``dagmc_install_test_file``: Install a file needed for unit testing.

* Add a ``FindFluka.cmake`` file to find the Fluka library. (PR#555)
* Add ``RPATH`` functionality so that all executables and libraries
  automatically know where their dependencies are located. This removes the need
  for users to add anything to their ``LD_LIBRARY_PATH``. This can be turned off
  by setting ``-DBUILD_RPATH=OFF``. (PR#555)
* Add ability to build with position-independent code (PIC). This can be turned
  on by setting ``-DBUILD_PIC=ON``. (PR#555)
* Add options to enable/disable building all optional functionality. The
  following options were added: (PR#555)
  
  * ``BUILD_BUILD_OBB``
  * ``BUILD_MAKE_WATERTIGHT``
  * ``BUILD_TESTS``

* Documentation explaining the new requirement that all PRs must include a file
  explaining what the PR does. (PR#545)
* Template for the news directory. (PR#545)
* Additional boundary condition options in the dagmcMetaData class (PR#690)


**Changed:**

* replaced calls for overlap tolerance and numerical precision in DagMC constructor (PR#708)
* Packages now use direct path to include one-another, instead of relative ones (PR#703)
* No longer require Fortran compiler unless building MCNP5/6 (PR#701)
* Update amalgamated PyNE version to v0.7.3 (PR#700)
* revamped Material management to leverage the PyNE::MaterialLibrary in place of the map<string, PyNE::Material> (PR#700)
* Add DAGMC guard around (PR#695):

  * changed if statement in history_neutral_high.F90
  * check for goto statement in charged_particle_history.F90

* Adding optional double-down dependency to enable ray tracing with Embree. (PR#693)
* Replacing Travis in favor of CircleCI (PR#692, PR#698)
* Splitting up the docker container building process into multiple files (PR#692, PR#697)
* dagmcMetaData (PR#688, PR#689, PR#690):

  * Behavior to ignore missing density assignments for more flexible integration with certain codes. (PR#688)
  * Updates to the coding style. (PR#689)
  * Allows boundary condition values, graveyard material assignments, and vacuum material assignments to be lowercase

* removed LAPACK dependency; replaced with Eigen3 for DAGMC (PR#686) and MOAB (PR#683) 
* Enabling testing for the shared object build of DAGMC (PR#674)
* Adding RPATH value for our build of Geant4 on CI (PR#674)
* Including additional test output on failure in CI (PR#674)
* PullRequest-Agent suggestionsmcnp (PR#665)

  * using std::err for errors
  * update to C++11 standards for converting ints to strings
  * removed unnecessary comments
  * moved Graveyard and Vacuum strings to variables

* PullRequest-Agent suggestions MakeWaterTight (PR#666):

  * remove commented code blocks that are either outdated or are debug statements
  * improvements to some logic for clarity
  * use of standard library containers to avoid potential memory leaks in Arc.cpp/Gen.cpp
  * improvements to struct/variable names
  * declared variables for "magic numbers"
  * passing by const reference where possible to avoid unnecessary memory allocation
  * removed an unused function (Arc::create_loops_from_oriented_edges_fast)

* PullRequest-Agent suggestions DagMC (PR#671, PR#676):

  * updated pointer management to RAII ("Resource Allocation Is Initialization") technique, MBI is now a shared_ptr unless passed 
    as a raw pointer in the DagMC constructor (can be returned as a shared_ptr if not provided as a raw pointer),  GTT is now a 
    shared_ptr, and can only be returned as such, GQT is now a uniq_ptr, (and can't be return - not change there)
  * tests: DagMC instance is now a shared_ptr, when used, MBI instance is now a shared_ptr

* PullRequest-Agent suggestions uwuw, tally, overlap_check, build_obb, misc/tests (PR#680)
* PullRequest-Agent suggestions Geat4 (PR#691)
* Updates to variable names in make_watertight files (PR#672)
* Changed name of overlap_check executable directory from "build" to
  "app". (PR#653)
* all directories named `build` are changed to `app` for clarity. (PR#645)
* ``dagmc/src/make_watertight``: now accepting output_filename. (PR#636)
* ``dagmc/src/check_watertight``: now accepting output_filename. (PR#636)
* Have the update_pyne script copy over the source.F90 files in pyne for MCNP5
  and MCNP6. (PR#626)
* Update amalgamated pyne. (PR#626)
* Updated amalgamated pyne. (PR#617)
* The `ASTYLE_ONLY` Travis variable has been replaced with a `HOUSEKEEPING_ONLY`
  variable. If this variable is on, DAGMC will not be built and it will instead
  only perform 3 housekeeping checks: (PR#610)

  * News file: the CI will fail if a news file with the correct filename is not
    included.

  * Astyle: the version of astyle we use on the CI has been upgraded to 3.1.
    This is the version that is default on Ubuntu 18.04.
  * Documentation: the CI will now attempt to build the DAGMC documentation and
    will fail if it finds any errors or warnings.

* The dockerfile has been modified so that it can be built with both Ubuntu
  16.04 and 18.04. (PR#610)
* The docker images have been moved from the cnerg dockerhub organization to the
  svalinn organization. (PR#610)
* The new build matrix for the non-housekeeping run is 2x2x2: (PR#610)

  * Ubuntu 16.04 vs. 18.04
  * gcc vs. clang
  * gcc-5.3 on 16.04; gcc-7.3 on 18.04
  * clang-3.8 on 16.04; clang-6.0 on 18.04
  * MOAB 5.1.0 vs. master vs. develop

* The builds that use MOAB master and develop are allowed to fail without the
  entire CI failing. The CI will show as having passed once the housekeeping
  build and the four MOAB 5.1.0 builds have passed. (PR#610)

* The CI will only build against MOAB master and develop during non-pull request
  builds; i.e. only during push builds and nightlies. (PR#610)

* MOAB 5.1.0 is now included in the docker image so it does not need to be built
  every time the CI is run. (PR#610)

  * This is to save time, since we expect that previous versions of MOAB will
    not change. If it does change, we can update the Docker images.
  * MOAB master is still built every time it is needed.

* MOAB is now built with pymoab support. This is for future-proofing in case
  DAGMC ever needs access to this functionality. (PR#610)
* MOAB is now built against both custom-built HDF5 (1.10.4, up from 1.8.13) and
  against system HDF5. (PR#610)

  * The MOAB built against system HDF5 is currently unused, however, as there is
    currently a bug that makes it so DAGMC cannot build static executables if
    using system HDF5. If/when this bug is fixed, then building DAGMC with
    system HDF5 can be added to the build matrix.

* Geant4 has been upgraded to version 10.5. (PR#610)
* Building the documentation will throw an error if it encounters any warnings
  or errors. The previous warnings and errors that were occurring have been
  fixed. (PR#610)
* Throw a fatal error if trying to build static executables but not static
  libraries, or shared executables but not shared libraries. (PR#605)
* Added measure and source_sampling to amalgamated pyne and removed the
  standalone files we used to use. (PR#604)
* Move keyword type to FC card in the document doc/userguide/tally.rst.
  (PR#600)
* A small change to a single line of the dag-mcnp model prep file. (PR#599)
* ``CMakeLists.txt`` (PR#597)
* ``src/mcnp/meshtal_funcs.cpp`` (PR#597)
* ``src/tally/KDEMeshTally.cpp`` (PR#597)
* ``src/tally/MeshTally.cpp`` (PR#597)
* ``src/tally/MeshTally.hpp`` (PR#597)
* ``src/tally/MeshTally.hpp`` (PR#597)
* ``src/tally/TallyData.cpp`` (PR#597)
* ``src/tally/TrackLengthMeshTally.cpp`` (PR#597)
* CMake commands for linking all DAGMC libraries s.t. they are added to the exported targets. (PR#662)
* Updated amalgamated pyne to match the main pyne repo. (PR#595)
* Travis CI no longer attempts to build DAGMC against moab master. (PR#584)
* When configuring MPI-enabled DAG-MCNP6, do not rely on
  ``MPI_Fortran_INCLUDE_PATH`` being set because this variable is not set when
  using CMake 3.10 or newer. Instead, use ``MPI_Fortran_COMPILER``. (PR#579)
* Use the values of ``MOAB_INCLUDE_DIRS`` and ``MOAB_LIBRARY_DIRS`` from
  ``MOABConfig.cmake`` instead of trying to determine them ourselves. Note that
  this change makes DAGMC incompatible with MOAB 5.0. (PR#578)
* Use MOAB 5.1.0 on CI instead of 5.0. (PR#578)
* CMakeFile for DAG-MCNP6 to accomodate MCNP6.2. (PR#569)
* Use bind(c) in fmesh_mod.F90 to avoid the need for name mangling on the C++
  side. (PR#556)
* Rename MCNP patch files to mcnpXXX.patch, where XXX is the version turned
  into a 3-digit number. (PR#556)
* Change pretty much every ``CMakeLists.txt`` file in the entire repo to use the
  new macros. Almost all the cmake files got much shorter because of this
  change. (PR#555)
* Change how we find HDF5. Previously, HDF5 was required to be in users'
  ``$PATH``. Now, the location of HDF5 is determined automatically by reading
  variables from ``MOABConfig.cmake``. (PR#555)
* Change how we find MOAB. Previously, MOAB was required to be in users'
  ``$LD_LIBRARY_PATH``. Now, users must specify ``-DMOAB_DIR`` when running
  cmake. (PR#555)

  * Note that the ``MOABConfig.cmake`` file is no longer used to find any MOAB
    files.

* Since users no longer need to change their ``$PATH`` or ``$LD_LIBRARY_PATH``,
  remove the changes to those variables in the CI scripts. (PR#555)
* Rename the cmake commands used to build DAG-MCNP5/6 with plotting and MPI
  support. The new commands are ``BUILD_MCNP_PLOT`` and ``BUILD_MCNP_MPI``.
  (PR#555)
* Rename the cmake command used to build static executables from
  ``BUILD_STATIC`` to ``BUILD_STATIC_EXE``. The old name was confusing because
  the option only controls the linking of executables, while libraries are
  always built both static and dynamic. (PR#555)
* Rename the ``test`` folders in ``src/dagmc`` and ``src/mcnp`` to ``tests`` to
  conform with other unit test directories. (PR#555)
* Move the source files for the make_watertight and uwuw_preproc executables
  into a new ``build`` directory, keeping the source files for the library where
  they are. This conforms with other DAGMC features that have both a library and
  an executable. (PR#555)
* Replace the mcnpfuncs internal library with an object library. (PR#555)
* For the pyne_dagmc library, only use ``-O0`` optimzation when building with
  Intel C++. (PR#555)
* Update documentation to reflect all changes. (PR#555)
* Moved all source code into the ``src`` directory. (PR#552)
* Fix download link to astyle 3.0.1 .deb file. (PR#549)
* Direct Travis to grab the docker image from the cnerg dockerhub account
  instead of Lucas's account. (PR#546)

**Deprecated:**

* DagMC: Deprecated constructor using a raw pointer for the MBI instance,
  prefered way uses shared_ptr for MBI instance. (PR#671)
* `DagMC::interface_revision` method (PR#707)

**Removed:**

* Remove the ``FindHDF5.cmake`` file as it is no longer needed. (PR#555)
* ``gtest/README`` and ``gtest/configure.sh``: no longer used; last commit in
  March 2014. (PR#544)
* ``tools/build/*``: no longer used; last commit in June 2014. (PR#544)
* ``cmake/FindPyne.cmake``: no longer used; last commit in June 2014. (PR#544)
* ``tools/finish_dagmc_geom*``: out of date; last commit in June 2014. (PR#544)
* ``tools/txcorp_bld/*``: no idea what this is; last commit in June 2014. (PR#544)
* ``tools/dagmc_tag_eg/*``: out of date; last commit in October 2014. (PR#544)
* ``tally/tools/boundary_correction/*``: broken; last commit in June 2016. (PR#544)

**Fixed:**

* pyne::MaterialLibrary calls in fluka_funcs.
* Tally::particle_name to Tally::particle_names update in fluka_funcs as well.
* eigen3:

  * remove bad flag in MOAB build (PR#684)
  * fixed use include directories (PR#694)

* Regenerate the DAGMC_LIBRARIES variable upon re-running cmake. (PR#643)
* Fix error in documentation where cmake was not pointing to the DAGMC source
  dir as it should. (PR#632)
* Updated links to OpenMC documentation. (PR#630)
* Make the MW_REG_TEST_MODELS_URL variable available to the docker image. (PR#621)
* The `make_watertight_regression_tests` should now be run if the CI is not
  doing a PR build. (PR#610)

  * I believe this was broken for an undetermined amount of time; I do not
    believe they were ever getting run regardless of whether the CI was doing a
    PR build or not. This is because intrinsic Travis variables like
    `$TRAVIS_PULL_REQUEST` are only available to `.travis.yml`; if they are
    needed in other scripts, they need to be passed manually, and this was not
    happening before.

* Fixes issue with unstructured mesh tallies. (PR#597)
* Now produces a vector tag of size num_groups instead of num_groups+2 scalar
  tags. (PR#597)
* Also produces a total tally tag. (PR#597)

**Security:** None
