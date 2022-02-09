Instructions for Creating a DAGMC Release
============================================

Because of a complex dependency between DAGMC and PyNE the following procedures
are recommended for navigating those dependencies for a new release.

DAGMC depends on PyNE only as a submodule from which C++ files are directly
incorporated into the build.  It may be desirable to perform a patch-level
release of DAGMC when any of the C++ files on PyNE change.

PyNE depends on DAGMC as a third-party library that must be installed prior to
building PyNE.  Changes in DAGMC do not necessarily require updates to PyNE, 
so may not result in a patch-level release of PyNE.

It is often desirable to perform patch-level releases of PyNE and DAGMC on
concert.  If so, the following procedure should reduce conflicts in the
dependencies.

Assumptions
------------

1. The DAGMC repository has a `stable` tag that points to the same SHA as the newest release.
2. The PyNE CI docker images rely on the DAGMC `stable` tag (this was incorporated in PyNE PR #1415)

Procedure
----------

1. Prepare PyNE for release
    a. Create release candidate (RC) branch in pyne/pyne
    b. Create a pull request (PR) *into the PyNE RC branch* in order to
        i. Change the version number in `pyne/pyne_version.py`
        ii. Update `CHANGELOG.txt`
    c. Merge that PR when approved
    d. Potentially perform other updates to the RC branch via PRs
    e. *DO NOT MERGE* RC branch
2. Prepare DAGMC for release
    a. Create RC branch in svalinn/DAGMC
    b. Create PR into DAGMC RC branch
        i. Update PyNE Submodule in DAGMC after all changes are made to PyNE RC
           branch. The PyNE RC branch is still not merged at this point.
        ii. Update the DAGMC version number in `CMakeLists.txt`
        iii. Update DAGMC's `CHANGELOG.rst`
    c. Merge that PR when approved
    d. Potentially perform other updates to the RC branch via PRs
    e. *DO NOT MERGE* RC branch
3. Update DAGMC `stable` tag to future release hash = DAGMC RC branch HEAD
4. Manually invoke Github action to rebuild PyNE Docker images.  Note that this
   will automatically build & test PyNE with the updated DAGMC RC that is
   pointed to by the `stable` tag
5. Publish the DAGMC release
    a. Merge DAGMC RC branch to develop
    b. Create reelase from develop
6. Publish PyNE release 
    a. Merge PyNE RC branch to develop
    b. Create reelase from develop
