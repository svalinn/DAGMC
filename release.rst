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

1. The DAGMC source code repository has a tag named `stable <stable_>`_ that
   points to the same SHA as the newest release.
2. The PyNE CI docker images rely on the DAGMC `stable` source code tag (this was incorporated in PyNE PR #1415)

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
3. Update DAGMC `stable <stable_>`_ tag to future release hash = DAGMC RC branch HEAD
   
   `git tag -d stable; git tag stable; git push -f --tags upstream`
4. Manually invoke `Github action to rebuild PyNE Docker images
   <https://github.com/pyne/pyne/actions/workflows/docker_publish.yml>`_.  Note
   that this will automatically build & test PyNE with the updated DAGMC RC that
   is pointed to by the `stable <stable_>`_ tag
5. Publish PyNE release - DAGMC depends on a specific hash of PyNE, while PyNE
   depends on a tag of DAGMC
    a. Merge PyNE RC branch to `develop`
    b. Create release from `develop`
6. Update PyNE Submodule in DAGMC after merge of PyNE RC
7. Publish the DAGMC release
    a. Merge DAGMC RC branch to `develop`
    b. Create release from `develop`

.. _stable: https://github.com/svalinn/DAGMC/releases/tag/stable