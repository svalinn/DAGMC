Building the DAGMC website
==========================

The `DAGMC website <DAGMC_>`_ is built using the Sphinx_ documentation
generator.

DAGMC uses a two-branch system to maintain a clean process of rebuilding the
website. The ``develop`` branch contains the source restructured text documents
and Sphinx configuration used to build the site. The ``gh-pages`` branch
contains the processed and published web content that is produced by Sphinx
using the source files in the ``develop`` branch. The files in this branch
should NOT be edited directly.

The DAGMC documentation build system relies on the ``Makefile`` located in the
top-level directory of the DAGMC repository. Here is a summary of the available
commands:

``make help``: Display available options and exit.

``make html``: Build the documentation for viewing on a local machine.

``make clean``: Remove the locally-built documentation.

``make publish``: Publish the documentation located in the ``develop`` branch to
the ``gh-pages`` branch. To prevent a situation where the wrong branch is used
to build the documentation, the git remote ``origin`` should be the main
repository and not a fork. Additionally, the branch used for building the
documentation should not contain any additional changes not present on Github.
In other words, in order to use ``make publish``, the result of
``git remote -v && git status`` should be
::

    origin  https://github.com/svalinn/DAGMC (fetch)
    origin  https://github.com/svalinn/DAGMC (push)
    On branch develop
    Your branch is up-to-date with 'origin/develop'.
    nothing to commit, working tree clean

..  _DAGMC: https://svalinn.github.io/DAGMC
..  _Sphinx: https://www.sphinx-doc.org
