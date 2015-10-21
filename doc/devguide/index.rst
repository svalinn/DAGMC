DAGMC Developers Guide
==========
DAGMC is an open source project that facilitates the running of CAD based particle transport problems.

Contributing
~~~~~~~~~~~~
Contributing to the DAGMC project is very straightforward, DAGMC is hosted on Github where issues and pull requests are 
discussed and merged. We use the git version control system which could be the most unfamiliar aspect of contributing for 
most people. The general workflow to contribute to DAGMC and many other open source projects involves steps like below.

.. image:: workflow.png
   :height: 300
   :width:  600
   :alt:    Image showing the github workflow

There are 5 main steps:

  1) Forking
  2) Cloning
  3) Branching
  4) Pushing
  5) Pull Requesting

These stages are outlined below.

Forking
--------
To start the repository must be forked. The default branch that is `checked out` is the develop branch, 
which contains all the most upto-date changes. If you wish to 
develop you should make a fork of the DAGMC repository in your Github account. The easiest way to do this is to click on the 
`fork` button from the `svalinn/dagmc` branch shown below.

.. image:: workflow_fork.png
   :height: 150
   :width:  600
   :alt:    Image showing how to fork a repo

This fork will be an exact snapshot of the `svalinn/dagmc` repository at the time you clicked `fork`. Any new features
that you wish to develop should be based from the develop branch of this repository, unless you know exactly what you're 
doing. 

Cloning
---------
You should now clone your fork of this repository to your local machine
::
   prompt %> git clone https://github.com/githubusername/dagmc
   prompt %> cd dagmc


Branching
---------
The base level of the repository contains folders for each of the supported monte carlo codes, the tools directory, and our
amalgamated PyNE build. First, you need to checkout a new branch in which to keep your changes.
::
   prompt %> git checkout -b "my_feature_branch"

To insert a new feature edit an existing file or add new ones as required, remember to update the 
CMakeList.txt files as required. Your new changes need to be added, commited and pushed.
::
   prompt %> git add <files needed to add>
   prompt %> git commit -m "This is a message that describes why we need these changes"

Pushing
---------------
Now that your changes are commited, you push the changes to your remote branch in your clone of DAGMC
::
   prompt %> git push myrepo my_feature_branch

Pull Requesting
----------------
If you immediately go to your fork on Github you should then see a message like that shown below, if you click on the green button
you will create a pull request against the develop branch of DAGMC. If you've waited a few tens of minutes between pushing and 
going to Github you may have to manually create a pull request. Your pull request will launch our continous integration tests and
at some point in the near future your changes will pass all the unit tests or indeed may break the tests.

Build System
~~~~~~~~~~
We exclusively use CMake as the build system for our tools, we have tried to adapt a fairly modular system where variables
are locally scoped where possible, the extent to which you change the build system as a developer will depend on the extent of 
your changes, for example adding new tests will require small changes but adding the support of another Monte Carlo code will
be more complex. 

Testing & Continuous Integration
~~~~~~~~~~

We use the `Google Test <https://code.google.com/p/googletest/>`_ gtest libraries to control testing of our code and we 
use the `Travis <https://travis-ci.org/>`_ continuous integration system to test all changes to the code. When you add 
features to the codebase, tests should always be added which prove the capabilities that have been added. 

When a developer
makes a pull request on GitHub, Travis detects this change and launches the build as specified in the .travis.yml file. Travis
pulls your feature branch, the MOAB libraries, HDF5, etc as required and then launches the tests. Each test is run in succession 
and failure is reported if any dependency fails to build or if any test fails, an example of a Travis report is shown below

.. image:: travis_example.png
   :height: 400
   :width:  600
   :alt:    Image showing the status of the an example Travis-CI run

Once the testing is complete and your changes have been verified not break any of the existing capabilities, a reviewer will check your pull request over and may suggest some modifications to meet the C++ style, good practice and then will approve or reject your
pull request. 

General Style
~~~~~~~~~~
Explicit namespacing is preferred, so rather than using the `using namespace xxx` command, you should prefix the variable with the
class name, i.e.
::
   pyne::Material new_material; // this is a new material

is preferred over, 
::
  using namepspace pyne;
  Material new_material; // this is a new material

C++ Style
~~~~~~~~~~

We conform to the C++ style guide, we have included a C++ style guide formatter to make a developers life much easier. When you
have added all the features you want to add, the style guide formatter should be run,
::
   prompt %> astyle --style=linux --indent=spaces=2

Then commit the changes.
