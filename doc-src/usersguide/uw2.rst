The University of Wisconsin Unified Workflow
----------------------------------------

The University of Wisconsin Unified Workflow (UW<sup>2</sup>) aims to solve the 
issue of running the same Monte Carlo problem using mutiple physics codes. Currently,
if you wish to run the same problem in multiple codes you must fully recreate the
input deck for each code, or maybe write a full syntax translator. UW<sup>2</sup>, allows users
to tag or associate groups of volumes or surfaces with a simple human readable syntax
that is translated and stored in the geometry file of a DAGMC problem.

The workflow uses the Python for Nuclear Engineering toolkit `PyNE <http://pyne.io>`_, we 
levereage the existing infrastructure in PyNE to allow a consistent transport problem to be
defined across all MC codes.

Materials
+++++++++++++++++++++++++++++++++++++++

Materials are the most painful transfer from code to code, since each MC code 
specifies materials in a different way. Instead, we tag groups of volumes
with a name and syntax that corresponds to material compositions in a predefined
material library.

The syntax for describing materials is:
::
     %> group "mat:<Name of Material>"

If you wish to specify a material that is present in the library with a different 
density, then 
::
     %> group "mat:<Name of Material>/den:<density>"

So for example, to specify a Stainless Steel at a density of 12.0 g/cc,
::
     %> group "mat:Stainless Steel/den:12.0"


Scoring
+++++++

Each MC code implements tallies or scores, in very specific ways such that there
is sometimes no equivlent to a tally you may be familiar with, code to code. However
, there is a syntax to allow you to request scores on geomemtric elments, for example,
::
     %> group "tally:flux/neutron"

or,
::
     %> group "tally:photon/current"

Using the underlying PyNE libraries, we can write out the appropriate MC code tally specification
snippet, this allows the number of codes the DAGMC supports to grow organically with those that
PyNE supports.

Source Defintion and other problem specific data
+++++
UW<sup>2</sup> is a work in progress and we would ultimately like to import all problem specific data into 
the geometry file along with material and tally metadata.
