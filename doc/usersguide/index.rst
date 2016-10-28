User's guide
============

..  |UW2| replace:: UW\ :sup:`2`

This document is intended for users of the DAGMC toolkit who want
to better understand all the potential workflows for each Monte Carlo
code. UW Madison have developed several workflows and exactly which 
workflow suits your needs depends on which codes you expect to use
and how much manual data transfer you wish to do.

+---------------------+----------------+----------------+
| Feature             | Basic Workflow | UW2 Workflow   | 
+---------------------+----------------+----------------+
|  MCNP5 Tallies      | Manual         |  Automated     |
+---------------------+----------------+----------------+
|  MCNP5 Materials    | Manual         |  Automated     |
+---------------------+----------------+----------------+
|  Fluka Tallies      | Manual         |  Automated     |
+---------------------+----------------+----------------+
|  Fluka Materials    | Manual         |  Automated     |
+---------------------+----------------+----------------+
|  Geant4 Tallies     | None           |  Automated     |
+---------------------+----------------+----------------+
|  Geant4 Materials   | None           |  Automated     |
+---------------------+----------------+----------------+
|  Tripoli4 Tallies   | Manual         |  None          |
+---------------------+----------------+----------------+
|  Tripoli4 Materials | Manual         |  None          |
+---------------------+----------------+----------------+

So, if you are interested in running your geometry in multiple physics 
packages with minimal effort, the UW2 workflow is your best option
as it does the most work for you, however if you are interested
in only a specific physics engine one of the basic workflows may
interest you more. We still would recommend the UW2 workflow as it 
offers robust and tested way of material creation.

..  toctree::
    :maxdepth: 1

    trelis_basics
    uw2
    codes/index
    trelis_workflow
    tally
    tools
    postprocessing

..  toctree::
    :hidden:

    mcnp2cad
