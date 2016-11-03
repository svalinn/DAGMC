User's guide
============

..  |UW2| replace:: UW\ :sup:`2`

This document is intended for users of the DAGMC toolkit who want to better
understand all the potential workflows for each Monte Carlo code. UW--Madison
has developed several workflows for the Monte Carlo codes it supports. Exactly
which workflow suits your needs depends on which codes you expect to use and how
much manual data transfer you wish to do.

+------------------------+---------------------+--------------------+
| **Feature**            | **|UW2| workflow**  | **Basic workflow** |
+------------------------+---------------------+--------------------+
| Assignment of metadata | Manual              |  Manual            |
+------------------------+---------------------+--------------------+
| MCNP5 tallies          | Automatic           |  Manual            |
+------------------------+---------------------+--------------------+
| MCNP5 materials        | Automatic           |  Manual            |
+------------------------+---------------------+--------------------+
| Fluka tallies          | Automatic           |  Manual            |
+------------------------+---------------------+--------------------+
| Fluka materials        | Automatic           |  Manual            |
+------------------------+---------------------+--------------------+
| Geant4 tallies         | Automatic           |  Manual            |
+------------------------+---------------------+--------------------+
| Geant4 materials       | Automatic           |  Manual            |
+------------------------+---------------------+--------------------+
| Tripoli4 tallies       | Not supported       |  None              |
+------------------------+---------------------+--------------------+
| Tripoli4 materials     | Not supported       |  None              |
+------------------------+---------------------+--------------------+

If you are interested in running your geometry in multiple physics packages with
minimal effort, the |UW2| workflow will be the best option as it will automate
much of the work. However, if you are interested in only a specific physics
engine, one of the basic workflows may interest you more. We still would
recommend the |UW2| workflow as it offers robust and tested methods for material
creation.

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
    UnderstandingRayHistoryState
