Amalgamated PyNE Build
==========================================================
The DAGMC toolkit now includes an amalgamated build of PyNE version 0.5, this helps
building of DAGMC on clusters and other systems that don't need nessesarilly need 
a full installation of PyNE.

Updating the PyNE Build
==========================================================
You need to generate the atomic and decay data from PyNE, issue the following commands from
the pyne/src directory

     python decaygen.py
     python atomicgen.py

To generate the PyNE amalgamated source files, we use the following command:

    ./amalgamate.py -s pyne.cpp -i pyne.h -f license.txt src/utils.* \
     src/extra_types.h src/h5wrap.h src/state_map.cpp src/nucname.*   \
     src/rxname.* src/particle.* src/data.* src/json-forwards.h       \
     src/json.h src/jsoncpp.cpp src/jsoncustomwriter.h                \
     src/jsoncustomwriter.cpp src/material.* src/enrichment_cascade.* \
     src/enrichment.* src/enrichment_symbolic.h                       \
     src/enrichment_symbolic20.cpp src/tally.* src/decay.* src/atomic_data.*

At the time of writing due to a bug in PyNE the produced source files will not
compile directly you must currently modifiy the line in pyne.cpp that says
``#include "particle.h"''

User update of PyNE
===========================================================
If a DAGMC users wishes to update the version of PyNE associated with 
DAGMC, they need only update the the amalgamated build by including the additional
source files they need relative to the above.


