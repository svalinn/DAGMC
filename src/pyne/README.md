Amalgamated PyNE Build
==========================================================
The DAGMC toolkit now includes an amalgamated build of PyNE version 0.5. This
helps the building of DAGMC on clusters and other systems that don't necessarily
need a full installation of PyNE.

Updating the PyNE Build
==========================================================
You will need to generate the atomic and decay data from PyNE. Issue the
following commands from the pyne/src directory:

    python decaygen.py
    python atomicgen.py

To generate the PyNE amalgamated source files, we use the following command:

    ./amalgamate.py -s pyne.cpp -i pyne.h -f license.txt src/utils.* \
     src/extra_types.h src/h5wrap.h src/state_map.cpp src/nucname.*  \
     src/rxname.* src/particle.* src/data.* src/json-forwards.h      \
     src/json.h src/jsoncpp.cpp src/jsoncustomwriter.h               \
     src/jsoncustomwriter.cpp src/material.* src/tally.* src/decay.* \
     src/atomic_data.*

The above represents the smallest build that is trivially compilable without
errors.

At the time of writing, the produced file ``pyne.cpp`` will need to be edited
slightly in order to make DAGMC compilable. You must comment out the line in
``pyne.cpp`` that says

    #include "particle.h"

and you must also replace the code inside the function
``pyne::Material pyne::Material::decay(double t)`` with

    pyne::Material pyne::Material::decay(double t)
    {
      Material rtn;
      std::cout << "--Warning--Warning--Warning--Warning--Warning--Warning--" << std::endl;
      std::cout << "  There is no decay function in the material object within" << std::endl;
      std::cout << "  this amalgamated pyne build" << std::endl;
      std::cout << "--Warning--Warning--Warning--Warning--Warning--Warning--" << std::endl;
      return rtn;
    }

User update of PyNE
===========================================================
If a DAGMC user wishes to update the version of PyNE associated with DAGMC, they
need only update the the amalgamated build by including the additional source
files they need relative to the above list.
