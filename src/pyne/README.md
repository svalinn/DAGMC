Amalgamated PyNE Build
======================
The DAGMC toolkit now includes an amalgamated build of PyNE version 0.5. This
helps the building of DAGMC on clusters and other systems that don't necessarily
need a full installation of PyNE.

Updating the PyNE Build
=======================
You will need to generate the atomic data from PyNE. Issue the following
command from the pyne/src directory:

    python atomicgen.py

To generate the PyNE amalgamated source files, we use the following command:

    ./amalgamate.py -f license.txt src/utils.* src/extra_types.h src/h5wrap.h \
    src/state_map.cpp src/nucname.* src/rxname.* src/particle.* src/data.* \
    src/json-forwards.h src/json.h src/jsoncpp.cpp src/jsoncustomwriter.* \
    src/material.* src/tally.* src/atomic_data.*

The above represents the smallest build that is trivially compilable without
errors, with the caveat that there are two functions that must be edited. These
functions are ``pyne::Material::decay``, which should be replaced with

    pyne::Material pyne::Material::decay(double t) {
      throw pyne::ValueError("Material::decay is not supported in this amalgamated"
                             "version of PyNE.");
    }

and ``pyne::Material::cram``, which should be replaced with

    pyne::Material pyne::Material::cram(std::vector<double> A,
                                        const int order) {
      throw pyne::ValueError("Material::cram is not supported in this amalgamated"
                             "version of PyNE.");
    }

User update of PyNE
===================
If a DAGMC user wishes to update the version of PyNE associated with DAGMC, they
need only update the the amalgamated build by including the additional source
files they need relative to the above list.
