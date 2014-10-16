#!/usr/bin/python

"""
Table of particles used in PyNE:

        =========  =========
        particle   name (z)
        =========  =========
        neutron    n
        proton     p
        deuterium  d
        tritium    t
        Helium-3   He3
        alpha      a
        gamma      gamma
        =========  =========
"""

"""
Define: pyne_particles as a disctionary of particles used in PyNE; 'keys' being PyNE convention for the particle and 'values' being possible names of the particle to be used on CAD group names.
"""

pyne_particles={'n':['Neutron', 'neutron', 'N', 'n'], 'p':['Proton','proton','P','p'], 'd':['Deuterium','deuterium','D','d'], 't':['Tritium','tritium','T','t'], 'He3':['Helium-3','Helium_3','helium-3','helium_3','He3','he3'], 'a':['Alpha','alpha','A','a'], 'gamma':['Gamma','gamma','G','g']}
         
