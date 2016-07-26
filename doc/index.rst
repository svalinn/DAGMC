.. DAGMC master file, created by
   sphinx-quickstart on Fri Aug 31 10:08:00 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

What is DAGMC?
======

The Direct Accelerated Geometry Monte Carlo (DAGMC) library, a 
component of the Mesh-Oriented datABase [MOAB_], provides
fundamental functions for the ray-tracing and related geometry
operations of Monte Carlo radiation transport directly on complex 3-D
geometries such as those created in modern solid modeling software.

DAGMC is designed in a modular fashion with the expectation that it
can be integrated into a variety of Monte Carlo radiation tools. The
CNERG_ group at UW-Madison has historically focussed development on the
MCNP5_ software developed at `Los Alamos National Laboratory
<http://www.lanl.gov>`_ and distributed by the `Radiation Safety
Information Computing Center <http://rsicc.ornl.gov>`_. However, recently
DAGMC has been integrated into the following Monte Carlo physics packages
: MCNP5_, Tripoli4_, Fluka_, Geant4_, and Shift_.

We have prior experience integrating DAGMC with MCNP/X, and planned
efforts to integrate DAGMC with other Monte Carlo physics packages
including: MCNP6_, Serpent2_, Phits_, OpenMC_ and Frensie.

While we don't have a complete GUI, we currently rely on the Trelis_ software
from Csimsoft or the Cubit_ software from Sandia.  It plays a role in our workflow 
that can include importing CAD-files from other tools such as SolidWorks,
CATIA, etc.  In addition to defining the geometry, we rely on
Cubit for material assignment and can also support some other aspects
of input definition that are tied to the geometry (tallies and
variance reduction parameters, for example).  Some knowledge of the 
specific Monte Carlo application is necessary for other 
parameters such as material definition, run control and source definition.

MOAB_ is developed by a team of collaborators at Argonne National Laboratory (ANL) and the 
University of Wisconsin-Madison (UW), who are also working on improving some of 
the GUI tools available for manipulating workflows.


.. toctree::
   :maxdepth: 1

   usersguide/index
   CNERG Support for DAGMC
   tools/index
   devguide/index
   gallery/gallery
   DagmcPublications
   upcoming

.. _MOAB: http://sigma.mcs.anl.gov/moab-library/
.. _MCNP5: http://mcnp.lanl.gov
.. _MCNP6: http://mcnp.lanl.gov/
.. _Cubit: http://cubit.sandia.gov
.. _Trelis: https://http://www.csimsoft.com/trelis.jsp
.. _CNERG: http://cnerg.engr.wisc.edu
.. _Fluka: http://www.fluka.org/fluka.php
.. _Geant4: http://geant4.cern.ch/
.. _Tripoli4: https://rsicc.ornl.gov/codes/ccc/ccc8/ccc-806.html
.. _Shift: http://web.ornl.gov/sci/nsed/rnsd/rt/code.shtml
.. _Serpent2: http://montecarlo.vtt.fi/development.htm
.. _OpenMC: https://mit-crpg.github.io/openmc/
.. _Phits: http://phits.jaea.go.jp/
