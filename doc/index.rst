..  DAGMC master file, created by
    sphinx-quickstart on Fri Aug 31 10:08:00 2012.
    You can adapt this file completely to your liking, but it should at least
    contain the root `toctree` directive.

DAGMC: Direct Accelerated Geometry Monte Carlo
==============================================

..  raw:: html
    :file: slideshow.html

|

..  image:: https://travis-ci.org/svalinn/DAGMC.svg?branch=develop
    :target: https://travis-ci.org/svalinn/DAGMC
    :height: 20
    :width: 90
    :align: left

|

Direct Accelerated Geometry Monte Carlo (DAGMC) is a software package that
allows users to perform Monte Carlo radiation transport directly on CAD models.

DAGMC has been integrated into a variety of Monte Carlo radiation codes
including MCNP5_, Geant4_, FLUKA_, Tripoli4_, and Shift_. There are also efforts
planned to integrate DAGMC into other codes such as MCNP6_, Serpent2_,
OpenMC_, Phits_, and Frensie.

DAGMC currently relies on using the solid modeling software Cubit_ or its
commercial counterpart, Trelis_, to prepare solid models. These packages can be
used to import CAD models from other tools such as SolidWorks, CATIA, etc., or
to create geometry from scratch. DAGMC also relies on Trelis/Cubit to assign
materials and other geometry-related information.

..  toctree::
    :maxdepth: 2

    install/index
    usersguide/index
    contribute/index

..  toctree::
    :hidden:

    gallery/gallery

..  _CNERG: http://cnerg.github.io
..  _MOAB: http://press3.mcs.anl.gov/sigma/moab-library
..  _Cubit: https://cubit.sandia.gov
..  _Trelis: http://www.csimsoft.com/trelis
..  _MCNP5: https://mcnp.lanl.gov
..  _MCNP6: https://mcnp.lanl.gov
..  _Geant4: http://geant4.cern.ch
..  _FLUKA: http://www.fluka.org/fluka.php
..  _Tripoli4: https://rsicc.ornl.gov/codes/ccc/ccc8/ccc-806.html
..  _Shift: http://web.ornl.gov/sci/nsed/rnsd/rt
..  _Serpent2: http://montecarlo.vtt.fi
..  _OpenMC: https://mit-crpg.github.io/openmc
..  _Phits: http://phits.jaea.go.jp
