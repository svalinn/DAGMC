..  DAGMC master file, created by
    sphinx-quickstart on Fri Aug 31 10:08:00 2012.
    You can adapt this file completely to your liking, but it should at least
    contain the root `toctree` directive.

DAGMC: Direct Accelerated Geometry Monte Carlo
==============================================

..  raw:: html
    :file: slideshow.html

|

..  image:: https://circleci.com/gh/svalinn/DAGMC.svg?style=shield
    :target: https://circleci.com/gh/svalinn/DAGMC
    :height: 20
    :width: 90
    :align: left

|

Direct Accelerated Geometry Monte Carlo (DAGMC) is a software package that
allows users to perform Monte Carlo radiation transport directly on CAD models.

DAGMC has been integrated into a variety of Monte Carlo radiation codes
including MCNP5_, MCNP6_, Geant4_, FLUKA_, Tripoli4_, OpenMC_, and Shift_. There are also
efforts planned to integrate DAGMC into other codes such as Serpent2_, Phits_, and FRENSIE_.

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
    CHANGELOG

..  toctree::
    :hidden:

    gallery/gallery

..  _CNERG: https://cnerg.github.io
..  _MOAB: https://press3.mcs.anl.gov/sigma/moab-library
..  _Cubit: https://cubit.sandia.gov
..  _Trelis: https://www.csimsoft.com/trelis
..  _MCNP5: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/mcnp5.shtml
..  _MCNP6: https://mcnp.lanl.gov
..  _Geant4: https://geant4.cern.ch
..  _FLUKA: http://www.fluka.org/fluka.php
..  _Tripoli4: https://rsicc.ornl.gov/codes/ccc/ccc8/ccc-806.html
..  _Shift: https://meitner.ornl.gov/doe-codes/shift
..  _Serpent2: http://montecarlo.vtt.fi
..  _OpenMC: https://docs.openmc.org/en/latest/index.html
..  _Phits: https://phits.jaea.go.jp
..  _FRENSIE: https://github.com/FRENSIE/FRENSIE
