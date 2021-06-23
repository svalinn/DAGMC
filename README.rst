DAGMC: Direct Accelerated Geometry Monte Carlo
==============================================

..  image:: https://circleci.com/gh/svalinn/DAGMC.svg?style=shield
    :target: https://circleci.com/gh/svalinn/DAGMC

Direct Accelerated Geometry Monte Carlo (DAGMC) is a software package that
allows users to perform Monte Carlo radiation transport directly on CAD models.

DAGMC has been integrated into a variety of Monte Carlo radiation codes
including MCNP5_, MCNP6_, Geant4_, FLUKA_, Tripoli4_, and Shift_. There are also
efforts planned to integrate DAGMC into other codes such as Serpent2_, OpenMC_,
Phits_, and FRENSIE_.

DAGMC currently relies on using the commercial solid modeling software Cubit_ (or its
`government-use counterpart <https://cubit.sandia.gov>`_ available from 
Sandia National Laboratories)
to prepare solid models. These packages can be
used to import CAD models from other tools such as SolidWorks, CATIA, etc., or
to create geometry from scratch. DAGMC also relies on Cubit to assign
materials and other geometry-related information.

For more information, please visit the `DAGMC website <DAGMC_>`_.

Quick links:

* `Install guide <https://svalinn.github.io/DAGMC/install/index.html>`_
* `Users guide <https://svalinn.github.io/DAGMC/usersguide/index.html>`_
* `Contributors guide <https://svalinn.github.io/DAGMC/contribute/index.html>`_

..  _DAGMC: https://svalinn.github.io/DAGMC
..  _Cubit: https://coreform.com/products/coreform-cubit/
..  _MCNP5: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/mcnp5.shtml
..  _MCNP6: https://mcnp.lanl.gov
..  _Geant4: https://geant4.cern.ch
..  _FLUKA: http://www.fluka.org/fluka.php
..  _Tripoli4: https://rsicc.ornl.gov/codes/ccc/ccc8/ccc-806.html
..  _Shift: https://meitner.ornl.gov/doe-codes/shift
..  _Serpent2: http://montecarlo.vtt.fi
..  _OpenMC: https://docs.openmc.org
..  _Phits: https://phits.jaea.go.jp
..  _FRENSIE: https://github.com/FRENSIE/FRENSIE
