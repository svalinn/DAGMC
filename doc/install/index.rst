DAGMC installation overview
===========================

This document describes the DAGMC installation process. Before you begin, you
should know that the following will be required:

1.  A basic understanding of Unix/Linux shell commands and how to navigate 
    the shell, for installing dependencies and DAGMC enabled codes
2.  MCNP5_ source code, if you intend to install DAG-MCNP5
3.  Fluka_, if you intend to install FluDAG
4.  Cubit_ or Trelis_, for the creation of geometry 

Once you have the basic pre-requisite peices in place you can proceed to the 
DAGMC installation. There are three main parts in the DAGMC installation:

1.  `Install the required dependencies <dependencies.html>`_ for DAGMC.
2.  `Install DAGMC <dagmc.html>`_  on top of the physics codes of your choice.
3.  `Install the CUBIT/Trelis plugin <plugin.html>`_  which will allow you to create and modify
    your own geometry.

Having completed the 3 stages of DAGMC installation, you will be ready to run
DAGMC based radiation transport calculations, the 
`user guide <../usersguide/index.html>`_ will guide you through the workflow options
available to you.

..  toctree::
    :hidden:
    :maxdepth: 1

    dependencies
    dagmc
    plugin

..  _MCNP5: https://mcnp.lanl.gov
..  _MCNP6: https://mcnp.lanl.gov
..  _Fluka: http://www.fluka.org/fluka.php
..  _Cubit: https://cubit.sandia.gov
..  _Trelis: http://www.csimsoft.com/trelis