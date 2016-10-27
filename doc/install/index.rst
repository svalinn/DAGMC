DAGMC installation overview
===========================

This document describes the DAGMC installation process. Before you begin, you
should know that the following will be required:

1. A basic understanding of shell commands and how to navigate the shell, for
   installing dependencies and DAGMC
2. MCNP5 source code, if you intend to install DAG-MCNP5
3. FLUKA, if you intend to install FluDAG
4. Cubit or Trelis, for the creation of geometry
5. A glass of wine, so you don't get as frustrated when something doesn't work
   (optional)

There are three main parts in the DAGMC installation.

1. Install the required dependencies for DAGMC.
2. Install DAGMC on top of the physics codes of your choice.
3. Install the CUBIT/Trelis plugin which will allow you to create and modify
   your own geometry.

.. toctree::
   :maxdepth: 1

   dependencies
   dagmc
   plugin
