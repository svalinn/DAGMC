.. DAGMC documentation master file, created by
   sphinx-quickstart. You can adapt this file
   completely to your liking, but it should at least
   contain the root `toctree` directive.

DAGMC:  The Direct Accelerated Geometry Monte Carlo Toolkit
===========================================================

DAGMC is an open-source toolkit for performing Monte Carlo radiation transport
simulations on CAD-based geometry models. It is built on top of the MOAB mesh
database. DAGMC is developed and maintained by the `Computational Nuclear
Engineering Research Group (CNERG) <https://cnerg.github.io/>`_ at the University of Wisconsin -- Madison.

**Please see the** :ref:`dagmc-intro` **section for more information on DAGMC and its role in Monte Carlo particle transport.**

.. carousel::
    :show_controls:
    :show_fade:
    :show_captions_below:

    .. figure:: assets/manifold-cad.png
       :height: 400px

       CAD

       A pipe manifold modeled in :term:`Coreform Cubit`.

    .. figure:: assets/manifold-tris.png
       :height: 400px

       Surface Mesh w/ Embedded Topology

       A DAGMC surface mesh of the piping manifold.

    .. figure:: assets/manifold-flux.png
       :height: 400px

       Flux Results

       A flux mapping generated using :term:`OpenMC`.

Table of Contents
=================

.. toctree::
    :maxdepth: 2

    intro/index
    usersguide/index
    contribute/index
    theoryguide/index
    filespec/index
    api/index
    gallery/index
    glossary