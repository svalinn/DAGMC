DAGMC tools
===========

The DAGMC infrastructure has a number of tools that are used to process, fix, diagnose and repair
processed models. This document serves to better inform and
educate what those tools are and how to use them.

make_watertight
~~~~~~~~~~~~~~~

The ``make_watertight`` tool has `its own page <watertightness>`_.

mbconvert
~~~~~~~~~

The mbconvert tool from the MOAB tool is used to translate MOAB meshes into some text based format. It is useful for converting
DAGMC models into some other visualization form for result post processing, such as stl files for viewing in Visit or Paraview. To
run mbconvert input option, input filename and output filename are specificed:
::

    $ mbconvert <input_filename> <output_filename> [options]

The options which control mbconvert are:
::

    -v  <int> - extract a specific volume or range of volumes
    -s  <int> - extract a specific surface or range of surfaces
    -c  <int> - extract a specific curve or range of curves
    -1  - extract edges only
    -2  - extract two dimensional entites only e.g. Tri, Quad, etc.
    -3  - extract three dimensional entities only, e.g. tet, hex, etc.
    -h  - print help

To use mbconvert one must already have a DAGMC h5m file ready to use, the following command will convert the file to an stl file
::

    $ mbconvert <dagmc.h5m> <dagmc.stl>

It is often the case that DAGMC models contain so many triangles that it is prohibitively slow to plot such a model even in
Visit or Paraview, in such a case it has proven useful to extract the faceted curves from the file and plot those. This can be
achieved with,
::

    $ mbconvert -1 dagmc.h5m curves.h5m
    $ mbconvert curves.h5m curves.vtk

Note that in the above example that first we must extract the curve information write to a new h5m file, and then convert that
file to a vtk file. When plotting the curves that define the boundary of each volume should be visible, an example of this is
shown below.

.. image:: fng_curves.png
   :height: 300
   :width:  300
   :alt:    Image showing the FNG curve information
.. image:: fng_facets.png
   :height: 300
   :width:  300
   :alt:    Image showing the FNG facet information

mklostvis
~~~~~~~~~

Sometimes either poor quality CAD, incorrect imprinting & merging or overlapping volumes; particles are regarded as lost
by the Monte Carlo code. It is therefore neccessary to be able to examine where the particles were lost and in which direction they
were travelling in. The tool `mklostvis <https://github.com/svalinn/meshtools/tree/master/lostparticles>`_ is designed for this
purpose, reading the output of the MCNP lost particle information and producing a `Cubit <https://cubit.sandia.gov/>`_ or
`Trelis <http://www.csimsoft.com/trelis.jsp>`_ journal file which will draw the lost particles as vertices and their directions as curves.
To run the script;
::

    $ mklostvis.pl [mcnp output filename] [vector length] > [journal file name]

The produced Cubit journal file can be 'played', and will plot these lost particles on top of whatever geometry is loaded into
your Cubit session, like that shown below.

.. image:: lost_p.png
   :height: 300
   :width:  300
   :alt:    Image showing lost particle information
.. image:: lost_p_zoom.png
   :height: 300
   :width:  300
   :alt:    Image showing lost particle information zoomed in
