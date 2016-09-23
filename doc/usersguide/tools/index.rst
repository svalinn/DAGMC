DAGMC Tools
===========

The DAGMC infrastructure has a number of tools that are used to process, fix, diagnose and repair
processed models. This document serves to better inform and
educate what those tools are and how to use them.

mcnp2cad
~~~~~~~~

UW has a tool, `mcnp2cad <https://github.com/svalinn/mcnp2cad>`_, which can translate a MCNP model
into ACIS format for use in future DAGMC simulations. The tool builds CAD solids from MCNP cell descriptions
often turning infinite bodies and planes into finite versions. At the time of writing, the only unsupported
MCNP surface descriptions are limited to GQ's and SQ's. To run mcnp2cad all that is needed is an MCNP input deck,
::

    $ mcnp2cad test.inp

Will result in a file called out.sat which will contain the CAD version of your MCNP input. Furthermore, mcnp2cad
automatically transfers material and importance assignments into the CAD model and will be translated to the
DAGMC file when processed.

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

make_watertight
~~~~~~~~~~~~~~~

When models are faceted by dagmc_preproc, facets are not guarenteed to be water tight, by which it is meant that
the edge of facets of one surface do not exactly align with those of another surface, such discrepancies can lead to numerical
gaps through which particles can leak. The make_watertight tool alleviates this problem by using the underlying faceted curve
representation to force triangle edges to be coincident along curves, this action closes any gaps and if succesful the model is
guarenteed to be water tight. The tool is run in the following way,
::

    $ make_wateright <filename>

The result of this step is a new file with the name, <filename>_zip.h5m, which can then be run with the rest of the workflow. The
degree of water tightness can be checked with the partner tool to make_wateright, check_wateright. When run it will give a summary
of how sealed the mode is. The check_watertight tool is run by:
::

    $ check_watertight <filename>

The make_watertight tool is built as part of the DAGMC build process.

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
