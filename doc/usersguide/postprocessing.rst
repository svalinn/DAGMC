Post-processing
===============

Visualizing DAGMC geometry
~~~~~~~~~~~~~~~~~~~~~~~~~~

To use mbconvert one must already have a DAGMC h5m file ready to use, the
following command will convert the file to an stl file
::

    $ mbconvert <dagmc.h5m> <dagmc.stl>

It is often the case that DAGMC models contain so many triangles that it is
prohibitively slow to plot such a model even in Visit or Paraview, in such a
case it has proven useful to extract the faceted curves from the file and plot
those. This can be achieved with,
::

    $ mbconvert -1 dagmc.h5m curves.h5m
    $ mbconvert curves.h5m curves.vtk

Note that in the above example that first we must extract the curve information
write to a new h5m file, and then convert that file to a vtk file. When plotting
the curves that define the boundary of each volume should be visible, an example
of this is shown below.

..  image:: fng_curves.png
    :height: 300
    :width:  300
    :alt:    Image showing the FNG curve information
..  image:: fng_facets.png
    :height: 300
    :width:  300
    :alt:    Image showing the FNG facet information

Visualizing DAGMC mesh tally output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``mbconvert`` can be used to convert the output mesh file to a .vtk file for
viewing or post-processing with VisIt or other plotting tools.
::

    $ mbconvert mesh_out.h5m mesh_out.vtk
