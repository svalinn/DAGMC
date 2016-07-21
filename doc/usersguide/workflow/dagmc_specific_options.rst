Additional Input Parameters for DAGMC
++++++++++++++++++++++++++++++++++++++

DAGMC introduces a number of new input variables that control the
behavior of DAGMC when it interacts with your geometry.  This section
describes the conceptual behavior of those parameters and later
sections will describe their specific implementation and syntax within
each code.

Geometry File (basic)
"""""""""""""""""""""
* required
* Default: none

This file contains the geometry that has been created and
pre-processed using the workflow previously defined.  This file 
is a MOAB facet file (usually with a ``.h5m`` file extension).

Facet File 
""""""""""""""""""
* optional
* Default: none

For some models, the initial processing can be time consuming.  When
reading a solid model geometry, this option causes a processed file to
be written that can be used on subsequent analyses.  This file will be
processed with the facet tolerance as defined above.  This facet
tolerance cannot be changed when the file is reused.

Overlap Thickness (advanced)
"""""""""""""""""""""""""""""

* optional
* Default: 0.0

Often CAD geometry has small overlaps between adjacent volumes due to
file translation or imprecise modeling. The particle tracking
algorithm can accommodate small overlaps if the maximum thickness of
overlaps is approximately known.

Source Cell Treatment (intermediate)
""""""""""""""""""""""""""""""""""""

* optional
* Default: on (same behavior as native code)

The implementation of this option is specific to the Monte Carlo code
being used.  Please refer to the documentation for your underlying
Monte Carlo code.

Use Distance Limit (experimental)
"""""""""""""""""""""""""""""""""

* optional
* Default: off

This option allows a previously determined distance to the next
collision to be used to accelerate the search for ray-surface
intersections.  Any candidate ray-surface intersection that is at a
distance beyond the known distance to collision will be rejected,
including bounding box tests in the OBB tree search.

