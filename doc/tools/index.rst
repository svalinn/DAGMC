DAGMC Tools
======

The DAGMC infrastructure has a number of tools that can help with issues that
arise when using the toolkit. This document serves to better inform and 
educate what those tools are and how to use them.

mcnp2cad
~~~~~~~~~~~~
UW has a tool, `mcnp2cad <https://github.com/svalinn/mcnp2cad>`_, which can translate a MCNP model
into ACIS format for use in future DAGMC simulations. The tool builds CAD solids from MCNP cell descriptions
often turning infinite bodies and plans into finite versions. At the time of writing the number of unsupported
MCNP surface descriptions is limited to GQ's and SQ's. To run mcnp2cad all that is needed is an MCNP input deck,
::
   prompt%> mcnp2cad test.inp

Will result in a file called out.sat which will contain the CAD version of your MCNP input. Furthermore, mcnp2cad
automatically transfers material and importance assignments into the CAD model and will be translated to the 
DAGMC file when processed.

dagmc_prepoc
~~~~~~~~~~~~
Previously it was typical of DAGMC users to input the acis file directly into DAGMC, the file would
be faceted and then immediately loaded into memory however, as models have grown in complexity, the time
required to facet and build the OBB tree has grown immensely and it is now recommended that faceting is
treated as a preprocess step.

When MOAB is build with CGM support (see `build instructions <get_install.html>`_) the tool 
dagmc_preproc is automatically built. You can confirm that you have MOAB in your path and that
the executable is built, by running
::
   prompt%> which dagmc_preproc
To which you should get a reply like
::
   prompt%> /data/opt/$USER/moab/bin/dagmc_preproc

If you get a reply like above then you are able to use dagmc_preproc. If not, you need to recompile with
moab support following the build instructions. Once the model is in the correct form as required
by our `workflow <workflow.html>`_ dagmc_preproc may be used. The dagmc_prepoc tool is invoked by,
::
   prompt%> dagmc_preproc <filename.sat> -o <output_filename.h5m> [options]
There are several options that are commonly used by dagmc_preproc to control the state of the output file.
::
   -f <float> - the facet distance tolerance to be used, 
                i.e. the maximum distance (cm) the facet 
                may lie from the curve 
   -l <float> - the facet length tolerance to be used, 
                i.e. the maximum length of the side of 
                the facet (cm)
   -a <int> - the angle tolerance to be used, i.e. 
               the maximum angle that the facet can 
                deviate from the surface normal
   -o <char> - the name of the output facet file to create
   --fatal_curves - fatal error when finding curves that cannot be faceted 
   -h - print help
For example, to facet a file called geom.sat into a file called facet_file.h5m with a facet tolerance of 1e-4 cm:
::
   prompt%> dagmc_preproc geom.sat -f 1.0e-4 -o facet_file.h5m

This file can then be used by DAGMC, MOAB and other down stream tools.

mbconvert
~~~~~~~~~
The mbconvert tool from the MOAB tool is used to translate MOAB meshes into some text based format. It is useful for converting
DAGMC models into some other visualization form for result post processing, such as stl files for viewing in Visit or Paraview. To
run mbconvert input option, input filename and output filename are specificed: 
::
   prompt%> mbconvert <input_filename> <output_filename> [options]

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
   prompt%> mbconvert <dagmc.h5m> <dagmc.stl>

It is often the case that DAGMC models to contain so many triangles that it is prohibitively slow to plot such a model even in 
Visit or Paraview, in such a case it has proven useful to extract the faceted curves from the file and plot those. This can be 
achieved with, 
::
   prompt%> mbconvert -1 dagmc.h5m curves.h5m
   prompt%> mbconvert curves.h5m curves.vtk

Note that in the above example that first we must extract the curve information write to a new h5m file, and then convert that
file to a vtk file. When plotting the curves that define the boundary of each volume should be visible, an example of this is 
shown below.

make_watertight
~~~~~~~~~~~~~~~
When models are faceted by dagmc_preproc, facets are not guarenteed to be water tight, by which it is meant that 
the edge of facets of one surface do not exactly align with those of another surface, such discrepancies can lead to numerical
gaps through which particles can leak. The make_watertight tool alleviates this problem by using the underlying faceted curve 
representation to force triangle edges to be coincident along curves, this action closes any gaps and if succesful the model is 
guarenteed to be water tight. The tool is run in the following way,
::
   prompt%> make_wateright <filename>

The result of this step is a new file with the name, <filename>_zip.h5m, which can then be run with the rest of the workflow. The
degree of water tightness can be checked with the partner tool to make_wateright, check_wateright. When run it will give a summary
of how sealed the mode is. The check_watertight tool is run by:
::
   prompt%> check_watertight <filename>


   
