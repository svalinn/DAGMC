DAGMC unstructured mesh tallies
===============================

DAGMC supports several advanced mesh tally options including tetrahedral mesh
tallies (tetmeshes) and kernel density estimator tallies. Both options require
an unstructured mesh.

Mesh production workflow
~~~~~~~~~~~~~~~~~~~~~~~~

Cubit/Trelis can be used to generate the unstructured meshes needed for tallies.
To do so, use the following steps.

1.  Load the geometry you wish to mesh into Cubit/Trelis.
2.  Use the mesh tools to produce the meshes you want.
3.  Save the file as a .trelis or .cub file. Remember to check the "Use Legacy .cub file Format" option in Trelis or Cubit.
4.  Use MOAB's ``mbconvert`` executable to convert from the Cubit/Trelis format
    to a faceted .h5m file that DAGMC can use.

Here is an example of how to use ``mbconvert``:
::

    $ mbconvert mesh.cub mesh.h5m

One can immediately view this mesh, by using ``mbconvert`` to convert to vtk:
::

    $ mbconvert mesh.h5m mesh.vtk

Tetrahedral mesh tallies
~~~~~~~~~~~~~~~~~~~~~~~~

Currently the tetmesh tally option is only available in DAG-MCNP5, but there are
plans to expand this to all the codes that DAGMC currently supports.

Once you have used the above instructions to produce one or more meshes in .h5m
format, you must specify to MCNP that your tally is a tetmesh tally. This is
done by putting special keywords on the tally comment (FC) line. For example,
adding the following lines to an MCNP input file would tell MCNP that tally 4 is
a tetmesh tally, that the input mesh file is ``mesh.h5m``, and the results
should be stored in ``mesh_out.h5m``.
::

    fmesh4:n geom=dag type=unstr_track
    fc4 dagmc inp=mesh.h5m out=mesh_out.h5m

Other standard MCNP options can also be used, such as energy bins:
::

    fmesh4:n geom=dag type=unstr_track
             emesh=1.0 2.0 15.0
    fc4 dagmc inp=mesh.h5m out=mesh_out.h5m

Or tally multipliers:
::

    fmesh4:p geom=dag type=unstr_track
    fc4 dagmc inp=mesh.h5m out=mesh_out.h5m
    fm4 -1 0 -5 -6

``mbconvert`` can be used to convert the output mesh file to a .vtk file for
viewing or post-processing with VisIt_ or other plotting tools.
::

    $ mbconvert mesh_out.h5m mesh_out.vtk

Kernel density estimator tallies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The kernel density estimator tallies are a suite of tallies that allow scoring
to be done using the vertices of a tetmesh or add additional mesh-based scoring
such as collision estimator tallies. The theoretical details behind KDE tallies
can be found in `Kerry Dunn's Ph.D. thesis
<http://digital.library.wisc.edu/1711.dl/OXDMBPODZJERF8A>`_.

To call a KDE collision tally, use:
::

    fmesh4:p geom=dag type=kde_coll
    fc4 dagmc inp=mesh.h5m out=mesh_out.h5m
        hx=0.198 hy=0.0663 hz=0.0662

To call a KDE track length tally, use:
::

    fmesh4:p geom=dag type=kde_track
    fc4 dagmc inp=mesh.h5m out=mesh_out.h5m
        hx=0.198 hy=0.0663 hz=0.0662

To call a KDE subtrack tally, use:
::

    fmesh4:p geom=dag type=kde_subtrack
    fc4 dagmc inp=mesh.h5m out=mesh_out.h5m
        hx=0.1042 hy=0.0833 hz=0.0833
        hx=0.1042 hy=0.0833 hz=0.0833
        subtracks=3 seed=11699913

.. _VisIt: https://wci.llnl.gov/simulation/computer-codes/visit
