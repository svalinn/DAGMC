DAGMC Tallies
----------------

DAGMC Mesh Tallies
----------------
There is support in DAGMC for a number of advanced mesh tally options, namely a Tetrahedral Mesh Tally and a Kernel Density
Esimator tally. Both mesh tally options require an unstructured mesh 

Mesh Production Workflow
===============
The meshes needed for tallying upon are generated using Cubit/Trelis, load the geometry that you wish to mesh, and using
the meshing tools produces meshes on the components on which you would like to determine neutronic results. Once meshing
is completed, save the file as a *.cub or *.trelis file. We use MOAB to convert from the Cubit/Trelis format to a
h5m that DAGMC can use.
::
   mbconvert file.cub file.h5m

This mesh is now ready for use in scoring in DAGMC.

Tetrahedral Mesh
===============
Currently the Tetrahedral Mesh (TM) tally option is only available in DAG-MCNP5, but there are plans to expand this to all
the codes that DAGMC currently supports. Having produced a mesh in the previous step, the following lines must be used to
allow the tally to be used in MCNP5:
::
    fmesh4:n geom=dag
    fc4 dagmc inp=mesh.h5m out=mesh_out.h5m

This will produce a TM tally with the mesh from mesh.h5m and tally results will be written to mesh_out.h5m. Other standard
MCNP options can be used, for example including energy bins:
::
    fmesh4:n geom=dag 
             emesh=1.0 2.0 15.0
    fc4 dagmc inp=mesh.h5m out=mesh_out.h5m
Or tally mulipliers:
::
    fmesh4:p geom=dag type=unstr_track
    fc4 dagmc inp=mesh.h5m out=mesh_out.h5m
    fm4 -1 0 -5 -6
The resultant mesh_out.h5m can be converted to a *.vtk file for post processing with Visit or other plotting tools.

Kernel Density Estimator
===============
The Kernel Density Estimator tallies are a suite of tallies that allow scoring to be done
using the vertices of a TM or add additional mesh based scoring such as collision estimator
tallies. The theoretical details behind KDE tallies can be found
in `Kerry Dunn's PhD thesis <http://digital.library.wisc.edu/1711.dl/OXDMBPODZJERF8A>`_. 

To call these tally methods from MCNP the following fmesh calls are used. To call a KDE
Collision tally, the following call is used.
::
    fmesh4:p geom=dag type=kde_coll
    fc4 dagmc inp=mesh.h5m out=mesh_out.h5m
       hx=0.198 hy=0.0663 hz=0.0662
To call a KDE track length tally, the following call is used.
::
    fmesh4:p geom=dag type=kde_track
    fc4 dagmc inp=mesh.h5m out=mesh_out.h5m
       hx=0.198 hy=0.0663 hz=0.0662 
To call a KDE subtrack tally, the following call is used
::
    fmesh4:p geom=dag type=kde_subtrack
    fc4 dagmc inp=mesh.h5m out=mesh_out.h5m
        hx=0.1042 hy=0.0833 hz=0.0833
        hx=0.1042 hy=0.0833 hz=0.0833 subtracks=3 seed=11699913 
