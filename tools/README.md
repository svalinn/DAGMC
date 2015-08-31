About the dagmc_get_materials and test scripts
====
# General:
'dagmc_get_materials.py' is a Python script to handle the group names from a CAD model by creating a list of all group names and then extract the materials group names which then are used to get the specifications/metadata of each material (for exapmle the composition, density, atoms per molecule and the like). It then creates an output h5m file containing a directory that contains all the materials from the group names on the original CAD model along with all its metadata. In order to be able to copy the material metadata we need a materials library to copy from so some dependencies are required.
- The script can be run as a python script and all the flags needed can be found using:
  ```python dagmc_get_materials.py --help```
=======
Contributed Scripts for Workflow Improvements
================================================

This directory contains a number of scripts that have been contributed
to help facilitate the workflow with DAGMC.

Automated Imprint/Merge & Graveyard
-------------------------------------

The script `finish_dagmc_geom.bash` also relies on the Cubit journal
file `finish_dagmc_geom.jou`, to perform the following steps
automatically on an ACIS file.  For a given geometry file, say
`geom.sat`, this do the following steps:

* imprint & merge all volumes/surfaces
* export an ACIS file without graveyard for visualization: `geom_ng.sat`
* add a graveyard volume with an inner surface that is a brick 10%
  larger than the geometry bounding box and an outer surface 1% larger
  than the inner surface
* assign that volume to the `graveyard` group
* export an ACIS file with graveyard for transport: `geom_g.sat`
* convert the non-graveyard geometry to STL: `geom.stl`
* convert the graveyard geometry to H5M using default settings: `geom.h5m`

# Sample output:
-  An output ('sample_output.h5m') obtained by running the script with  test.h5m file. 
  
- test #2:
('test_output.py') tests the output h5m file created by running the script ans is run as a python script or using nosetests 

- test #3:
('test\_model.py') tests the group names on the CAD model with the output group names list from 'get\_tag\_values' function.
 
