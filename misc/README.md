About uwuw_preproc and test scripts
====
General:
'uwuw_preproc' is a Python script that handles group names from a CAD model
by creating a list of all group names included then extracts materials group
names which are then used to get the specifications/metadata of each material 
(for exapmle the composition, density, atoms per molecule, etc.). 
It then creates an output h5m (or writes to the geometry h5m) file a material 
library object that contains all materials from the original CAD model along 
with all its metadata. In order to be able to copy the material metadata we 
need a materials library to copy from so some dependencies are required.

The material Library provided is expected to have a specified structure for the 
script to effectively copy all materials metadat from. This can be achieved by 
using PyNE functions to create the material library.
e.g. write_hdf5("filename",nucpath="/material_library/nucid",datapath="/material_library/materials") 
All needed info can be found here: 
http://pyne.io/pyapi/material.html#pyne.material.Material

The script can be run as a python script and all the flags needed can be found using:
  ```python uwuw_preproc -h```
=======
Contributed Scripts for Workflow Improvements
================================================

This directory contains a number of scripts that have been contributed
to help facilitate the workflow with DAGMC.

Automated Imprint/Merge & Graveyard
-------------------------------------

The script `finish_dagmc_geom.bash` also relies on the Cubit journal file `finish_dagmc_geom.jou`,
to perform the following steps automatically on an ACIS file.  For a given geometry file, say `geom.sat`, 
this do the following steps:

* imprint & merge all volumes/surfaces
* export an ACIS file without graveyard for visualization: `geom_ng.sat`
* add a graveyard volume with an inner surface that is a brick 10% larger than the geometry 
  bounding box and an outer surface 1% larger than the inner surface
* assign that volume to the `graveyard` group
* export an ACIS file with graveyard for transport: `geom_g.sat`
* convert the non-graveyard geometry to STL: `geom.stl`
* convert the graveyard geometry to H5M using default settings: `geom.h5m`

# tests:
In tools/tests directory a sample model was provided for testing purposes. 
- test #1:
('test_uwuw_preproc_functiond.py') tests all functions in the uwuw_preproc script and is run as a python script or using nosetests 
- test #2:
('test_uwuw_preproc_model.py') tests the group names on the CAD model with the output group names list from 'get_tag_values' function.
 
