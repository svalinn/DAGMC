#!/usr/bin/python
import subprocess
import argparse
import string
try:
    from itaps import iMesh, iBase
except:
    raise ImportError("The PyTAPS dependencies could not be imported.")    
try:
   from pyne import material
   from pyne.material import Material, MaterialLibrary
except:
    raise ImportError("The PyNE dependencies could not be imported.")    

import hzetrn as one_d_tool


"""
Load the pyne material library
filename : string of the name of the material library
returns  : PyNE MaterialLibrary instance
ref      : DAGMC/tools/parse_materials/dagmc_get_materials.py
"""

def load_ray_tuples(filename):
    # do I need to declare this outside the if?
    ray_tuples = []

    if not filename:
       ray_tuples = [(1.0, 0.0, 0.0),
                     (0.0, 1.0, 0.0),
	   	     (0.0, 0.0, 1.0)]
    else:
    # placeholder for code
    # In future, this will be a randomly generated list
       print('The ray tuples file is')
       print(filename)
       ray_tuples = [(1.0, 0.0, 0.0),
                     (0.0, 1.0, 0.0),
	   	     (0.0, 0.0, 1.0)]
    return ray_tuples

def 
"""
Dummy implementation of argument parsing
returns : args: -f for the input geometry file, -r for the input ray tuple file
ref     : DAGMC/tools/parse_materials/dagmc_get_materials.py
"""
def parsing():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', action='store', dest='datafile', 
	help='The path to the .h5m file')
    parser.add_argument(
        '-r', action='store', dest='ray_dir_file', 
	help='The path to the file with ray direction tuples')

    args = parser.parse_args

    if not args.datafile:
        raise Exception('h5m file path not specified. [-f] not set')
    if not args.nuc_data:
        raise Exception('nuc_data file path not specified. [-d] not set')
    if not args.ray_dir_file:
        args.ray_dir_file = FALSE

    # ToDo: look at what this is used for
    if not args.output:
       args.output = 'dummy'
    return args

def main():
    # parse the the command line parameters
    args = parsing()
    
    # get list of rays
    ray_tuples = load_ray_tuples(args.ray_dir_file) 

    # Prepare xs_input.dat
    input_filename = "xs_input.dat"
    # ToDo: Start the file with header lines which contain the names of 
    #       some folders the cross_section processing will need
    #       These may be hard-coded to start with; they can always be
    #       edited later
    # create_and_insert_header_lines(input_filename)

    # load the material library from the geometry file
    mat_lib = material.MaterialLibrary()
    mat_lib.from_hdf5(args.datafile)

    # For each material create an entry and append to input_filename
    for key in mat_lib.iterkeys():
    	material_obj = mat_lib.get(key)
	material_entry = get_xs_info(material_obj)
	# ToDo: append material_entry to input_filename
	# append_entry(input_filename, material_entry)

    # This method will make a subprocess call
    1_d_tool.cross_section_process(input_filename)


    
    #
    # load geometry in pydagmc
    dag_geom=iMesh.Mesh()
    dag_geom.load(args.datafile)
    # PyTAPS!
    dag_geom.getEntities()
    # from get_tag_values in dagmc_get_materials.py
    # maybe I need these params?
    # dag_geom.getEntities(iBase.Type.all, iMesh.Topology.triangle)


