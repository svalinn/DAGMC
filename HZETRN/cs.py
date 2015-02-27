#! /usr/bin/env python
import subprocess
import argparse
import string
import os
import ConfigParser
import shutil


try:
   from pyne import material
   from pyne.material import Material, MaterialLibrary
except:
    raise ImportError("The PyNE dependencies could not be imported.")    

try:
    from pyne import dagmc
except:
    raise ImportError("The dagmc dependency could not be imported.")    

from pyne import data
from pyne import nucname
import hzetrn as one_d_tool

class CS_CFG_Error(Exception):
    pass

# 'cross_dir' is a subdirectory for storing named input files
cross_dir   = 'cross'

def write_cs_input(header, mat_lib, path):
    """ write_cs_input
    Create cross_section input files in the current directory with the
    modified default name cs_input_FLUKA_NAME.dat for each material
    """
    cs_file_mats = {}
    num_materials = len(mat_lib.keys())
    for key in mat_lib:
        print 'in write_cs_input:', key, mat_lib[key]
        material_obj = mat_lib[key]
	# This can be removed after an update
	if 'Water' in key:
           h2o_atoms = {'H1': 2.0, 'O16': 1.0}
	   material_obj.from_atom_frac(h2o_atoms)
	# ToDo:  what materials do we not want to collapse?
	coll = material_obj.collapse_elements([])
	(fluka_name, xs_material_entry) = xs_create_entry(coll)
        xs_input_filename = 'cs_input_' + fluka_name + '.dat'
	xs_input_path = path + xs_input_filename
	# Create a dictionary of mat_name/input filename pairs.
	cs_file_mats[fluka_name] = xs_input_filename
        f = open(xs_input_path, 'w')
        f.write(header)
        f.write(xs_material_entry)
        f.close()
    return cs_file_mats
	

def xs_create_entry(coll):
    """ xs_create_entries_from_lib
    From the MaterialLibrary taken from the geometry file extract all the
    information needed for the cross-section input file.
    ToDo:  formatting of species line members
    """
    print coll
    fluka_name = coll.metadata['fluka_name']
    density = coll.density
    num_species = len(coll.comp)
    print "       ", fluka_name, ', ', coll.density, ', ', num_species
    material_entry = fluka_name + '\n' + "{0:.1f}".format(density) + '\n' + str(num_species) + '\n'
    for key in coll.comp:
        compname = nucname.name(key)
        str_comp_atomic_mass = "{0:.0f}".format(data.atomic_mass(key)) + '.'
        str_comp_charge = "{0:.0f}".format(nucname.znum(key)) + '.'
        str_comp_atoms_per_g = "{0:.2E}".format(coll.comp[key]*data.N_A/data.atomic_mass(key))
        print compname, str_comp_atomic_mass, str_comp_charge, str_comp_atoms_per_g
        material_entry += str_comp_atomic_mass + '  ' + str_comp_charge + '  ' + str_comp_atoms_per_g + '\n'
    return fluka_name, material_entry

def parsing():
    """
    Argument parsing
    returns : args: -f for the input geometry file, -r for the input ray tuple file
    ref     : DAGMC/tools/parse_materials/dagmc_get_materials.py
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', action='store', dest='uwuw_file', 
	help='The relative path to the .h5m file')

    parser.add_argument(
        '-d', action='store', dest='run_dir',
	help='The name of the holding directory for all hzetrn runs')
	
    args = parser.parse_args()

    if not args.uwuw_file:
        raise Exception('h5m file path not specified. [-f] not set')

    if not args.run_dir:
        raise Exception('Need a directory containing HZETRN executable structure.  [-d] not set.')

    return args

def main():
    """
    For each material in the material library created from the uwuw geometry 
    file create a separate input file to be used to calculate the required
    cross-sections for said material
    """
    global run_directory, cross_dir

    # Setup: parse the the command line parameters
    args = parsing()
    geom_filepath = os.path.join(os.path.dirname('__file__'), args.uwuw_file)

    cur_path = os.path.dirname(os.path.abspath(__file__)) + '/'
    
    # Create some generic full paths 
    run_path = cur_path + args.run_dir + '/'
    cross_path = cur_path + cross_dir + '/'

    # Make the cross_path subdir if necessary
    if not os.path.isdir(cross_path):
        print 'Creating path', cross_path
        os.mkdir(cross_path)
    for pathname in (run_path, cross_path):
        assert os.path.isdir(pathname), "Path '{0}' not found.".format(pathname)

    header = one_d_tool.xs_create_header()

    # Cross-Section: load the material library from the uwuw geometry file
    mat_lib = material.MaterialLibrary()
    print geom_filepath
    mat_lib.from_hdf5(geom_filepath)
    # Create a cs_input file for each material in the problem
    cs_input_for_name = write_cs_input(header, mat_lib, cross_path) 

    for fname, cs_input in cs_input_for_name.iteritems():
        one_d_tool.cross_section_process(cross_path + cs_input, run_path, fname)

    return
    ##########################################################
   
if __name__ == '__main__':
    main()

