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

class DagmcError(Exception):
    pass

""" 
xs_create_header creates a string conting the first few material-
independent lines of the input file for the cross-section call

returns string contining the lines
Note: this is not strictly a header, as what1 and what 2 have to match
the actual directories.
"""
def xs_create_header(common_data_dir, cs_out_dir):

    comment1 = '# Name of folder where static data is stored'
    comment2 = '# Name of folder to put data in (folder must already exist)'
    # 84 dashes
    divider  = '{:-<84}'.format('-')
    # right-align, pad with spaces
    lines  = '{:<36}'.format(common_data_dir) + comment1 + '\n'
    lines += '{:<36}'.format(cs_out_dir) + comment2 + '\n'
    lines += '\n'
    lines += divider + '\n'
    lines += '\n'
    return lines

""" one_cs_for_each
Create cross_section input files in the current directory with the
modified default name cs_input_MATERIAL.dat for each material

ToDo:  Think of a better name for this method
"""
def one_cs_for_each(header, mat_lib):
    cs_file_mats = {}
    num_materials = len(mat_lib.keys())
    for key in mat_lib:
        material_obj = mat_lib[key]
	# This can be removed after an update
	if 'Water' in key:
           h2o_atoms = {'H1': 2.0, 'O16': 1.0}
	   material_obj.from_atom_frac(h2o_atoms)
	# ToDo:  what materials do we not want to collapse?
	coll = material_obj.collapse_elements([])
	name, xs_material_entry = xs_create_entry(coll)
        xs_input_filename = 'cs_input_' + name + '.dat'
	# Create a dictionarey of mat_name/input filename pairs.
	# Slightly redundant, but separate access will be useful.
	cs_file_mats[name] = xs_input_filename
        f = open(xs_input_filename, 'w')
        f.write(header)
        f.write(xs_material_entry)
        f.close()
    return cs_file_mats
	

""" xs_create_entries_from_lib
From the MaterialLibrary taken from the geometry file extract all the
information needed for the cross-section input file.
ToDo:  formatting of species line members
"""
def xs_create_entry(coll):
    name1 = coll.metadata['fluka_name']
    density = coll.density
    num_species = len(coll.comp)
    print "       ", name1, ', ', coll.density, ', ', num_species
    material_entry = name1 + '\n' + "{0:.1f}".format(density) + '\n' + str(num_species) + '\n'
    for key in coll.comp:
        compname = nucname.name(key)
        str_comp_atomic_mass = "{0:.0f}".format(data.atomic_mass(key)) + '.'
        str_comp_charge = "{0:.0f}".format(nucname.znum(key)) + '.'
        str_comp_atoms_per_g = "{0:.2E}".format(coll.comp[key]*data.N_A/data.atomic_mass(key))
        print compname, str_comp_atomic_mass, str_comp_charge, str_comp_atoms_per_g
        material_entry += str_comp_atomic_mass + '  ' + str_comp_charge + '  ' + str_comp_atoms_per_g + '\n'
    return name1, material_entry

"""
Argument parsing
returns : args: -f for the input geometry file, -r for the input ray tuple file
ref     : DAGMC/tools/parse_materials/dagmc_get_materials.py
"""
def parsing():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', action='store', dest='uwuw_file', 
	help='The relative path to the .h5m file')

    parser.add_argument(
        '-c', action='store', dest='config_file',
	help='The name of the cs config file')

    args = parser.parse_args()

    if not args.uwuw_file:
        raise Exception('h5m file path not specified. [-f] not set')
    # Should cs.cfg be in the protected directory?
    if not args.config_file:
        args.config_file = 'cs.cfg'

    return args


def main():
    # Setup: parse the the command line parameters
    args = parsing()
    path = os.path.join(os.path.dirname('__file__'), args.uwuw_file)
    config = ConfigParser.ConfigParser()
    config.read(args.config_file)


    # Start the file with header lines which contain the names of 
    # some folders the cross_section processing will need
    cs_common = config.get("cs_files", "common_data")
    cs_outdir = config.get("cs_files", "cs_out")
    xs_header = xs_create_header(cs_common, cs_outdir)
    print xs_header

    # Cross-Section: load the material library from the uwuw geometry file
    mat_lib = material.MaterialLibrary()
    mat_lib.from_hdf5(path)
    cs_file_dict = one_cs_for_each(xs_header, mat_lib) 

    # Using the input file just created, prepare the materials subdirectory
    # This method will make a subprocess call
    curdir = os.path.dirname(os.path.abspath(__file__))
    subdir = curdir + '/' + config.get("cs_files", "run_directory")
    for name in cs_file_dict:
	src = curdir + '/' + cs_file_dict[name]
        one_d_tool.cross_section_process(subdir, src, name, cs_outdir)
    ##########################################################
    return
   
if __name__ == '__main__':
    main()

