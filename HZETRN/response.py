import subprocess
import argparse
import string
import os

# Mesh specific imports                                                                                                             
try:
    from itaps import iMesh, iBase
except:
    raise ImportError("The iMesh dependencies could not be imported.")

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
import numpy as np
import hzetrn as one_d_tool
import tag_utils

class DagmcError(Exception):
    pass

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

def get_rand_dir():
    # [0.0,1.0)
    rnum = np.random.uniform()
    # map rnum onto [0.0,1.0)
    z = 2*rnum - 1
    theta = 2*np.pi*rnum
    norm_fac = np.sqrt(1 - z*z)
    y = norm_fac*np.sin(theta)
    x = norm_fac*np.cos(theta)
    return [x, y, z]

"""
Dummy implementation of argument parsing
returns : args: -f for the input geometry file, -r for the input ray tuple file
ref     : DAGMC/tools/parse_materials/dagmc_get_materials.py
"""
def parsing():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', action='store', dest='uwuw_file', 
	help='The relative path to the .h5m file')

    parser.add_argument(
        '-r', action='store', dest='ray_dir_file', 
	help='The path to the file with ray direction tuples')

    args = parser.parse_args()

    if not args.uwuw_file:
        raise Exception('h5m file path not specified. [-f] not set')
    if not args.ray_dir_file:
        args.ray_dir_file = False

    return args


""" Find the volume of the geometry that contains the ref_point
"""
def find_ref_vol(ref_point):
    # code snippet from dagmc.pyx:find_graveyard_inner_box
    volumes = dagmc.get_volume_list() 
    ref_vol = 0 
    for v in volumes:
        if dagmc.point_in_volume(v, ref_point):
            ref_vol = v

    return ref_vol

""" 
xs_create_header creates a string conting the first few material-
independent lines of the input file for the cross-section call

returns string contining the lines
Note: this is not strictly a header, as what1 and what 2 have to match
the actual directories.
"""
def xs_create_header():
    what1 = 'common_static_data'
    what2 = 'cross_sections_out'
    comment1 = '# Name of folder where static data is stored'
    comment2 = '# Name of folder to put data in (folder must already exist)'
    # 84 dashes
    divider  = '{:-<84}'.format('-')
    # right-align, pad with spaces
    lines  = '{:<36}'.format(what1) + comment1 + '\n'
    lines += '{:<36}'.format(what2) + comment2 + '\n'
    lines += '\n'
    lines += divider + '\n'
    lines += '\n'
    return lines

""" create_entries_from_lib
From the MaterialLibrary taken from the geometry file extract all the
information needed for the cross-section input file.
ToDo:  formatting of species line members
"""
def create_entries_from_lib(mat_lib):
    # For each material create an entry and append to input_filename
    material_entries = []
    num_materials = len(mat_lib.keys())
    for key in mat_lib.iterkeys():
   	material_obj = mat_lib.get(key)
	# ToDo:  need exclusion list
	coll = material_obj.collapse_elements([])
	name1 = coll.metadata['name']
	density = coll.density
	num_species = len(coll.comp)
	# Don't know if I need this
	number_density = coll.number_density()
	material_entry = name1 + '\n' + str(density) + '\n' + str(num_species) + '\n'
	for key in coll.comp:
	    compname = nucname.name(key)
	    str_comp_atomic_mass = str(data.atomic_mass(key))
	    str_comp_charge = str(nucname.znum(key))
	    str_comp_atoms_per_g = str(coll.comp[key]*data.N_A/data.atomic_mass(key))
	    # print compname, str_comp_atomic_mass, str_comp_charge, str_comp_atoms_per_g
	    material_entry += str_comp_atomic_mass + '  ' + str_comp_charge + '  ' + str_comp_atoms_per_g
        material_entries.append(material_entry)
    return material_entries

"""
For a given volume get the name and the density of the meterial in it
"""
def get_name_density(vol, mat_assigns, mat_lib):
    # Get the material in each vol
    name = ''
    density = -1.0
    foundMat = False
    
    for pair in mat_assigns:
        if foundMat:
            break
	if vol == pair[0]:
	    if 'graveyard' in pair[1]:
	        name = 'graveyard'
	        foundMat = True
	    else:
                for mat in mat_lib:
		    if foundMat:
		        break
		    if pair[1] == mat:
	    	        name = pair[1]
		        density = mat_lib[mat].density
			foundMat = True
                    
    if not foundMat:
        name = 'implicit_complement'
    return name, density

def main():
    # Setup: parse the the command line parameters
    args = parsing()
    path = os.path.join(os.path.dirname('__file__'), args.uwuw_file)
    # ToDo: Start the file with header lines which contain the names of 
    #       some folders the cross_section processing will need
    #       These may be hard-coded to start with; they can always be
    #       edited later
    xs_header = xs_create_header()

    # Cross-Section: load the material library from the uwuw geometry file
    mat_lib = material.MaterialLibrary()
    mat_lib.from_hdf5(path)
    xs_material_entries = create_entries_from_lib(mat_lib)

    # Prepare xs_input.dat
    xs_input_filename = 'xs_input.dat'
    f = open(xs_input_filename, 'w')
    f.write(xs_header)
    f.write("\n".join(xs_material_entries))
    f.close()

    # Using the input file just created, prepare the materials subdirectory
    # This method will make a subprocess call
    one_d_tool.cross_section_process(xs_input_filename)
    # End material cross-section part
    ##########################################################
    # Get the DAG object for this geometry
    rtn = dagmc.load(path)
    # dagmc.load gets us only so far, so get list of tag values
    (tag_values,mat_assigns) = tag_utils.get_mat_tag_values(path)
    print ('#####################################################################')
    print ('mat_assigns: ', mat_assigns)
    """ Produces:
    ('mat_assigns: ', [(1, 'mat:Lead'), (2, 'mat:Lead'), (3, 'mat:Lead/rho:12.3'), (4, 'mat:Lead/rho:12.3'), (7, 'mat:graveyard')])
    #####################################################################
    """ 
    ##########################################################
    # get list of rays
    ray_tuples = load_ray_tuples(args.ray_dir_file) 

    # Use 0,0,0 as a reference point for now
    ref_point = [0.0, 0.0, 0.0]
    xyz = ref_point
    start_vol = find_ref_vol(ref_point)
    surf = 0
    dist = 0
    huge = 1000000000
    slab_length = []
    slab_mat_name = []
    slab_density = []
    vols_traversed = []

    v = start_vol
    for dir in ray_tuples:
	print ('#################### dir = ', dir, '###################################') 
	last_vol = start_vol
        for (vol, dist, surf, xyz) in dagmc.ray_iterator(start_vol, ref_point, dir, yield_xyz=True):
            # print('  next intersection at distance', dist, 'on surface', surf)
	    # print('  new xyz =', xyz, 'proceeding into volume', vol)
	    if (dist < huge) and (surf != 0):
	        slab_length.append(dist)
		# need data from vol we are leaving
		name, density = get_name_density(last_vol, mat_assigns, mat_lib)
		slab_mat_name.append(name)
		slab_density.append(density)
		vols_traversed.append(last_vol)
	        last_vol = vol

	print ('lengths ',  slab_length)
	print ('materials', slab_mat_name)
	print ('densities', slab_density)
	print ('vols traversed', vols_traversed)
      
        slab_length = []
        slab_mat_name = []
        slab_density = []
	vols_traversed = []

    ###################################### 
    #
    # Write out the file that will be the input for the transport step
    # ToDo: Note seeing yow we need a density array for this file
    #       Unless maybe we have materials of same name and different densities
    one_d_tool.write_spatial_transport_file(slab_length, slab_mat_name, slab_density)
    one_d_tool.transport_process()
    one_d_tool.response_process()
    outcome = one_d_tool.extract_results()
    # append_csv_file(filename, outcome)

    return
   
if __name__ == '__main__':
    main()

