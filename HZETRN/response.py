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
Also check and make sure there is a graveyard
"""
def find_ref_vol_and_check_graveyard(ref_point):
    # code snippet from dagmc.pyx:find_graveyard_inner_box
    volumes = dagmc.get_volume_list() 
    # ToDo: should -1 be used for next two lines?
    graveyard = 0
    ref_vol = 0 
    for v in volumes:
	# Actually, it doesn't matter if there isn't one
        if dagmc.volume_is_graveyard(v):
            graveyard = v
	    break
        if dagmc.point_in_volume(v, ref_point):
            ref_vol = v

    if graveyard == 0:
        print ('Could not find a graveyard volume')
        # raise DagmcError('Could not find a graveyard volume')

    # ToDo: error checking
    return (ref_vol, graveyard)

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
	# print material_obj
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
function to transform the tags into strings
tag : string of the tag to add to tag_list
tag_list : vector of tags in the problem
returns tag_list
"""
def tag_to_script(tag):
    a = []
    # since we have a byte type tag loop over the 32 elements
    for part in tag:
        # if the byte char code is non 0
        if (part != 0):
            # convert to ascii
            a.append(str(unichr(part)))
            # join to end string
            test = ''.join(a)
            # the the string we are testing for is not in the list of found
            # tag values, add to the list of tag_values

    return test

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
    material_entries = create_entries_from_lib(mat_lib)

    # Prepare xs_input.dat
    xs_input_filename = 'xs_input.dat'
    f = open(xs_input_filename, 'w')
    f.write(xs_header)
    f.write("\n".join(material_entries))
    f.close()

    # Using the input file just created, prepare the materials subdirectory
    # This method will make a subprocess call
    one_d_tool.cross_section_process(xs_input_filename)
    # End material cross-section part
    ##########################################################
    # Get the DAG object for this geometry
    rtn = dagmc.load(path)
    # get list of tag values
    (tag_values,tally_assigns) = tag_utils.get_tag_values(path)
    print ('#####################################################################')
    print tag_values, tally_assigns 
    print ('#####################################################################')
    """
    ##########################################################
    dag_vol_names = [] # list of dag volume names (Cubit id)
    dag_properties = set()
    # material tags
    dag_material_tags = []
    # tally tags
    dag_tally_tags=[]
    # tally assignments
    tally_assigns=[]

    # create imesh instance  
    dag_geom = iMesh.Mesh()
    # load the file
    dag_geom.load(path)

    # get all entities      
    ents = dag_geom.getEntities()
    # create a  mesh set         
    mesh_set = dag_geom.getEntSets()
    # list of volume ent handles 
    mat_list = []
    # get all geom_dimension ents
    geom_list = []
    cat_list  = []

    vol_tag  = dag_geom.getTagHandle('GEOM_DIMENSION')
    name_tag = dag_geom.getTagHandle('GLOBAL_ID')
    mat_tag  = dag_geom.getTagHandle('NAME')
    cat_tag  = dag_geom.getTagHandle('CATEGORY')

    # get the list we need    
    for i in mesh_set:
        tags = dag_geom.getAllTags(i)
        for tag in tags:
            if tag == vol_tag:
                geom_list.append(i)
            if tag == mat_tag:
                mat_list.append(i)
            if tag == cat_tag:
                cat_list.append(i)

    # for the 3d entities     
    for entity in geom_list:
        if vol_tag[entity] == 3:
             dag_vol_names.append(str(name_tag[entity]))

    # loop over all the volumes 
    for entity in geom_list:
        # loop over the material sets
        for meshset in mat_list:
            # if volume in set       
            if meshset.contains(entity):
                mat_name = mat_tag[meshset]
                volume_name = name_tag[entity]
                # dag_materials[volume_name]="".join( chr( val ) for val in mat_name )
                # dag_properties.add("".join( chr( val ) for val in mat_name))
                dag_properties.add(tag_to_script(mat_name))
                if 'tally:' in tag_to_script(mat_name):
                    pair = (volume_name,tag_to_script(mat_name))
                    tally_assigns.append(pair)

    # now we have dag properties, create one with materials and one with tallies
    for tag in dag_properties:
        if 'mat:' in tag:
            dag_material_tags.append(tag)
        # if 'tally:' in tag:
        #     dag_tally_tags.append(tag)
    print ('#####################################################################')
    print dag_vol_names, dag_properties
    print ('#####################################################################')
    """
    ##########################################################
    # get list of rays
    ray_tuples = load_ray_tuples(args.ray_dir_file) 

    # Use 0,0,0 as a reference point for now
    ref_point = [0.0, 0.0, 0.0]
    xyz = ref_point
    (start_vol, graveyard) = find_ref_vol_and_check_graveyard(ref_point)
    surf = 0
    dist = 0
    huge = 1000000000

    v = start_vol
    for dir in ray_tuples:
	# dagmc.tell_ray_story(ref_point, dir)
        slab_length = []
        slab_mat_name = []
        slab_density = []
        last_vol = start_vol
        for (vol, dist, surf, xyz) in dagmc.ray_iterator(start_vol, ref_point, dir, yield_xyz=True):
            print('  next intersection at distance', dist, 'on surface', surf)
	    print('  new xyz =', xyz)
	    print('proceeding into volume', vol)
	    if (dist < huge) and (surf != 0):
	        slab_length.append(dist)
		# need data from vol we are leaving
		md = dagmc.volume_metadata(last_vol)

		if 'material' in md:
		    slab_mat_name.append(md['material'])
		if 'rho' in md:
		    slab_density.append(md['rho'])
	        last_vol = vol
	    else:
	        vol = graveyard
        print('No more intersections for direction')
	print dir
	print slab_length
	print slab_mat_name
	print slab_density
      
    # End ray dir
    # now have filled
    # slab thicknesses (slab_length[])
    # slab materials (slab_mat_name[])
    # slab_densities (slab_density[]
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

