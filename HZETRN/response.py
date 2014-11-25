import subprocess
import argparse
import string

try:
   from pyne import material
   from pyne.material import Material, MaterialLibrary
except:
    raise ImportError("The PyNE dependencies could not be imported.")    

from pyne import dagmc
import numpy as np
import hzetrn as one_d_tool


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
	help='The path to the .h5m file')
    parser.add_argument(
        '-r', action='store', dest='ray_dir_file', 
	help='The path to the file with ray direction tuples')

    args = parser.parse_args

    if not args.uwuw_file:
        raise Exception('h5m file path not specified. [-f] not set')
    if not args.nuc_data:
        raise Exception('nuc_data file path not specified. [-d] not set')
    if not args.ray_dir_file:
        args.ray_dir_file = FALSE

    # ToDo: look at what this is used for
    if not args.output:
       args.output = 'dummy'
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
        if dagmc.volume_is_graveyard(v):
            graveyard = v
	    break
        if dagmc.point_in_volume(v, ref_point):
            ref_vol = v

    if graveyard == 0:
        raise DagmcError('Could not find a graveyard volume')

    # ToDo: error checking
    return (ref_vol, graveyard)
    
def xs_create_header():
    what1 = 'common_static_data'
    what2 = 'cross_sections_out'
    comment1 = '# Name of folder where static data is stored'
    comment2 = '# Name of folder to put data in (folder must already exist)'
    divider  = '{:-<84}'.format('-')
    lines = '{:<36}'.format(what1) + comment1 + '\n'
    lines = '{:<36}'.format(what2) + comment2 + '\n'
    lines += '\n'
    lines += divider + '\n'
    lines += '\n'

    return lines

def create_entries_from_lib(mat_lib):
    # Consider appending to string frist, for ease of testing
    # For each material create an entry and append to input_filename
    material_entry = []
    for key in mat_lib.iterkeys():
    	material_obj = mat_lib.get(key)
	# print material_obj.metadata['name']
	coll = material_obj.collapse_elements([])
        print coll
	print coll.metadata['name']
	material_entry.append[coll.metadata['name']]
	
	# material_entry = get_xs_info(material_obj)
	# ToDo: append material_entry to input_filename
	# append_entry(input_filename, material_entry)
    return material_entry

def main():
    # parse the the command line parameters
    args = parsing()

    # ToDo: Start the file with header lines which contain the names of 
    #       some folders the cross_section processing will need
    #       These may be hard-coded to start with; they can always be
    #       edited later
    xs_header = xs_create_header()

    # load the material library from the uwuw geometry file
    mat_lib = material.MaterialLibrary()
    mat_lib.from_hdf5(args.uwuw_file)
    material_entry = create_entries_from_lib(mat_lib)

    # Prepare xs_input.dat
    input_filename = 'xs_input.dat'
    f = open(input_filename, 'w')
    f.write(xs_header)
    # for item in material_entry
    #    f.write(item)
    f.close()

    # Using the input file just created, prepare the materials subdirectory
    # This method will make a subprocess call
    one_d_tool.cross_section_process(input_filename)
    # End material cross-section part
    ##########################################################
    # Get the DAG object for this geometry
    rtn = dagmc.load(uwuw_file)
   
    # get list of rays
    ray_tuples = load_ray_tuples(args.ray_dir_file) 

    # Use 0,0,0 as a reference point for now
    ref_point = [0.0, 0.0, 0.0]
    cur_point = ref_point
    (start_vol, graveyard) = find_ref_vol_and_check_graveyard(ref_point)
    surf = 0
    dist = 0
    huge = 1000000000
    slab_length = []
    slab_mat_name = []
    slab_density = []

    v = start_vol
    for dir in ray_tuples:
    	print dir	
	while v != graveyard:
            (surf,dist) = dagmc.ray_fire(v, cur_point, dir)
	    if (dist < huge) and (surf != 0):
	        slab_length.append(dist)
		md = dagmc.get_volume_metadata(v)
		slab_mat_name.append(md['material'])
		slab_density.append(md['rho'])
		# reassign v to next vol
		v = dagmc.next_vol(surf,v)
		# sqrt((dir[0]*dist)^2+dir[1]*dist)^2+ dir[2]*dist)^2)
		cur_point += dir*dist
	    else:
	        v = graveyard
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

