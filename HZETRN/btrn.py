#########################################################
# btrn.py
# Backscatter Transport and Respons
#########################################################

import subprocess
import argparse
import string
import os
import pprint

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
import logging

"""
To Run:

python btrn.py -f {4bricks.h5m, sphere_Wat_Al.h5m} [-r 10ThPhi.txt] [-e spe*/gcr] -d mdir
# for n random directions
python btrn.py -f {4bricks.h5m, sphere_Wat_Al.h5m} -n 50 -d kdir

"""

class DagmcError(Exception):
    pass

###############################
# globals
config_log  = 'config.log'
# 'spatial_dir' is a subdir for storing input files as they are created
spatial_dir = 'spatial'

def load_ray_start(filename):
    ray_start = []
    if not filename:
        ray_start = [0.0, 0.0, 0.0]
    else:
        with open(filename) as f:
	    # ToDo:  just do f.readline(), no for
	    for line in f:
	        ray_start = map(float, line.split())
    # pprint.pprint(ray_start)
    return ray_start

"""
Load or create tuples representing dirs
In file with 10,000 tuples, here are the ranges
('min_theta', 0.00117047, 'max_theta', 3.136877355, 'min_phi', -3.141219364, 'max_phi', 3.14087921)

filename name of file with rays in it, checked to exist
         before entering this function
"""
def load_ray_tuples(filename):

    ray_tuples = []
    
    with open(filename) as f:
        for line in f:
	    nums = map(float, line.split())
	    z = np.cos(nums[0])
	    x = np.sin(nums[0])*np.cos(nums[1])
	    y = np.sin(nums[0])*np.sin(nums[1])
	    ray_tuples.append((x, y, z))
    # print(filename, 'ray_tuples  x             y                z')
    # pprint.pprint(ray_tuples)
    return ray_tuples

"""
Load a subset of tuples representing directions
This loads directions within a tolerance of the z-plane
"""
def subset_ray_tuples(filename):
    side  = 10.0
    xlate = 20.0
    adj = xlate - side/2
    opp   = side/2
    limit = np.arctan(opp/adj)/2
    limitfac = 32.0

    ray_tuples = []

    # print('The ray tuples file is', filename)
    with open(filename) as f:
        for line in f:
	    # Get the list of numbers from the line
            nums = map(float, line.split())
	    theta = nums[0]
	    phi   = nums[1]
	    # limit limit by another factor of 2 to get a smaller subset
	    if (np.pi/2 - limit/limitfac < theta) and (theta < np.pi/2 + limit/limitfac):
	        z = np.cos(nums[0])
	        x = np.sin(nums[0])*np.cos(nums[1])
	        y = np.sin(nums[0])*np.sin(nums[1])
	        ray_tuples.append((x, y, z))
        # print ('num of rays within 1/', limitfac, 'of x-y plane', len(ray_tuples))
    return ray_tuples

"""
Return an array of directions
"""
def get_directions(args):
    # if a filename is given get the ray directions or a subset from it
    if args.ray_dir_file:
        if args.ray_subset:
            # Get all the rays from a large ray file that are
            # within a few degrees of the xy plane 
            ray_tuples = np.array(subset_ray_tuples(args.ray_dir_file))
	else:
            ray_tuples = np.array(load_ray_tuples(args.ray_dir_file))

    elif args.rand_dirs > 0:
	# print 'Getting ', args.rand_dirs, ' random directions'
        ray_tuples = tag_utils.get_rand_dirs(args.rand_dirs)
 
    else:
        ray_tuples = np.array( [(1.0, 0.0, 0.0),
                                (0.0, 1.0, 0.0),
	   	                (0.0, 0.0, 1.0)] )
    return ray_tuples


""" 
Find the vol-id of the geometry that contains the ref_point
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
Do all the work of finding the slab distances and material names
start_vol  The volume containing the ref_point
ref_point  the origin of the ray
dir 	   the direction of the ray
fname_for_vol a mapping of vol id to material object

returns    arrays  for length and material name as the ray
           travels from the start point to the graveyard

This depends on dagmc's ray_iterator, which calls ray_fire
"""
def slabs_for_ray(start_vol, ref_point, dir, fname_for_vol):
    is_graveyard = False
    foundMat = False
    surf = 0
    dist = 0
    huge = 1000000000

    slab_lengths = []
    slab_mat_names = []
    vols_traversed = []

    # print 'start_vol is', start_vol
    # Could be starting in the implicit_complement
    if start_vol in fname_for_vol:
        #  Done before we've started
        if fname_for_vol[start_vol] == 'graveyard':
	    # print 'Starting volume is graveyard'
            return slab_length, slab_mat_name 
    # else:
        # print 'We are starting in the implicit complement'
    is_graveyard = False    
    last_vol = start_vol
    
    for (next_vol, dist, surf, xyz) in dagmc.ray_iterator(last_vol, ref_point, dir, yield_xyz=True):
        # print 'next_vol', next_vol
        if not is_graveyard and ((dist < huge) and (surf != 0)):
	    # We have already checked last_vol for being a graveyard
	    if last_vol in fname_for_vol:
	        slab_lengths.append(dist)
		# Get a name associated with volume
	        slab_mat_names.append(fname_for_vol[last_vol])
		# Optional documentation
                vols_traversed.append(last_vol)

	    if next_vol in fname_for_vol:
   	        if fname_for_vol[next_vol] == 'graveyard':
	            is_graveyard = True

            # If last_vol is the implicit complement, this is the only line executed
            last_vol = next_vol

	# We hit the graveyard
        else:
	    break

    #########################################################################
    sdir = [format(dir[0],'.6f'),format(dir[1], '.6f'), format(dir[2], '.6f')]
    # print ('for dir', sdir, 'vols traversed', vols_traversed)
    #########################################################################

    return slab_lengths, slab_mat_names
      
"""
Argument parsing
returns : args: -f for the input geometry file, -r for the input ray tuple file
ref     : DAGMC/tools/parse_materials/dagmc_get_materials.py
"""
def parsing():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', action='store', dest='uwuw_file', 
	help='The relative path to the .h5m geometry file')

    parser.add_argument(
        '-d', action='store', dest='run_dir',
	help='The name of the holding directory for all hzetrn runs')

    parser.add_argument(
        '-e', action='store', dest='rad_env', 
	help='The radiation environment: spe (default) or gcr')
    parser.set_defaults(rad_env='spe')

    parser.add_argument(
        '-r', action='store', dest='ray_dir_file', 
	help='The path to the file with ray direction tuples')
    parser.set_defaults(ray_dir_file='')

    parser.add_argument(
        '-n', action='store', dest='rand_dirs', 
	help='The number of random ray directions to use', type=int)

    parser.add_argument(
        '-p', action='store', dest='ray_start', 
	help='File containing Cartesion coordinate of starting point of all rays')

    parser.add_argument(
        '--ray_subset', dest='ray_subset', action='store_true',
        help='Use only rays within a small angle of the xy-plane')
    parser.add_argument(
        '--no_ray_subset', dest='ray_subset', action='store_false',
        help='Use all rays (default)')
    parser.set_defaults(ray_subset=False)

    parser.add_argument(
        '-target', action='store', dest='target',
        help='The target material in which to calculate the response.  \
	      NB: the target material cross-section must exist'    )
    parser.set_defaults(target='')

    args = parser.parse_args()

    if not args.uwuw_file:
        raise Exception('h5m file path not specified. [-f] not set')
    
    return args

""" 
    btrn.py
    - Log this call
    - Ensure needed template files exist
    - Create a data directory in the run path
    - Create a temporary, named directory in the current directory
           to hold a collection of in put files.
    - Load the given geometry file and get all its FLUKA names
    - Get all the ray directions we are going to use
    - Get the reference point and determine the start volume
    - For each ray direction, starting from the start vol
      o get the list of materials and distances it passes through
      o create the contents of the transport input file
"""
def main():
    global spatial_dir
    global config_log

    # Setup: parse the the command line parameters
    args = parsing()
    uwuw_file = os.path.join(os.path.dirname('__file__'), args.uwuw_file)

    cur_path = os.path.dirname(os.path.abspath(__file__)) + '/'
    run_path = cur_path + args.run_dir + '/'
    config_path = run_path + config_log

    logging.basicConfig(filename=config_path, level=logging.INFO, format='%(asctime)s %(message)s')
    message = 'python btrn.py -f ' + args.uwuw_file + \
              ' -d ' + args.run_dir + \
              ' -e ' + args.rad_env + \
              ' -r ' + ('not given' if not args.ray_dir_file else args.ray_dir_file)
    logging.info(message)

    # Ensure needed template files exist
    one_d_tool.check_temp_files(run_path)

    # Transport results for each direction will be placed in this directory:
    # Ensure it exists before proceeding.
    data_path = run_path + 'data/'
    if not os.path.isdir(data_path):
        # print 'Creating data path', data_path
	os.mkdir(data_path)

    rad_env_file = data_path + 'rad_env.txt'
    with open(rad_env_file, 'w') as f:
        f.write(args.rad_env)

    spatial_path = cur_path + spatial_dir + '/'
    if not os.path.isdir(spatial_path):
        # print 'Creating path', spatial_path
	os.mkdir(spatial_path)

    # Load the DAG object for this geometry
    rtn = dagmc.load(uwuw_file)

    vol_fname_dict = tag_utils.get_fnames_for_vol(uwuw_file)
    # print 'fnames_dict', vol_fname_dict
    ray_tuples = get_directions(args)

    # The default starting point is 0,0,0
    ref_point = load_ray_start(args.ray_start)
    # print 'ref_point', ref_point
    start_vol = find_ref_vol(ref_point)
    
    for dir in ray_tuples:
	###################################
	# For backscatter:  
	# Create two lists 'a' and 'b'.  
	#    'a' is the list of slab materials and lengths in the given direction.  
	#    'b' is the list of same items going OPPOSITE (180 degrees)
	slab_lens_a, slab_mat_names_a = slabs_for_ray(start_vol, ref_point,  dir, vol_fname_dict)
	slab_lens_b, slab_mat_names_b = slabs_for_ray(start_vol, ref_point, -dir, vol_fname_dict)

	# reverse order and add backward dir
        # slab_lens = slab_lens_a[::-1] + slab_lens_b
        slab_lens      = slab_lens_b[::-1]      + slab_lens_a
	slab_mat_names = slab_mat_names_b[::-1] + slab_mat_names_a
	num_slabs_to_ref = len(slab_lens_a)
	# print 'slab_lens', slab_lens, 'names', slab_mat_names, 'num to ref', num_slabs_to_ref
	num_mats = len(slab_mat_names)
        #############################################
	# Create the transport geometry file contents
        # Need to have a non-zero number of mats
	transport_input = []
	if 0 != num_mats:
	    effective_num_mats = num_mats
	    # last_mat_name = ''
	    # last_slab_points = []
	    last_mat_name = slab_mat_names[0]
	    last_slab_points = [0.0, slab_lens[0]]

	    # print 0, last_mat_name, slab_lens[0]
	    for n in range(1,num_mats):
                if slab_mat_names[n] == last_mat_name:
	            # print n, 'SAME\t', last_mat_name, slab_lens[n]
		    #  Need to ADD the two lenghs together: they are cumulative
		    last_slab_points.append(slab_lens[n] + slab_lens[n-1])
		    effective_num_mats = effective_num_mats - 1
		else:
		    # print n, 'NEW\t', slab_mat_names[n], slab_lens[n]
	            # print '\tWriting\t', last_mat_name, last_slab_points
	            transport_input.append(last_mat_name)
	            transport_input.append(str(len(last_slab_points)))
		    slab_points_list = ["{0:.1f}".format(x) for x in last_slab_points]
		    transport_input.append(", ".join(slab_points_list))
                    # Now re-initialize to the current mat/length
		    last_mat_name = slab_mat_names[n]
		    last_slab_points = [0.0, slab_lens[n]]
	    # The last slab got set up but not written
	    transport_input.append(last_mat_name)
	    transport_input.append(str(len(last_slab_points)))
	    slab_points_list = ["{0:.1f}".format(x) for x in last_slab_points]
	    transport_input.append(", ".join(slab_points_list))
	        
	    transport_input.insert(0,(str(effective_num_mats)))

	# No materials at all for this or the opposite direction: make a bogus spatial file 
	else:
	    transport_input.append('0')
	    transport_input.append('No material traversed, either direction')
	    transport_input.append('2')
	    transport_input.append('0.0 0.0')

	# Need a newline at end of file
	transport_input.append('\n')
        # print 'ray_tuple', dir, 'transport_input', transport_input

	dir_string = "{0:.4f}".format(dir[0]) + '_' + \
	             "{0:.4f}".format(dir[1]) + '_' + \
		     "{0:.4f}".format(dir[2]) 

	# ToDo: change the base name here and on the other side
	if len(ray_tuples) < 20:
	    spatial_filename = 'spatial_b' + dir_string + '.dat'
	    # sslab = []
            # for d in slab_lengths:
	    #     sslab.append(format(d,'.6f'))
	    # print ('dist, mats ', sslab, slab_mat_names)
            #############################################
	else:
	    spatial_filename = 'spatial.dat'

	# Write out the spatial file even for no materials traversed
	# A local subdirectory is used.  This is currently 
	# hardcoded to ./spatial, but a temporary directory could also be used.
        spatial_filepath = spatial_path + spatial_filename
	f = open(spatial_filepath, 'w')
	f.write("\n".join(transport_input))
	f.close()

	if 0 != num_mats:
            one_d_tool.transport_process(run_path, spatial_filepath, args.rad_env)
            one_d_tool.response_process(run_path, args.target)

	one_d_tool.collect_results_for_dir(run_path, data_path, dir_string, num_mats, num_slabs_to_ref) 
	# if num_slabs_to_ref == 1:
	#     print "For direction", dir_string, ", one slab to graveyard"
    return 
    ###################################### 

   
if __name__ == '__main__':
    main()

