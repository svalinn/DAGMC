import subprocess
import argparse
import string
import os
import pprint
import ConfigParser

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

# 'spatial_dir' is a subdir for storing input files as they are created
spatial_dir = 'spatial'

def load_ray_start(filename):
    ray_start = []
    if not filename:
        ray_start = [0.0, 0.0, 0.0]
    else:
        with open(filename) as f:
	    for line in f:
	        ray_start = map(float, line.split())
    pprint.pprint(ray_start)
    return ray_start

"""
Load or create tuples representing dirs
In file with 10,000 tuples, here are the ranges
('min_theta', 0.00117047, 'max_theta', 3.136877355, 'min_phi', -3.141219364, 'max_phi', 3.14087921)
"""
def load_ray_tuples(filename, number):

    ray_tuples = []
    if number:
        ray_tuples = get_rand_dirs(number)
	return ray_tuples
    
    # ToDo: randomly generate dirs
    if not filename:
       ray_tuples = [(1.0, 0.0, 0.0),
                     (0.0, 1.0, 0.0),
	   	     (0.0, 0.0, 1.0)]
    else:
       with open(filename) as f:
           for line in f:
	       nums = map(float, line.split())
	       z = np.cos(nums[0])
	       x = np.sin(nums[0])*np.cos(nums[1])
	       y = np.sin(nums[0])*np.sin(nums[1])
	       ray_tuples.append((x, y, z))
       print(filename, 'ray_tuples  x             y                z')
       pprint.pprint(ray_tuples)
    return ray_tuples

"""
Load a subset of tuples representing directions
This loads directions within a tolerance of the z-plane
"""
def subset_ray_tuples(filename):
    side = 10.0
    xlate = 20.0
    adj = xlate - side/2
    opp   = side/2
    limit = np.arctan(opp/adj)/2
    limitfac = 32.0

    ray_tuples = []

    print('The ray tuples file is', filename)
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
        print ('num of rays within 1/', limitfac, 'of x-y plane', len(ray_tuples))
    return ray_tuples

"""
Produce a length-one tuple randomly oriented on a sphere
"""
def get_rand_dirs(number):
    rays = []
    if number > 0:
        for i in range(number):
            # [0.0,1.0)
            rnum = np.random.uniform()
            # map rnum onto [0.0,1.0)
            z = 2*rnum - 1
            theta = 2*np.pi*rnum
            norm_fac = np.sqrt(1 - z*z)
            y = norm_fac*np.sin(theta)
            x = norm_fac*np.cos(theta)
	    rays.append([x, y, z])
    return rays


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

    print 'start_vol is', start_vol
    # Could be starting in the implicit_complement
    if start_vol in fname_for_vol:
        #  Done before we've started
        if fname_for_vol[start_vol] == 'graveyard':
	    print 'Starting volume is graveyard'
            return slab_length, slab_mat_name 
    else:
        print 'We are starting in the implicit complement'
    is_graveyard = False    
    last_vol = start_vol
    
    for (next_vol, dist, surf, xyz) in dagmc.ray_iterator(last_vol, ref_point, dir, yield_xyz=True):
        print 'next_vol', next_vol
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
    print ('for dir', sdir, 'vols traversed', vols_traversed)
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
	help='The relative path to the .h5m file')
#    parser.add_argument(
#        '-c', action='store', dest='config_file',
#	help='The name of the config file')
    parser.add_argument(
        '-d', action='store', dest='run_dir',
	help='The name of the holding directory for all hzetrn runs')

    parser.add_argument(
        '-r', action='store', dest='ray_dir_file', 
	help='The path to the file with ray direction tuples')

    parser.add_argument(
        '-n', action='store', dest='rand_dirs', 
	help='The number of randome ray directions to use', type=int)

    parser.add_argument(
        '-p', action='store', dest='ray_start', 
	help='Cartesion coordinate of starting point of all rays')

    args = parser.parse_args()

    if not args.uwuw_file:
        raise Exception('h5m file path not specified. [-f] not set')
#    if not args.config_file:
#        args.config_file = 'hze.cfg'
    if not args.ray_dir_file and not args.rand_dirs:
        args.ray_dir_file = False
       
    return args

def main():
    global spatial_dir

    # Setup: parse the the command line parameters
    args = parsing()
    path = os.path.join(os.path.dirname('__file__'), args.uwuw_file)

    curdir = os.path.dirname(os.path.abspath(__file__)) + '/'
    run_path = curdir + args.run_dir + '/'

    spatial_path = curdir + spatial_dir + '/'
    if not os.path.isdir(spatial_path):
        print 'Creating path', spatial_path
	os.mkdir(spatial_path)

    # Load the DAG object for this geometry
    rtn = dagmc.load(path)

    vol_fname_dict = tag_utils.get_fnames_for_vol(path)
    print 'fnames_dict', vol_fname_dict
    # get list of rays
    # ray_tuples = load_ray_tuples(args.ray_dir_file, args.rand_dirs) 
    ray_tuples = subset_ray_tuples(args.ray_dir_file) 

    # The default starting point is 0,0,0
    ref_point = load_ray_start(args.ray_start)
    print 'ref_point', ref_point
    start_vol = find_ref_vol(ref_point)
    
    response_filepath = run_path + 'collect_doseq_table.csv' 
    with open(response_filepath,'a') as f:
        f.write('\'' + path+'\'\n')
       
    # i=1
    for dir in ray_tuples:
	slab_lengths, slab_mat_names = slabs_for_ray(start_vol, ref_point, dir, vol_fname_dict)
        #############################################

	num_mats = len(slab_mat_names)
	print 'num_mats', num_mats

        # Need to have a non-zero number of mats
	transport_input = []
	transport_input.append(str(num_mats))
	if 0 != num_mats:
	    for n in range(num_mats):
	        transport_input.append(slab_mat_names[n])
	        transport_input.append('2')
	        transport_input.append('0.0 ' + "{0:.1f}".format(slab_lengths[n]))
	# Make a bogus spatial file, for the record
	else:
	    transport_input.append('No material traversed')
	    transport_input.append('2')
	    transport_input.append('0.0 0.0')

	# Need a newline at end of file
	transport_input.append('\n')
        print 'ray_tuple', dir, 'transport_input', transport_input
	
	dir_string = "{0:.4f}".format(dir[0]) + '_' + \
	             "{0:.4f}".format(dir[1]) + '_' + \
		     "{0:.4f}".format(dir[2]) 


	if len(ray_tuples) < 20:
	    spatial_filename = 'spatial_' + dir_string + '.dat'
	    sslab = []
            for d in slab_lengths:
	        sslab.append(format(d,'.6f'))
	    print ('dist, mats ', sslab, slab_mat_names)
            #############################################
	else:
	    spatial_filename = 'spatial.dat'

	# if i < 50:
	# Write out the spatial file even for no materials traversed
        trn_path = spatial_path + spatial_filename
	f = open(trn_path, 'w')
	f.write("\n".join(transport_input))
	f.close()
	# else:
	if 0 != num_mats:
            one_d_tool.transport_process(run_path, trn_path)
            one_d_tool.response_process(run_path)

	data_line = one_d_tool.collect_results_for_dir(run_path, dir_string, trn_path)

        print 'response_filepath:', response_filepath
        with open(response_filepath,'a') as f:
            f.write(data_line)

    with open(response_filepath) as f:
        f.close()
    return 
	
	    
# 	i = i + 1
    ###################################### 
    #
    # Write out the file that will be the input for the transport step
    # one_d_tool.write_spatial_transport_file(slab_length, slab_mat_name)
    # append_csv_file(filename, outcome)

    return
   
if __name__ == '__main__':
    main()

