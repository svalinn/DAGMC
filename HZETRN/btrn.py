#########################################################
# btrn.py
# Backscatter Transport and Response
#########################################################

import subprocess
import argparse
import string
import os
import sys

try:
   from pyne.material import Material, MaterialLibrary
except:
    raise ImportError("The PyNE dependencies could not be imported.")    

try:
    from pyne import dagmc, data, nucname
except:
    raise ImportError("The DAGMC dependency could not be imported.")    

import numpy as np
import hzetrn as one_d_tool
import tag_utils


###############################
# globals
# 'spatial_dir' is a subdir for storing input files as they are created
spatial_dir = 'spatial'

def load_ray_tuples(filename):
    """ Load float tuples representing uvw directions from a file.
    
    filename : string
        filename of file containing pairs of angles (radians) representing polar and
	azimuthal angles. 

    Returns
    -------
    A list of float pairs representing the unit vector angles

    Notes
    ------
    In file with 10,000 tuples, here are the ranges
    ('min_theta', 0.00117047, 'max_theta', 3.136877355, 'min_phi', -3.141219364, 'max_phi', 3.14087921)

    """

    ray_tuples = []
    
    with open(filename) as f:
        for i, line in enumerate(f):
	    nums = map(float, line.split())
	    if len(nums) < 2:
	        msg = "Ray direction file {} requires two floats on each line; line {} is ill-formed".format(filename, i)
		sys.exit(msg)
	    w = np.cos(nums[0])
	    u = np.sin(nums[0])*np.cos(nums[1])
	    v = np.sin(nums[0])*np.sin(nums[1])
	    ray_tuples.append((u, v, w))
    return ray_tuples

def load_directions(ray_dir_file, num_ray_dirs):
    """ Load a numpy array with all the ray directions to be used.
    
    Parameters
    ----------
    ray_dir_file : string, filename 
        filename of file containing pairs of angles (radians) representing polar and
	azimuthal angles.  Does not have to be set.

    num_ray_dirs : integer
        Number of ray directions to use.  
	o If a filename is given and num_ray_dirs is
	  nonzero and less than the number of directions in the file, the return array 
	  is truncated to be of length num_ray_dirs.
	o If no filename is given and num_ray_dirs > 0, num_ray_dirs random
	  directions are generated
	o If no filename is given and num_ray_dirs = 0, three directions are
	  generated to be the Cartesian axes.
        

    Returns
    -------
    The set of directions is returned as a numpy row array
    """
    # Return the ray directions: use a cutoff if provided
    if ray_dir_file:
        ray_tuples = load_ray_tuples(ray_dir_file)
	if num_ray_dirs > 0:
	    ray_tuples = ray_tuples[:num_ray_dirs]

    elif num_ray_dirs > 0:
        ray_tuples = tag_utils.get_rand_dirs(num_ray_dirs)
 
    else:
        ray_tuples = [(1.0, 0.0, 0.0),
                      (0.0, 1.0, 0.0),
	   	      (0.0, 0.0, 1.0)]
    return np.array(ray_tuples)


def find_special_vols(ref_point, material_name_dict):
    """ Combine findng the vol-id of the geometry that contains the ref_point 
    and the vol-id whose name is 'graveyard'.  

    Requirements
    --------
    - The DAGMC python library must be loaded
    - Graveyard exists in the geometry

    Parameters
    ----------
    ref_point : float triplet
        Coordinate of point of interest.

    material_name_dict : dictionary, { string : int }
        A mapping of vol id to material object
	Key   = vol_id : int
	Value = unique material name

    Returns
    -------
    ref_vol : int
        The vol_id of the volume containing the ref_point.

    graveyard_vol : int
        The vol_id of the volume whose unique name is 'graveyard'.
	
    Exceptions
    ----------
        Raises exception if no volume is named 'graveyard'
	Raises exception if none of the volumes contain the reference point

    Notes
    -----
    This function has the drawback of going through every volume.  
    It could leave as soon as all (both) special volumes are found.
    """

    ref_vol = -1
    graveyard_vol = -1
    # code snippet from dagmc.pyx:find_graveyard_inner_box
    volumes = dagmc.get_volume_list() 
    for v in volumes:
        if dagmc.point_in_volume(v, ref_point):
	    ref_vol = v
	# Must check for inclusion in case it's the implicit complement, which
	# is not in the dict
	elif v in material_name_dict and material_name_dict[v] == 'graveyard':
	    graveyard_vol = v

    if ref_vol == -1:
        raise Exception("Reference volume not found!")
    if graveyard_vol == -1:
        raise Exception("Graveyard vol not found!  Please use a geometry file tagged with a graveyard.")

    return ref_vol, graveyard_vol

def find_slabs_for_ray(graveyard_vol, start_vol, ref_point, uvw, material_name_dict):
    """ Find the material distances and names along a direction
    in a tagged geometry.

    Requirements
    ------------
    graveyard_vol != start_vol

    Parameters
    ----------
    graveyard_vol : int
        The id of the graveyard.

    start_vol : int
        The id of the volume containing the ref_point.
    
    ref_point : float triplet
         The starting point from which to measure the slab lengths.

    uvw : float triplet
         Unit direction vector in the direction of travel

    material_name_dict : dictionary, { string : int }
        A mapping of vol id to material object
	Key   = vol_id : int
	Value = unique material name

    Returns
    -------
    slab_lengths : list of floats
        Distances traveled through materials ('slabs') encountered
    
    slab_mat_names : list of strings
        Unique names of slabs in the same order as they are encountered,
	meaning in the same order as the slab_lengths

    Notes
    -----
    This depends on DAGMC's ray_iterator, which calls ray_fire
    """

    surf = 0
    dist = 0
    # This geometrical hugeness, not integer hugeness
    huge = 1000000000

    slab_lengths = []
    slab_mat_names = []

    last_vol = start_vol
    for (next_vol, dist, surf, xyz) in dagmc.ray_iterator(last_vol, ref_point, uvw, yield_xyz=True):

        if last_vol in material_name_dict:
            slab_lengths.append(dist)
   	    # Get the unique name associated with volume
	    slab_mat_names.append(material_name_dict[last_vol])

	if next_vol == graveyard_vol or dist >= huge or surf == 0:
	    break;

        # If last_vol is the implicit complement, this is the only line executed
        last_vol = next_vol	    

	# Previous method: leave in for comparison for now
	"""
        is_graveyard = False    
        if not is_graveyard and ((dist < huge) and (surf != 0)):
	   
	    # We have already checked last_vol for being a graveyard
	    if last_vol in fname_for_vol:
	        slab_lengths.append(dist)
		# Get a name associated with volume
	        slab_mat_names.append(fname_for_vol[last_vol])

	    if next_vol in fname_for_vol:
   	        if fname_for_vol[next_vol] == 'graveyard':
	            is_graveyard = True

            # If last_vol is the implicit complement, this is the only line executed
            last_vol = next_vol

	# We hit the graveyard
        else:
	    break
	"""

    return slab_lengths, slab_mat_names
      
def parse_command_line_arguments():
    """Perform command line argument parsing
    Notes
    -----
    - Only one argument is absolutely required, the geometry file.
    - A run directory must have been previously set up with 
      HZETRN directories and whatever cross-sections are needed

    Questions
    ---------
    1. action = 'store' is default -- ok to remove?
    2. long-form keyword is also name of argparse parameter. 
       Underscores or camelCase?
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', action='store', dest='uwuw_file', 
	help='The relative path to the .h5m geometry file')

    parser.add_argument(
        '-d', '--rundir', action='store', 
	help='The name of the holding directory for all hzetrn runs')
    parser.set_defaults(rundir='rundir')

    parser.add_argument(
        '-e', '--environment',
	help='The radiation environment: spe (default) or gcr')
    parser.set_defaults(environment='spe')

    parser.add_argument(
        '-r', action='store', dest='ray_dir_file', 
	help='The path to the file with ray direction tuples')
    parser.set_defaults(ray_dir_file='')

    parser.add_argument(
        '-n', action='store', dest='num_ray_dirs', 
	help='The number of ray directions to use from file, \
	      or randomly create if no file given.  \
	      If neither file nor n is given, use Cartesion axes.', type=int)
    parser.set_defaults(num_ray_dirs=0)

    parser.add_argument(
        '-p', '--ref_point', nargs=3, type=float,
	help='Cartesion coordinates of the reference point for measurements')
    parser.set_defaults(ref_point=[0.0,0.0,0.0])

    parser.add_argument(
        '-s', '--response', action='store', 
        help='The material in which to calculate the response.  \
	      NB: the response material cross-section must exist')
    parser.set_defaults(response='WATERLIQ')

    args = parser.parse_args()

    # With no other parameters set byt the geometry file, a single randomly 
    # chosen direction will be run, using the defaults listed.
    if not args.uwuw_file:
        raise Exception('h5m file path not specified. [-f] not set')
    
    return args

def generateTransportInfo(slab_lens, slab_mat_names):
        """Extract and refine information from two lists.

	Parameters
	----------
        slab_lens : list of floats
	    Distances traveled through all non-graveyard, 
	    non-vacuum slabs in the geometry
        
	slab_mat_names : list of strings
	    Ordered list names of all non-vacuum, non-graveyard
	    fluka (unique) names from the tagged geometry.  
	    Matches the slab_lens list.

	Requirements
	------------
	slab_lens is not empty
	"""

	# Initialize list
	last_mat_name    = slab_mat_names[0]
	last_slab_points = [0.0, slab_lens[0]]

	for n in range(1,num_slabs):
            if slab_mat_names[n] == last_mat_name:
                last_slab_points.append(slab_lens[n] + slab_lens[n-1])
	    else:
	        tmp_mat_name    = last_mat_name
		tmp_slab_points = last_slab_points
                # update internal state before leaving
		last_mat_name = slab_mat_names[n]
		last_slab_points = [0.0, slab_lens[n]]
		yield tmp_mat_name, tmp_slab_points
	 
	    yield last_mat_name, last_slab_points

def createTransportDictionary(uvw, graveyard_vol, start_vol, ref_point, material_name_dict):
    """ Create a dictionary containing material names and distances for use in 
        one-dimensional transport code.

	Note
	----
	This is the first place we can know if the rays did not go through a material.
    """

    trn_dict = []
    # Rays fired here!
    slab_lens_a, slab_mat_names_a = find_slabs_for_ray(graveyard_vol, start_vol, ref_point,  uvw, material_name_dict)
    slab_lens_b, slab_mat_names_b = find_slabs_for_ray(graveyard_vol, start_vol, ref_point, -uvw, material_name_dict)

    if len(slab_lens_b) + len(slab_lens_a) == 0:
        return trn_dict, -1

    slab_lens      = slab_lens_b[:-1]      + slab_lens_a
    slab_mat_names = slab_mat_names_b[:-1] + slab_mat_names_a

    # If the ray does not go through any materials the yield is ('', [])
    for (mat_name, slab_points) in generateTransportInfo(slab_lens, slab_mat_names):
        trn_dict[mat_name] = slab_points

    return trn_dict, len(slab_lens_b)

def call_transport_response(run_path, uvw, trn_dict, env='spe', response='WATERLIQ'):
    """ Wrapper for call to one_d_tool module. 
    """
    one_d_tool.transport_response(run_path, uvw, trn_dict, env, response)
    return

def call_collect_results_for_dir(run_path, data_path, uvw, num_slabs_to_ref, index=0 ):
    one_d_tool.collect_results_for_dir(run_path, data_path, uvw, num_slabs_to_ref, index) 
    return

def main():
    """Main routine for btrn.py

    Notes
    -----
    The sequence of actions:
    - Ensure needed template files exist
    - Create a data directory in the run path
    - Create a temporary, named directory in the current directory
           to hold a collection of input files.
    - Load the given geometry file and get all its unique material (FLUKA) names
    - Get all the ray directions that are going to be used
    - Get the reference point and determine the reference volume
    - For each ray direction, starting from the reference vol
      o get the list of materials and distances it passes through in the given direction, up to the graveyard
      o get the list of materials and distances it passes through in the opposite direction, up to the graveyard
      o combine the materials and distances as if that of a ray passing through the reference point,
        which starts and ends at the graveyard
      o create the contents of the transport input file for this combined ray
    """

    # Setup: parse the the command line parameters
    args = parse_command_line_arguments()
    uwuw_file = os.path.join(os.getcwd(), args.uwuw_file)

    cur_path = os.getcwd() + '/'
    run_path = cur_path + args.rundir + '/'

    # Ensure needed template files exist
    one_d_tool.check_template_files(run_path)

    # Transport results for each direction will be placed in this directory:
    # Ensure it exists before proceeding.
    data_path = run_path + 'data/'
    if not os.path.isdir(data_path):
	os.mkdir(data_path)
    elif os.listdir(data_path) != []:
        raise Exception("Data directory {} is not empty! Please empty or delete it.".format(data_path))
        
    rad_env_file = data_path + 'rad_env.txt'
    with open(rad_env_file, 'w') as f:
        f.write(args.environment)

    # spatial_path = cur_path + spatial_dir + '/'
    # if not os.path.isdir(spatial_path):
    #	os.mkdir(spatial_path)

    # Load the DAG object for this geometry
    rtn = dagmc.load(uwuw_file)

    # ToDo:  Use a different function
    material_name_dict = tag_utils.get_fnames_for_vol(uwuw_file)
    ray_tuples = load_directions(args.ray_dir_file, args.num_ray_dirs)

    # The default starting point is 0,0,0
    ref_point = args.ref_point
    start_vol, graveyard_vol = find_special_vols(ref_point, material_name_dict)
    if start_vol == graveyard_vol:
        msg = "Exiting program: The reference point is in the graveyard!"
        sys.exit(msg)
    
    # Test function version
    index = 0
    for uvw in ray_tuples:
        index = index + 1
	trn_dict, number_slabs_to_reference_point =  \
	             createTransportDictionary(uvw, graveyard_vol, start_vol, \
		                               ref_point, material_name_dict)

        call_transport_response(run_path, trn_dict, args.environment, args.response)
        call_collect_results_for_dir(run_path, data_path, uvw, num_slabs_to_ref, index)

	
    """
    for uvw in ray_tuples:
	###################################
	# For backscatter:  
	# Create two lists 'a' and 'b'.  
	#    'a' is the list of slab materials and lengths in the given direction.  
	#    'b' is the list of same items going OPPOSITE (180 degrees)
	slab_lens_a, slab_mat_names_a = find_slabs_for_ray(graveyard_vol, start_vol, ref_point,  uvw, unique_material_name)
	slab_lens_b, slab_mat_names_b = find_slabs_for_ray(graveyard_vol, start_vol, ref_point, -uvw, unique_material_name)

        # slab_lens = slab_lens_a[::-1] + slab_lens_b
        slab_lens      = slab_lens_b[:-1]      + slab_lens_a
	slab_mat_names = slab_mat_names_b[:-1] + slab_mat_names_a
	num_slabs_to_ref = len(slab_lens_b)
	num_mats = len(slab_mat_names)
        #############################################
	# Create the transport geometry file contents
        # Need to have a non-zero number of mats
	transport_input = []
	if 0 != num_mats:
	    effective_num_mats = num_mats
	    last_mat_name = slab_mat_names[0]
	    last_slab_points = [0.0, slab_lens[0]]

	    for n in range(1,num_mats):
                if slab_mat_names[n] == last_mat_name:
		    #  Need to ADD the two lengths together: they are cumulative
		    last_slab_points.append(slab_lens[n] + slab_lens[n-1])
		    effective_num_mats = effective_num_mats - 1
		else:
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

	# No materials for this or the opposite direction: make bogus spatial file 
	else:
	    transport_input.append('0')
	    transport_input.append('No material traversed, either direction')
	    transport_input.append('2')
	    transport_input.append('0.0 0.0')

	# Need a newline at end of file
	transport_input.append('\n')

	dir_string = "{0:.4f}".format(uvw[0]) + '_' + \
	             "{0:.4f}".format(uvw[1]) + '_' + \
		     "{0:.4f}".format(uvw[2]) 

	# ToDo: change the base name here and on the other side
	if len(ray_tuples) < 20:
	    spatial_filename = 'spatial_ba' + dir_string + '.dat'
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
            one_d_tool.transport_process(run_path, spatial_filepath, args.environment)
            one_d_tool.response_process(run_path, args.response)

	one_d_tool.collect_results_for_dir(run_path, data_path, dir_string, num_mats, num_slabs_to_ref) 
        """
    # end for uvw in ray_tuples: 
    return 
    ###################################### 

   
if __name__ == '__main__':
    main()

