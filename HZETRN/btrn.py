#########################################################
# btrn.py
# Backscatter Transport and Response
#########################################################

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

import uwuw_preproc
import numpy as np
from hzetrn import transportresponse


def get_rand_dirs(number):
    """ Produce length-one tuples randomly oriented on a sphere.

    Parameter
    ---------
    number : int
        The number of tuples to be returned.

    Return
    ------
    A two-dimensional numpy array with rows consisting of the 
    ray.
    """
    rays = []
    if number > 0:
        for i in range(number):
            # [0.0,1.0)
            rnum = np.random.uniform()
            # map rnum onto [0.0,1.0)
            z = 2*rnum - 1

	    # Call a new randome number
            theta = 2*np.pi*np.random.uniform()
            norm_fac = np.sqrt(1 - z*z)
            y = norm_fac*np.sin(theta)
            x = norm_fac*np.cos(theta)
	    rays.append([x, y, z])
    return rays

def load_ray_tuples(filename):
    """ Load float tuples representing uvw directions from a file.
    
    Parameters
    ----------
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
                msg = "Ray direction file {} requires two floats on each line; \
                       line {} is ill-formed".format(filename, i)
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
        print 'num_ray_dirs', num_ray_dirs
        ray_tuples = load_ray_tuples(ray_dir_file)
        if num_ray_dirs > 0:
            ray_tuples = ray_tuples[:num_ray_dirs]
        elif num_ray_dirs < 0:
	    ray_tuples = ray_tuples[-num_ray_dirs:]
    elif num_ray_dirs > 0:
        ray_tuples = get_rand_dirs(num_ray_dirs)
    else:
        ray_tuples = [(1.0, 0.0, 0.0),
                      (0.0, 1.0, 0.0),
                      (0.0, 0.0, 1.0)]
    return np.array(ray_tuples)


def find_ref_vol(ref_point):
    """ Find the vol-id of the geometry that contains the ref_point.

    Requirements
    ------------
    - The DAGMC python library must be loaded

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

    Exceptions
    ----------
    Raises exception if none of the volumes contain the reference point
    """

    ref_vol = -1
    # code snippet from dagmc.pyx:find_graveyard_inner_box
    volumes = dagmc.get_volume_list() 
    for v in volumes:
        if dagmc.point_in_volume(v, ref_point):
            ref_vol = v
	    break

    if ref_vol == -1:
        raise Exception("Reference volume not found!")

    return ref_vol 

def find_slabs_for_ray(graveyard_vol, start_vol, ref_point, uvw, vol_name_dens):
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
    slab_areal_densities = []
    slab_ids = []

    last_vol = start_vol
    for (next_vol, dist, surf, xyz) in dagmc.ray_iterator(last_vol, ref_point, uvw, yield_xyz=True):
        # if last_vol in material_name_dict:
        if last_vol in vol_name_dens:
            slab_lengths.append(dist)
	    name = vol_name_dens[last_vol][0]
	    dens = vol_name_dens[last_vol][1]
            slab_mat_names.append(name)
	    slab_areal_densities.append(dist*dens)

	    # Get the list of vol-id's in order. 
	    slab_ids.append(last_vol)

        if next_vol == graveyard_vol or dist >= huge or surf == 0:
            break;

        # If last_vol is the implicit complement, this is the only line executed
        last_vol = next_vol        

    return slab_mat_names, slab_areal_densities, slab_lengths, slab_ids
    # End of find_slabs_for_ray
      
def spatialTuples(slab_mat_names, slab_depthq):
    """ Implement an algorithm to create tuples from the material/distances
    lists that are in a convenient form for creating transport execution.

    Parameters
    ----------
    slab_mat_names : list of strings
        Unique material names from the geometry, in the order they are 
    encountered for a particular ray.
    
    slab_depthq : list of floats
        Areal density: Distances traveled through each material by the ray multiplied 
    by the density of that material, in the same order as slab_mat_names

    Return
    ------
    spatial : list of tuples
        Tuples are (string, list of floats)
    First element of list of floats is always zero
    List of floats is length 2 at a minimum
    """
    spatial = []
    last_mat_name = slab_mat_names[0]
    current_slab_pts = [0.0]
    # for mat_name, slab_thick in zip(slab_mat_names, slab_lens):
    for mat_name, slab_depthq in zip(slab_mat_names, slab_depthq):
        if mat_name != last_mat_name:
            # Every time we get a new mat, store the name/slab points tuple...
            spatial.append((last_mat_name, current_slab_pts))
            # ... and re-initialize
            last_mat_name = mat_name
            current_slab_pts = [0.0]
        # Do always; first time through, this is all that is done
        current_slab_pts.append(current_slab_pts[-1] + slab_depthq)

    # Once more for the trailing tuple 
    spatial.append((last_mat_name, current_slab_pts))

    return spatial

def createTransportDictionary(uvw, graveyard_vol, start_vol, ref_point, vol_name_dens, info_type):
    """ Create a dictionary containing material names and equivalent 
    distances (length times density) for use in one-dimensional transport code.

    Parameters
    ----------
    uvw : list of 3 floats
        The direction vector of a unit vector.  The final ray going through 
    the material points in this direction.

    graveyard_vol : int
        Vol id of the graveyard, a convenience value

    start_vol : int
        Vol id of the volume containing the reference point.

    ref_point : list of 3 floats
        Triplet with the location of the point at which data is desired

    vol_name_dens : dictionary
        Mapping from each volume in the geometry to a tuple consisting
	of (unique "fluka" name created by uwuw_preproc, density of that material)
       
    Return
    ------
    spatial : list of tuples
        Tuples are (string, list of floats)
    First element of list of floats is always zero
    List of floats is length 2 at a minimum

    Number of slabs the ray goes through to get to the reference point.
    
    A structured array record giving the volume id, name, and density
        at the reference point, plus the cumulative depth and areal
	density to the reference point

    Note
    ----
    This is the first place we can know if the rays did not go through a 
    material.
    """

    # Rays fired here!
    ###################################
    # Create two lists 'a' and 'b'.  
    #    'a' is the list of slab materials and lengths in the given direction.  
    #    'b' is the list of same items going OPPOSITE (180 degrees)
    slab_mat_names_a, areal_densities_a, slab_lens_a, slab_ids_a = \
        find_slabs_for_ray(graveyard_vol, start_vol,  ref_point,   \
                           uvw, vol_name_dens)

    slab_mat_names_b, areal_densities_b, slab_lens_b, slab_ids_b = \
        find_slabs_for_ray(graveyard_vol, start_vol, ref_point,    \
                          -uvw, vol_name_dens)

    # The ray did not go through any material
    info_at_ref = np.zeros(1, dtype=info_type)
    if len(slab_lens_b) + len(slab_lens_a) == 0:
        return [], -1, info_at_ref

    slab_lens       = slab_lens_b[::-1]       + slab_lens_a
    slab_mat_names  = slab_mat_names_b[::-1]  + slab_mat_names_a
    areal_densities = areal_densities_b[::-1] + areal_densities_a
    slab_ids        = slab_ids_b[::-1]        + slab_ids_a
    

    ref_slab_depth = len(slab_lens_b)
    ref_id      = slab_ids[ref_slab_depth]
    ref_name    = slab_mat_names[ref_slab_depth]
    ref_density = vol_name_dens[ref_id][1]
    ref_depth   = np.sum(slab_lens[:ref_slab_depth])
    print "slab_lens to ref", slab_lens[:ref_slab_depth], "sum", ref_depth
    ref_depthq  = np.sum(np.around(areal_densities[:ref_slab_depth],1))
    print "areal_densities to ref", areal_densities[:ref_slab_depth], "sum", ref_depthq
    
    info_at_ref[0] = (ref_id, ref_name, ref_density, ref_depth, ref_depthq)
    return spatialTuples(slab_mat_names, areal_densities), ref_slab_depth,  info_at_ref
    # End of createTransportDictionary

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

def main():
    """Main routine for btrn.py

    Notes
    -----
    The sequence of actions:
    - Parse the command line
    - Create a data directory in the run path
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
    
    args = parse_command_line_arguments()
    uwuw_file = os.path.join(os.getcwd(), args.uwuw_file)
    run_path  = os.path.join(os.getcwd(), args.rundir)

    # Transport results for each direction will be placed in this directory:
    # Ensure it exists before proceeding, and does not contain data 
    # subdirectories from previous runs (files are ok)
    data_path = os.path.join(run_path,'data')
    # If a negative number of ray_dirs is given we are in continuatio nmode
    if 0 <= args.num_ray_dirs:
        if not os.path.isdir(data_path):
            os.mkdir(data_path)
        else:
            for el in os.listdir(data_path):
                el_path = os.path.join(data_path, el)
                if os.path.isdir(el_path):
                    raise Exception("Data directory {} has subdirectories! \
		                     Remove these and continue.". \
				     format(data_path))
    
    # Initialize the transportresponse class
    one_d_tool = transportresponse(run_path, data_path, \
                                   args.environment, args.response)

    #  ToDo: create this in transportresponse
    info_type = [('id','i4'),('name','a12'),('density', 'f4'), ('depth','f4'), ('depthq','f4')]
    # Load the DAG object for this geometry
    dagmc.load(uwuw_file)

    name_density_dict, graveyard_vol = uwuw_preproc.get_fnames_and_densities(uwuw_file)
    rays_file  = os.path.join(os.getcwd(), args.ray_dir_file)
    ray_tuples = load_directions(rays_file, args.num_ray_dirs)

    # The default starting point is 0,0,0
    ref_point = args.ref_point
    start_vol = find_ref_vol(ref_point)

    if start_vol == graveyard_vol:
        msg = "Exiting program: The reference point is in the graveyard!"
        sys.exit(msg)
    
    if 0 <= args.num_ray_dirs:
        index = 1
    else:
        index = -args.num_ray_dirs

    for uvw in ray_tuples:
        spatial_tuples, number_slabs_to_reference_point, info_tuple =  \
                 createTransportDictionary(uvw, graveyard_vol, start_vol, \
                                           ref_point, name_density_dict, info_type)
        one_d_tool.execute(spatial_tuples)
	
	for key in name_density_dict:
	    print "vol, name, density", key, name_density_dict[key][0], name_density_dict[key][1]
        one_d_tool.collect_results(index, uvw, number_slabs_to_reference_point, info_tuple) 
        index = index + 1
    return 
   
if __name__ == '__main__':
    main()

