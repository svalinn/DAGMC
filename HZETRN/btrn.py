#########################################################
# btrn.py
# Transport and Response
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

##############################################################################
#                  Alltran
##############################################################################
class alltran(object):
    """ Class that manages the data and attributes for a geometry/ray-direction
    problem
    """
    uwuw_filepath = ''
    runpath  =  ''
    datapath = '' 

    numdirs = -1
    # Default attributes help with testing
    environment = 'spe'
    reference_point = np.array([0.0,0.0,0.0])
    response    = 'WATERLIQ'

    # int : The vol_id of the volume containing the reference_point.
    reference_vol   = -1
    # Vol id of the graveyard, a convenience value
    graveyard_vol = -1

    slab_info_type =  [('id', 'i4'), ('name','a12'),('depth', 'f4'), \
                       ('density', 'f4'), ('areal_density', 'f4')]

    # A mapping of vol id to  Value = (unique_material_name, density)
    name_for_vol = {}
    # dtype for the returned info record array
    info_type = []

    def __init__(self, uwuw, rundir, n):
        """ Instantiate the class and set up the data directories
        """

	self.numdirs = n
        self.set_paths(uwuw, rundir)

	return

    def set_paths(self, uwuw, rundir):
        """ Create a data directory in the run path
	"""
        run_path = os.path.join(os.getcwd(), rundir)

        # Transport results for each direction will be placed in this directory:
        # Ensure it exists before proceeding, and does not contain data 
        # subdirectories from previous runs (files are ok)
        data_path = os.path.join(run_path,'data')
        # If a negative number of ray_dirs is given we are in continuatio nmode
        if 0 <= self.numdirs:
            if not os.path.isdir(data_path):
                os.mkdir(data_path)
            else:
                for el in os.listdir(data_path):
                    el_path = os.path.join(data_path, el)
                    if os.path.isdir(el_path):
                        raise Exception("Data directory {} has subdirectories! \
		                         Remove these and continue.". \
				         format(data_path))

	self.uwuw_filepath = os.path.join(os.getcwd(), uwuw)
	self.runpath   = run_path
	self.datapath  = data_path
        return 

    def get_rays_and_reference_vol(self, ray_dir_file): 
        """ Read the ray directions from the file and perform checking

	Parameters
	----------
	ray_dir_file : string
	    The name of the file

        Returns
	-------
	index
	    Normally 1, but if given a negative input parameter for 
	    numdirs, will start at -numdirs.  This is to restart without
	    having to clear previously good data directories

	ray_tuples
	    Floating point triplet for the ray direction vectors, read
	    from the file
        """

        rays_file  = os.path.join(os.getcwd(), ray_dir_file)
        ray_tuples = self.load_directions(rays_file, self.numdirs)

        # The default starting point is 0,0,0
        self.reference_vol = self.find_reference_vol()
	print "references: point & vol", self.reference_point, self.reference_vol

        if self.reference_vol == self.graveyard_vol:
            msg = "Exiting program: The reference point is in the graveyard!"
            sys.exit(msg)
    
        if 0 <= self.numdirs:
            index = 1
        else:
            index = -self.numdirs
        return index, ray_tuples

    def get_rand_dirs(self, number):
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

    def load_ray_tuples(self, filename):
        """ Load float tuples representing uvw directions from a file.
    
        Parameters
        ----------
        filename : string
            filename of file containing pairs of angles (radians) 
	    representing polar and azimuthal angles. 

        Returns
        -------
        A list of float pairs representing the unit vector angles

        Notes
        ------
        In file with 10,000 tuples, here are the ranges
        ('min_theta', 0.00117047, 'max_theta', 3.136877355, 
	 'min_phi', -3.141219364, 'max_phi', 3.14087921)

        """

        ray_tuples = []
    
        with open(filename) as f:
            for i, line in enumerate(f):
                nums = map(float, line.split())
                if len(nums) < 2:
                    msg = "Ray direction file {} requires two floats on \
		           each line; \
                           line {} is ill-formed".format(filename, i)
                    sys.exit(msg)
                w = np.cos(nums[0])
                u = np.sin(nums[0])*np.cos(nums[1])
                v = np.sin(nums[0])*np.sin(nums[1])
                ray_tuples.append((u, v, w))
        return ray_tuples

    def load_directions(self, ray_dir_file, num_ray_dirs):
        """ Load a numpy array with all the ray directions to be used.
    
        Parameters
        ----------
        ray_dir_file : string, filename 
            filename of file containing pairs of angles (radians) 
	    representing polar and
        azimuthal angles.  Does not have to be set.

        num_ray_dirs : integer
            Number of ray directions to use.  
        o If a filename is given and num_ray_dirs is
          nonzero and less than the number of directions in the file, 
	  the return array 
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
            ray_tuples = self.load_ray_tuples(ray_dir_file)
            if num_ray_dirs > 0:
                ray_tuples = ray_tuples[:num_ray_dirs]
            elif num_ray_dirs < 0:
	        ray_tuples = ray_tuples[-num_ray_dirs:]
        elif num_ray_dirs > 0:
            ray_tuples = self.get_rand_dirs(num_ray_dirs)
        else:
            ray_tuples = [(1.0, 0.0, 0.0),
                          (0.0, 1.0, 0.0),
                          (0.0, 0.0, 1.0)]
        return np.array(ray_tuples)

    def find_reference_vol(self):
        """ Find which volume in the geometry contains the reference point.
	
        Requirements
        ------------
        - The DAGMC python library must be loaded

	Returns
	----------
	vol_id of volume containing the reference point

        Error Condition
        ----------
        Exits if none of the volumes contain the reference point
        """

        # code snippet from dagmc.pyx:find_graveyard_inner_box
        volumes = dagmc.get_volume_list() 
	print "volume list is", volumes
	ref_vol = -1
        for v in volumes:
	    print "checking v", v, "has", self.reference_point
            if dagmc.point_in_volume(v, self.reference_point):
		print "it does!"
		ref_vol = v
        	break

        if ref_vol == -1:
            msg = "Exiting program: Reference volume not found!"
            sys.exit(msg)

        return ref_vol

    def find_slabs_for_ray(self, uvw):
        """ Find the material distances and names along a direction
        in a tagged geometry.

        Requirements
        ------------
        graveyard_vol != start_vol

        Parameters
        ----------
        uvw : float triplet
             Unit direction vector in the direction of travel

	Returns
	-------
	slab_info, with dtype:
	slab_info_type =  [('id', 'i4'), ('name','a12'),('depth', 'f4'), \
	                   ('density', 'f4'), ('areal_density', 'f4')]
	where:
        id   : vol_ids of the slabs traveled through
        name : Unique names of slabs in the same order as they are 
	       encountered, meaning in the same order as the slab_lengths
        depth : Distances traveled through materials ('slabs') encountered
        density : Densities of the materials traveled through
        areal_density : Accumulated distance*density for each slab 

        Notes
        -----
        This depends on DAGMC's ray_iterator, which calls ray_fire
        """

        surf = 0
        dist = 0
        # This geometrical hugeness, not integer hugeness
        huge = 1000000000

        ids = []
        mat_names = []
        lengths = []
        densities = []
        areal_densities = []

        prev_vol = self.reference_vol
        for (next_vol, dist, surf, xyz) in dagmc.ray_iterator(prev_vol, \
                                             self.reference_point, uvw, \
					     yield_xyz=True):
            if prev_vol in self.name_for_vol:
    	        # Get the list of vol-id's in order. 
    	        ids.append(prev_vol)
                lengths.append(dist)

	        name = self.name_for_vol[prev_vol][0]
	        dens = self.name_for_vol[prev_vol][1]
                mat_names.append(name)
	        densities.append(dens)
	        areal_densities.append(dist*dens)

            if next_vol == self.graveyard_vol or dist >= huge or surf == 0:
                break;

            # If prev_vol is the implicit complement, 
	    # this is the only line executed
            prev_vol = next_vol        

	if ( len(ids) != len(mat_names) or 
	     len(ids) != len(lengths) or 
	     len(ids) != len(densities) or
	     len(ids) != len(areal_densities) ):
	    raise Exception("find_slabs_for_ray unequal list lengths")

        slab_info = np.zeros(len(ids), dtype=self.slab_info_type)
        slab_info['id']      = ids
	slab_info['name']    = mat_names
	slab_info['depth']   = lengths
	slab_info['density'] = densities
	slab_info['areal_density'] = areal_densities

	return slab_info
        # End of find_slabs_for_ray
      
    def spatialTuples(self, slab_mat_names, slab_depthq):
        """ Implement an algorithm to create tuples from the material/distances
        lists that are in a convenient form for creating transport execution.

        Parameters
        ----------
        slab_mat_names : list of strings
            Unique material names from the geometry, in the order they are 
            encountered for a particular ray.
    
        slab_depthq : list of floats
            Areal density: Distances traveled through each material by the 
            ray multiplied by the density of that material, in the same order 
            as slab_mat_names

        Return
        ------
        spatial : list of tuples
           Tuples are (string, list of floats)
        First element of list of floats is always zero
        List of floats is length 2 at a minimum

        Note
        ----
        This functions makes no use of class attributes
        """
        spatial = []
        last_mat_name = slab_mat_names[0]
        current_slab_pts = [0.0]
        # for mat_name, slab_thick in zip(slab_mat_names, slab_lens):
        for mat_name, mat_depthq in zip(slab_mat_names, slab_depthq):
            if mat_name != last_mat_name:
                # Foreach new mat, store the name/slab points tuple...
                spatial.append((last_mat_name, current_slab_pts))
                # ... and re-initialize
                last_mat_name = mat_name
                current_slab_pts = [0.0]
            # Do always; first time through, this is all that is done
            current_slab_pts.append(current_slab_pts[-1] + mat_depthq)

        # Once more for the trailing tuple 
        spatial.append((last_mat_name, current_slab_pts))

        return spatial

    def createTransportDictionary(self, uvw):
        """ Create a dictionary containing material names and areal density 
        (length times density) for use in one-dimensional transport code.

        Parameters
        ----------
        uvw : list of 3 floats
            The direction vector of the unit ray. The final ray going through 
            the material points in this direction.

        Return
        ------
        spatial : list of tuples
            Tuples are (string, list of floats)
            First element of list of floats is always zero
            List of floats is length 2 at a minimum
       	    The tuples are calculated by spatialTuples(..) 
	        using the unique material names and the 
		areal densities *rounded to 1 digit after the decimal place*.

        info_at_ref : structured array with dtype=self.info_type
            A structured array record giving information useful when recording 
	    the transport-response data.  
	    It includes the volume id, name, and density at the reference point,
	    plus the cumulative depth and areal density to the reference point
	    slab_info_type =  [('id', 'i4'), ('name','a12'),('depth', 'f4'), 
	                   ('density', 'f4'), ('areal_density', 'f4')]

        Notes
        ----
        This is the first place we can know if the rays did not go through a 
        material.

        For the giiven ray direction, starting from the reference vol
        o get the list of materials and distances it passes through in the 
          given direction, up to the graveyard
        o get the list of materials and distances it passes through in the 
          opposite direction, up to the graveyard
        o combine the materials and distances as if that of a ray passing 
          through the reference point, which starts and ends at the graveyard
        o create the contents of the transport input file for this combined ray
        """

        # Rays fired here!
        ###################################
        # Create two lists 'a' and 'b'.  
        #    'a' is the list of slab materials and lengths in the 
	#    given direction.  
        #    'b' is the list of same items going OPPOSITE (180 degrees)
        slab_info_a = self.find_slabs_for_ray(uvw)
	num_a = len(slab_info_a['id'])

        slab_info_b = self.find_slabs_for_ray(-uvw)
	num_b = len(slab_info_b['id'])

        # The ray did not go through any material
        info_at_ref = np.zeros(1, dtype=self.info_type)
        if num_a + num_b == 0:
            return [], -1, info_at_ref

        slab_mat_names  = np.concatenate((slab_info_b['name'][::-1], \
	                  slab_info_a['name']))
        areal_densities = np.concatenate((slab_info_b['areal_density'][::-1], \
	                  slab_info_a['areal_density']))
        
	# The a and b rays start at the same point, however if the 
	# start point is in the implicit complement, the first 
	# volume encountered by the ray and it's opposite may differ.
	# Therefore it is not an error if the id's differ.
	ref_id = slab_info_b['id'][0]
	ids_match = True
	if slab_info_a['id'][0] != ref_id:
	    ids_match = False

	ref_name    = slab_info_b['name'][0]
	ref_density = slab_info_b['density'][0]
        ref_depth   = np.sum(slab_info_b['depth'])
	ref_depthq  = np.around(np.sum(slab_info_b['areal_density']), decimals=1)
	if not ids_match:
	    print "BTRN:  ids a and b (ref):", slab_info_a['id'][0], ref_id
	    print "     names a and b (ref):", slab_info_a['name'][0], ref_name
	    # print "     ray b areal_densities:", slab_info_b['areal_density']
	    # print "     Sum b, then round:", ref_depthq

        info_at_ref[0] = (num_b, ref_id, ref_name, ref_density, ref_depth, ref_depthq)
        return self.spatialTuples(slab_mat_names, areal_densities), info_at_ref
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
        '-f', action='store', dest='uwuw_filename', 
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
        '-p', '--reference_point', nargs=3, type=float,
    help='Cartesion coordinates of the reference point for measurements')
    parser.set_defaults(reference_point=[0.0,0.0,0.0])

    parser.add_argument(
        '-s', '--response', action='store', 
        help='The material in which to calculate the response.  \
          NB: the response material cross-section must exist')
    parser.set_defaults(response='WATERLIQ')

    args = parser.parse_args()

    # With no other parameters set byt the geometry file, a single randomly 
    # chosen direction will be run, using the defaults listed.
    if not args.uwuw_filenam:
        raise Exception('h5m filename not specified. [-f] not set')
    
    return args

def test_tuples(uwuw, rundir, ray_dir_file):
    """ Test function returns class instance and ray directions.
    This is all that's needed to run the transport
    """
    at = alltran(uwuw, rundir, 0)
    odt = transportresponse(at.runpath, at.datapath, \
                            at.environment, at.response)
    at.info_type = odt.info_type
    dagmc.load(uwuw)

    # Get name/vol_id information from the geometry file
    at.name_for_vol, at.graveyard_vol = \
        uwuw_preproc.get_fnames_and_densities(uwuw)

    index, ray_tuples = at.get_rays_and_reference_vol(ray_dir_file)

    return at, ray_tuples

def main():
    """Main routine for btrn.py

    Notes
    -----
    The sequence of actions:
    - Parse the command line
    - Initialize the alltran class with paths and number of directions
      o Separately initialize attributes that have defaults
    - Initialize the transportresponse class with paths and settings
    - Load the geometry file and get all its unique material (FLUKA) names
    - Determine/check the start index
    - For each ray direction
      o get geometry information necessary to create the input files
      o execute the transport with the input file
      o collect results
    """
    
    args = parse_command_line_arguments()

    # Initialize the local transport class
    at = alltran(args.uwuw_filename, args.rundir, args.num_ray_dirs)

    # Don't set these with initializer since there are default values
    at.environment     = args.environment
    at.reference_point = args.reference_point
    at.response        = args.response

    # Initialize the transportresponse class
    one_d_tool = transportresponse(at.runpath, at.datapath, \
                                   at.environment, at.response)
    # Ensure agreement on a record array dtype for information sharing
    at.info_type = one_d_tool.info_type

    # Load the DAG object for this geometry
    dagmc.load(at.uwuw_filepath)

    # Get name/vol_id information from the geometry file
    at.name_for_vol, at.graveyard_vol = \
        uwuw_preproc.get_fnames_and_densities(at.uwuw_filepath)

    index, ray_tuples = at.get_rays_and_reference_vol(args.ray_dir_file)
    if at.reference_vol in at.name_for_vol:
        print "The reference volume is material", at.name_for_vol[at.reference_vol]
    else:
        print "The reference volume is in the implicit complement"

    for uvw in ray_tuples:
        spatial_tuples, info_record =  at.createTransportDictionary(uvw)
	# print index, ".", "spatial_tuples for ", info_record
	# for item in spatial_tuples:
	#     print "\t", item
        one_d_tool.execute(spatial_tuples)

        one_d_tool.collect_results(index, uvw, info_record) 
        index = index + 1
    return 
   
if __name__ == '__main__':
    main()
