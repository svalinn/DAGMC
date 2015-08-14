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

    # int -- The vol_id of the volume containing the ref_point.
    reference_vol   = -1
    # Vol id of the graveyard, a convenience value
    graveyard_vol = -1

    slab_info_type =  [('id', 'i4'), ('name','a12'),('depth', 'f4'), \
                       ('density', 'f4'), ('areal_density', 'f4')]

    # name_for_vol :  dictionary, { int : (string , float) }
    # A mapping of vol id to  Value = (unique_material_name, density)
    name_for_vol = {}
    # dtype for the returned info record array
    info_type = []

    def __init__(self, uwuw, rundir, n):
        """
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

    def setup_rays(self, ray_dir_file): 

        rays_file  = os.path.join(os.getcwd(), ray_dir_file)
        ray_tuples = self.load_directions(rays_file, self.numdirs)

        # The default starting point is 0,0,0
        ref_vol = self.find_reference_vol()

        if ref_vol == self.graveyard_vol:
            msg = "Exiting program: The reference point is in the graveyard!"
            sys.exit(msg)
    
        if 0 <= self.numdirs:
            index = 1
        else:
            index = -self.numdirs
        return index, ray_tuples

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

    def find_reference_vol(self):
        """ Find which volume in the geometry contains the reference point.
	
        Requirements
        ------------
        - The DAGMC python library must be loaded

	Attributes
	----------
	self.reference_vol = vol_id of volume containing the reference point

        Exceptions
        ----------
        Raises exception if none of the volumes contain the reference point
        """

        # code snippet from dagmc.pyx:find_graveyard_inner_box
        volumes = dagmc.get_volume_list() 
        for v in volumes:
	    print "checking v", v, "has", self.reference_point
            if dagmc.point_in_volume(v, self.reference_point):
		print "it does!"
		self.reference_vol = v
        	break

        if self.reference_vol == -1:
            raise Exception("Reference volume not found!")

        return self.reference_vol

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
        ids : list of ints
            vol_ids of the slabs traveled through

        slab_mat_names : list of strings
            Unique names of slabs in the same order as they are encountered,
            meaning in the same order as the slab_lengths

        slab_lengths : list of floats
            Distances traveled through materials ('slabs') encountered
    
        slab_densities : list of floats
            Densities of the materials traveled through

        slab_areal_densities : list of floats
            Accumulated distance*density for each slab travelled through
	

        Notes
        -----
        This depends on DAGMC's ray_iterator, which calls ray_fire
        """

        surf = 0
        dist = 0
        # This geometrical hugeness, not integer hugeness
        huge = 1000000000

        ids = []
        slab_mat_names = []
        lengths = []
        densities = []
        slab_areal_densities = []

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
                slab_mat_names.append(name)
	        densities.append(dens)
	        slab_areal_densities.append(dist*dens)

            if next_vol == self.graveyard_vol or dist >= huge or surf == 0:
                break;

            # If prev_vol is the implicit complement, 
	    # this is the only line executed
            prev_vol = next_vol        

	"""
	slab_info_type =  [('id', 'i4'), ('name','a12'),('depth', 'f4'), \
	                   ('density', 'f4'), ('areal_density', 'f4')]
	"""
	if ( len(ids) != len(slab_mat_names) or 
	     len(ids) != len(lengths) or 
	     len(ids) != len(densities) or
	     len(ids) != len(slab_areal_densities) ):
	    raise Exception("find_slabs_for_ray unequal list lengths")

	# info_len = len(ids)
        slab_info = np.zeros(len(ids), dtype=self.slab_info_type)
	print "ids", ids, "slab_mat_names", slab_mat_names, \
	      "lengths", lengths, "densities", densities, \
	      "slab_areal_densities", slab_areal_densities
        slab_info['id'] = ids
	slab_info['name'] = slab_mat_names
	slab_info['depth'] = lengths
	slab_info['density'] = densities
	slab_info['areal_density'] = slab_areal_densities
        # return ids, slab_mat_names, lengths, densities, slab_areal_densities
	return slab_info
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
        for mat_name, slab_depthq in zip(slab_mat_names, slab_depthq):
            if mat_name != last_mat_name:
                # Foreach new mat, store the name/slab points tuple...
                spatial.append((last_mat_name, current_slab_pts))
                # ... and re-initialize
                last_mat_name = mat_name
                current_slab_pts = [0.0]
            # Do always; first time through, this is all that is done
            current_slab_pts.append(current_slab_pts[-1] + slab_depthq)

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


        start_vol : int
            Vol id of the volume containing the reference point.

        Return
        ------
        spatial : list of tuples
            Tuples are (string, list of floats)
            First element of list of floats is always zero
            List of floats is length 2 at a minimum
       	    The tuples are calculated by spatialTuples(..) 
	        using the unique material names and the 
		areal densities *rounded to 1 digit after the decimal place*.

        ref_slab_depth : int
            Number of slabs the ray goes through to get to the reference point.
    
        info_at_ref : structured array with dtype=self.info_type
            A structured array record giving information useful when recording 
	    the transport-response data.  
	    It includes the volume id, name, and density at the reference point,
	    plus the cumulative depth and areal density to the reference point

        Note
        ----
        This is the first place we can know if the rays did not go through a 
        material.
    - For each ray direction, starting from the reference vol
      o get the list of materials and distances it passes through in the 
        given direction, up to the graveyard
      o get the list of materials and distances it passes through in the 
        opposite direction, up to the graveyard
      o combine the materials and distances as if that of a ray passing 
        through the reference point,
        which starts and ends at the graveyard
      o create the contents of the transport input file for this combined ray
        """

        # Rays fired here!
        ###################################
        # Create two lists 'a' and 'b'.  
        #    'a' is the list of slab materials and lengths in the 
	#    given direction.  
        #    'b' is the list of same items going OPPOSITE (180 degrees)
	"""
	slab_info_type =  [('id', 'i4'), ('name','a12'),('depth', 'f4'), \
	                   ('density', 'f4'), ('areal_density', 'f4')]
	"""
        # slab_ids_a, slab_mat_names_a, slab_lens_a, slab_densities_a, areal_densities_a = \
	print "Slab A"
        slab_info_a = self.find_slabs_for_ray(uvw)
	num_a = len(slab_info_a['id'])

        # slab_ids_b, slab_mat_names_b, slab_lens_b, slab_densities_b, areal_densities_b = \
	print "Slab B"
        slab_info_b = self.find_slabs_for_ray(-uvw)
	num_b = len(slab_info_b['id'])

        print "a id", slab_info_a['id'], "a name", slab_info_a['name'], \
	      "a depth", slab_info_a['depth'], "a density", \
	      slab_info_a['density'], "a areal_density", \
	      slab_info_a['areal_density']

        print "b id", slab_info_b['id'], "b name", slab_info_b['name'], \
	      "b depth", slab_info_b['depth'], "b density", \
	      slab_info_b['density'], "b areal_density", \
	      slab_info_b['areal_density']

        # The ray did not go through any material
        info_at_ref = np.zeros(1, dtype=self.info_type)
        if num_a + num_b == 0:
            return [], -1, info_at_ref

        num_slabs_to_ref = num_b

        slab_mat_names  = np.concatenate((slab_info_b['name'][::-1], \
	                  slab_info_a['name']))
        areal_densities = np.concatenate((slab_info_b['areal_density'][::-1], \
	                  slab_info_a['areal_density']))
        
	ref_id = slab_info_b['id'][0]
	if slab_info_a['id'][0] != ref_id:
	    raise Exception("Ids at reference point differ!")

	ref_name = slab_info_b['name'][0]
	if slab_info_a['name'][0] != ref_name:
	    raise Exception("Names at reference point differ!")

        # ref_density = slab_densities[slabs_to_ref]
	ref_density = slab_info_b['density'][0]
	if slab_info_a['density'][0] != ref_density:
	    raise Exception("Densitis at reference point differ!")

        # ref_depth   = np.sum(slab_lens[:slabs_to_ref])
        ref_depth = np.sum(slab_info_b['depth'])
        print "sum slab_info_b", ref_depth

        # ref_depthq  = np.sum(np.around(areal_densities[:slabs_to_ref],1))
        ref_depthq = np.sum(np.around(slab_info_b['areal_density'],1))
        print "areal_densities to ref", ref_depthq
    
        info_at_ref[0] = \
	    (slabs_to_ref, ref_id, ref_name, ref_density, ref_depth, ref_depthq)
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
    if not args.uwuw_filename:
        raise Exception('h5m filename not specified. [-f] not set')
    
    return args

def test(uwuw, rundir):
    at = alltran(uwuw, rundir, 0)
    odt = transportresponse(at.runpath, at.datapath, \
                            at.environment, at.response)
    at.info_type = odt.info_type
    dagmc.load(uwuw)

    # Get name/vol_id information from the geometry file
    at.name_for_vol, at.graveyard_vol = \
        uwuw_preproc.get_fnames_and_densities(uwuw)

    # index = at.setup(args.ray_dir_file)

    return at

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
    at.environment = args.environment
    at.ref_point   = args.reference_point
    at.response    = args.response

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

    index, ray_tuples = at.setup_rays(args.ray_dir_file)

    for uvw in ray_tuples:
        spatial_tuples, info_record =  at.createTransportDictionary(uvw)
        one_d_tool.execute(spatial_tuples)
	
	"""
	for key in name_density_dict:
	    print "id, name, density", key, name_density_dict[key][0], 
	                               name_density_dict[key][1]
	"""
        one_d_tool.collect_results(index, uvw, info_record) 
        index = index + 1
    return 
   
if __name__ == '__main__':
    main()
