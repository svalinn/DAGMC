import argparse
import numpy as np
try:
    from scipy.spatial import ConvexHull
except ImportError:
    raise ImportError("scipy.spatial.ConvexHull could not be imported.")


class resultrec(object):

    # Number of columns devoted to direction and other
    # non-response data
    num_meta = 6
    num_let_groups  = 700
    num_flux_groups = 100
    neutflux_names = ['forward', 'backward', 'total']

    # Format for meta data
    meta_type=[('UVW', 'f4', 3), ('density','f4'), ('depth','f4'), ('depthq','f4')]
    #meta_type=[('u','f4'), ('v','f4'), ('w','f4'), 
    #           ('density','f4'), ('depth','f4'), ('depthq','f4')]

    # Initialized by c'tor
    nrows = -1
    response_type    = []
    response_data    = []
    tr_meta  = np.zeros([])
    tr_all   = np.zeros([])
    tr_let_group  = np.zeros([])
    tr_flux_flux  = np.zeros([])
    vertices = np.zeros([])

    # set by process_dose_header()
    num_species = 0
    particle_names = []


    def __init__(self, infile):
        """ Constructor for the resultrec class
        """
	row_data, meta_dose_header, group_header = self.get_data(infile)
        dose_header = meta_dose_header[self.num_meta:]
        self.process_dose_header(dose_header)

        self.nrows = row_data.shape[0]
	self.set_tr_meta(row_data)

        self.tr_let  = set_group_struct(0, num_let_groups, group_header)
	self.tr_flux = set_group_struct(0, num_flux_groups, group_header[2*num_let_groups:])

	# Set up for VTK file creation
	# ToDo: this is embedded in tr_meta, it should be deleted
	ray_dirs      = np.array(row_data)
	self.vertices = ray_dirs[:,:3]

        # remove meta_columns from alldata
        response_data = row_data[:,self.num_meta:]
        self.set_tr_all(response_data)
        return
	

    def get_data(self, infile):
        """ Read the lines of a database file, stripping 
        off the header and statistics lines.
    
        Parameters
        ----------
        infile : file with a header, 2 statistics rows, and otherwise
                 row by column floats

	Throws
	______
	Exception if sum of columns in the headers does not match the number
	of data columns.

        Returns
        -------
        allrows -- a 2D array of floats containing all the meta-data
                   AND responses.
        meta_dose_header -- list of meta-data, dose and doseq headers,
                            all of which are ascii
        group_header     -- list of header fields which are numeric,
                            e.g. energy boundaries for let data
    			This list has many repeated sequences
    			Note that meta_dose_header + group_header 
    			a) is the entire header and
    			b) has a length related to the row data
    
        Notes
        -----
    	- len(allrows) = len(meta_dose_header) + len(group_header)

        - File is the result of calling encode.py on an HZETRN dataset
    
        end of get_data(infile)
        """ 
        allrows = np.array([])
        with open(infile) as fp:
            for index, line in enumerate(fp):
    	        i = index + 1
    	        if i == 1:
    	            header = line.split()
    	        else:
    	            try: 
    	                row = map(float, line.split())
    	                if len(allrows) == 0:
    	                    allrows = row
    	                else:
    	                    allrows = np.vstack((allrows, row))
		    # Statistics row will not map to floats
    	            except ValueError:
    		        pass
    	           
        # Provide a subset of the header that is just the ascii column headers.
        # The rest of the headers can be converted to energy group floats
        max_non_numeric_col = -1    
	meta_dose_header = []
        group_header     = []
        for field in header:
            try:
    	        num = float(field)
    	    # it doesn't convert, must be ascii
    	    except ValueError:
    	        meta_dose_header.append(field)
    	    else:
    	        group_header.append(num)

        if len(meta_dose_header) + len(group_header) != allrows.shape[1]:
            raise Exception('The header did not split nicely into non-numeric and numeric sections')
        
	return allrows, meta_dose_header, group_header
    
    def process_dose_header(self, dose_header):
        """ Process the dose_header attribute to get information about
	    particle names and number

        Parameters
	----------
	    dose_header - array containing the header columns for dose and doseq

        Attributes
        -------
	self.num_species i--  gotten from the dose_header
	self.particle_names -- list of names that matches the input data

	Raises
	------
	Exception if the calculated number of species does not equal the number of names
	"""
        num = int((len(dose_header) - 2) / 2)

        # To fill in the PARTICLE field, extract name list depending
        # on the number of species, i.e. the radiation environment
        if 6 == num:
           names = [x.split('_')[1] for x in dose_header[2:2+num]] 
        else:
           names = ['particle_' + str(i) for i in range(num)]

	if len(names) != num:
	    raise Exception('The number of particle names does not match the number of species!')
	
	self.num_species    = num
	self.particle_names = names
        return 

    def set_tr_meta(self, data):
        """ Set a structured array consisting of the first num_meta columns
        of the passed in array.
    
        Parameters
        ----------
        data : numpy array with at least num_meta columns
            a two-dimensional array of floats 
    
        Attributes
        -------
	self.tr_meta -  A structured array with the same number of tuples as 
	                the rows in the input data array
    
        """ 
        self.tr_meta = np.zeros(self.nrows, dtype=self.meta_type) 
        for i in range(self.nrows):
            # self.tr_meta[i] = tuple(data[i,:self.num_meta])
            self.tr_meta[i] = tuple( (data[i,0:2]), data[3:self.num_meta] )
        return 
    
    # ToDo - make this work fo r let and flux groups
    def set_group_struct(start, group_size, header):
        """ Get a structured array for the energy groups. 

        Parameters
        ----------
        start : int 
            The first column of the let groups in the array

        group_size : int
            700 for let, 100 for flux

        header : array of floats 

        Note
        ----
        These were read in as ascii values from a file and have
        previously been converted to floats
        """
        end = start + group_size
        group_header = header[start:end]
        tr_group = np.zeros(group_size, \
               dtype = [('id', 'i4'), ('energy', 'f4')])

        for i in range(group_size):
            tr_group[i] = (i, group_header[i])
        return tr_group
    
    def init_tr_all(self, response_data):
        """ Put just the response data into a recarray

	Parameters
	----------
        response_data : np array, shape = (nrows, repsonse_cols)
         
	Attributes
	----------
	self.response_type
	self.response_data
	self.tr_all

	"""

        self.response_type = [('RESPONSE', 'a12'), ('PARTICLE', 'a12'), \
	                      ('GROUP', 'i4'), ('VALUES', 'f4', self.nrows)]

	response_cols = self.response_data.shape[1]

	# Create a record array for 
	# a) all the incoming data cols PLUS
	# b) the flux for each particle (i.e. summed over energy groups) PLUS 
	# c) the total neutron backward flux
        total_cols = response_cols + self.num_species + 1
        self.tr_all = np.zeros(response_cols, dtype=self.response_type)

	# Make one record with placeholders
        for col in range(response_cols):
            self.tr_all[col] = ('',  'All', -1, self.response_data[:, col])
        return

    def fill_dose_tuples(self, start, response_name):
	""" Set the non-value portion of the tuples for dose and doseq

	Parameters
	----------
	start : integer
	    The next column whose recarray tuple needs to be completed

        response_name : string
	    The type of response, dose or doseq

	Returns
	-------
        The index after the last column whose tuple was set
	"""
        for particle, d in zip(self.particle_names, range(start, start + self.num_species)):
	    self.tr_all[d]['RESPONSE'] = response_name
	    self.tr_all[d]['PARTICLE'] = particle
        return start + self.num_species

    def fill_let_tuples(self, start, response_name):
	""" Set the non-value portion of the tuples for let

	Parameters
	----------
	start : integer
	    The next column whose recarray tuple needs to be completed

        response_name : string
	    The type of response, dif_let or int_let

	Returns
	-------
        The index after the last column whose tuple was set
	"""
	for group_number, g in zip(range(self.num_let_groups), range(start, start + self.num_let_groups)):
	    self.tr_all[g]['RESPONSE'] = response_name
	    self.tr_all[g]['PARTICLE'] = 'All'
	    self.tr_all[g]['GROUP'] = group_number 
	return start + num_let_groups

    def fill_flux_tuples(self, start, response_name):
	""" Set the non-value portion of the tuples for flux

	Parameters
	----------
	start : integer
	    The next column whose recarray tuple needs to be completed

        response_name : string
	    The type of response, flux or intflux

	Returns
	-------
        The index after the last column whose tuple was set
	"""
        n_start = start

	for particle in self.particle_names:
	    # DUP1 - almost identical to fill_let_tuples
	    for group_number, g in zip(range(self.num_flux_groups), range(n_start, n_start + self.num_flux_groups)):
	        self.tr_all[g]['RESPONSE'] = response_name
	        self.tr_all[g]['PARTICLE'] = particle
	        self.tr_all[g]['GROUP']    = group_number 
	n_start += self.num_flux_groups
	return start + self.num_species * self.num_flux_groups

    def fill_neutflux_tuples(self, start):
	""" Set the non-value portion of the tuples for neuflux

	Parameters
	----------
	start : integer
	    The next column whose recarray tuple needs to be completed

	Returns
	-------
            The index after the last column whose tuple was set

	Notes
	-----
	    This also sums all columns for backwards neutron flux and sets
	    an additional record tuple with it
	"""
        total_backward_col = np.zeros([])
	# Some convenience parameters; note len(neutflux_names) == 3
        num_neutflux = len(self.neutflux_names)
        n_start = start

	for name in self.neutflux_names:
	    print n_start
	    # DUP2 - almost identical to fill_let_tuples
	    for group_number, g in zip(range(self.num_flux_groups), range(n_start, n_start + self.num_flux_groups)):
	        self.tr_all[g]['RESPONSE'] = 'neutflux'
	        self.tr_all[g]['PARTICLE'] = name
	        self.tr_all[g]['GROUP']    = group_number
		# sum the total neutron columns over group
	        if name == 'backward':
		    total_backward_col += self.tr_all[g]['VALUES']
	    n_start += self.num_flux_groups

	last = start + num_neutflux_names*self.num_flux_groups + 1
        self.tr_all[last] = ('totalbackneut', 'neutron', -1, total_backward_col)
        return last + 1

    def set_tr_all(self, data):
        """ Set all the tuples in tr_all.  

	Parameters
	----------
	data : two dimentional numpy array
	       columns are particles or energy groups
	       rows are different ray directions/meta data

	Attributes
	----------
	self.tr_all
	"""
        self.init_tr_all(data)
        ######################### Dose #################################
        # Fill in the dose and doseq RESPONSE and PARTICLE headers
	self.tr_all[0]['RESPONSE'] = 'Dose'
	self.tr_all[1]['RESPONSE'] = 'Doseq'

	start = self.fill_dose_tuples(2, 'Dose')
	start = self.fill_dose_tuples(start, 'Doseq')
        ######################### LET  #################################
	start = self.fill_let_tuples(start, 'dif_let')
	start = self.fill_let_tuples(start, 'int_let')
        ######################### Flux #################################
	# reference for summing over groups
	flux_start = start
        start = self.fill_flux_tuples(start, 'flux')
        start = self.fill_flux_tuples(start, 'intflux')
	start = fill_neutflux_tuples(start) 

	# Add in the summed flux cols.  Need to go back to previous columns
	flux_p = flux_start
	for particle, i in (self.particle_names, range(start, start + self.num_species)):
	    self.tr_all[i]['RESPONSE'] = 'totalflux'
	    self.tr_all[i]['PARTICLE'] = particle
	    self.tr_all[i]['GROUP']    = -1 
	    for flux_j in range(flux_p, flux_p + self.num_flux_groups):
	        value_col += self.tr_all[flux_j]['VALUES']
	    self.tr_all[i]['VALUES'] = value_col
	    flux_p += self.num_flux_groups        

        last = start + self.num_species
	print last
	return
    

    def init_mesh():
        """Create an iMesh for the directions of this problem.
	Notes
	-----
	    No tagging yet

	Returns
	-------
            The mesh, and other info necessary to tag it
	"""
	# ToDo - get the vertices from self.tr_meta
	# Double check this is what is needed
	# WAS self.vertices
	vtcs = self.tr_meta['UVW']
        # Triangulate 
        hull       = ConvexHull(vtcs)
        indices    = hull.simplices
        facets     = vtcs[indices]
        num_facets = facets.shape[0]

        # Create the mesh, and a data tag on it
        msph = iMesh.Mesh()
	return msph, num_facets, indices


    def tag_mesh(mesh, num_facets, indices, tuple_set):
        """ Create the writable VTK object from the tuple set and the directions

        """
        tag_map = {}

	# Redo for-loop from original tag_mesh for recarrays
	for stuff in tuple_set:
	    data = stuff['VALUES']
	    di = data[indices]
	    data_name = stuff['PARTICLE'] + '_' + stuff['RESPONSE'] 
	    # NOTE: The default set will not have group data
	    if stuff['GROUP'] != -1:
		# Wrong, use recarray to get the group value.  Also, make sure
		# the correct group is identified (let or flux)
	        data_name += '_' + str(stuff['GROUP'])
	    # Associate the hashable Tag array object with the ascii data_name
            tag_map[data_name] = mesh.createTag(data_name, 1, float)
	    """
            for col in range(start, end): 
            # get the desired column of data
	    ######### GET FROM RECARRAY ###########
	    data = vtxd[:,col]
	    # Reorder the data column in terms of vertices
            di = data[indices]
            # creat a header-field-based name for the current column
            data_name = header[col]+"_{}".format(col)
	    # Associate the hashable Tag array object with the ascii data_name
            tag_map[data_name] = mesh.createTag(data_name, 1, float)
            """
            # Use indexing to get matching data at each vertex
            for i in range(num_facets):
                facet        = facets[i]
    	        data_at_vtcs = di[i]
    	        # Create an entity handle for each point in the facet
    	        verts = mesh.createVtx(facet)
    	        # Tag each entity handle (lhs key) with  corresponding data value (rhs)
    	        tag_map[data_name][verts[0]] = data_at_vtcs[0]
    	        tag_map[data_name][verts[1]] = data_at_vtcs[1]
    	        tag_map[data_name][verts[2]] = data_at_vtcs[2]
    
    	        # Tell the mesh about this triangular facet
    	        tri, stat = msph.createEntArr(iMesh.Topology.triangle, verts)
	
        # All that's left is to mesh.save(outfile)
        return mesh

    def get_default_tuples():
        """ Prepare default inputs for meshing, nameing and writing to VTK

	Notes
	-----
        See results_vis.py:tag_mesh

        HZETRN Default output on sphere and surface 
        1 TotalDose, 
        1 TotalDoseEq, 
        N Dose(per particle), 
        N DoseEq(per particle), 
        N TotalFlux(per particle), 
        1 TotalFluxBackNeutrons 
        This is 2+N*3+1 plots.  For SPE N is 6 so we have 21, for GCR N is 59 so we have 180.
        """
 
        default_tuples = []
	for item in tr_all:
	   if item['RESPONSE'] == 'Dose' or item['RESPONSE'] == 'Doseq' or \
	      item['RESPONSE'] == 'totalflux' or item['RESPONSE'] == 'totalbackneut':
	       default_tuples.append(item)
	
	return default_tuples

def parse_arguments():
    """ Argument parsing
    """
    
    min_col = self.num_meta
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', action='store', dest='infile',
	help='The name of the .dat file with one line per direction.')


    parser.add_argument(
        '-r', action='store', dest='response',
	help='Response, from dose, doseq, let, flux, int_flux, or neutron_flux.')

    parser.add_argument(
        '-p', action='store', dest='particle',
	help="Particle: 'all', name-list (dose, doseq); 'int', 'dif' (let); name-list (flux, int_flux); 'forward', 'backward', 'total' (neutron_flux).")

    parser.add_argument(
        '-g', action='store', dest='group',
	help='Group: 1-700 (let); 1-100 (flux, int_flux, neutron_flux). Default = entire range. Not used for (dose, doseq).')

    # parser.add_argument(
    #     '-c', action='store', dest='data_column', 
    #	help='The column, column range, or "a" for all columns.  Must be {} or greater'.format(num_non_data))

    parser.add_argument(
        '-o', action='store', dest='outfile',
	help='The name of the .vtk file, with extension.')

    args = parser.parse_args()
    if not args.infile:
        raise Exception('Input file not specified. [-f] not set.')

    if not args.outfile:
        args.outfile = 'out.dat'
    
    return args



# ToDo Delete
def get_meta_from_file(infile):
    """ Get a structured array from data in a file.

    Parameters
    ----------
    filename : file with a header, 2 statistics rows, and otherwise
               row by column floats

    Reference
    ---------
    See get_meta_struct

    """ 
    dm_header, group_header, row_data = get_data(infile)
    return get_meta_struct(row_data[:,:4])

"""
"""
def list_particle_names(infile):
    dm_header, group_header, row_data = get_data(infile)
    return get_particle_names(dm_header[num_meta:])
    
def write_dose(infile, name_list):
    dm_header, group_header, row_data = get_data(args.infile)
        
def data_structs_from_file(infile):
    global response_type
    dm_header, group_header, row_data = get_data(infile)
    nrows = row_data.shape[0]
    if not response_type:
        set_response_type(nrows)
        print "response_type", response_type

    tr_meta = get_meta_struct(row_data[:,:num_meta])
    return get_data_structs(dm_header, group_header, row_data)


def main():
    args = parse_arguments()

    db = resultrec(args.infile)
    def_set = db.get_default_tuples()
    mesh, num_facets, indices = db.init_mesh()
    self.tagged_mesh = db.tag_mesh(mesh, num_facets, indices, def_set)
    tagged_mesh.save(args.outfile)

    ####################################################################
    """
    with open('myfile.dat', 'w') as f:
	for row in tr_meta:
	    for el in row:
	        f.write(str(el) + ' ')
	    f.write('\n')
    """
if __name__ == '__main__':
    main()
