import argparse
import numpy as np

class resultrec(object):

    # Number of columns devoted to direction and other
    # non-response data
    num_meta = 6
    num_let_groups  = 700
    num_flux_groups = 100
    neutflux_names = ['forward', 'backward', 'total']

    # Format for meta data
    meta_type=[('u','f4'), ('v','f4'), ('w','f4'), 
               ('density','f4'), ('depth','f4'), ('depthq','f4')]

    # Initialized by c'tor
    nrows = -1
    dose_header      = []
    group_header     = []
    response_type    = []
    response_data    = []
    tr_meta  = np.zeros([])
    tr_all   = np.zeros([])
    tr_group = np.zeros([])

    # set by get_particle_names()
    num_species = 0
    particle_names = []

    def __init__(self, infile):
        # set group_header
	row_data, meta_dose_header = self.get_data(infile)
        self.dose_header = meta_dose_header[self.num_meta:]
        self.process_dose_header()

        self.nrows = row_data.shape[0]
	self.set_tr_meta(row_data)

        # remove meta_columns from arrays
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
    
	Attributes
	----------
	self.meta_dose_header, self.group_header

        group_header     -- list of header fields which are numeric,
                            e.g. energy boundaries for let data
    			This list has many repeated sequences
    			Note that meta_dose_header + group_header 
    			a) is the entire header and
    			b) has a length related to the row data

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
    	            except ValueError:
    		        pass
    	           
        # Provide a subset of the header that is just the ascii column headers.
        # The rest of the headers can be converted to energy group floats
        max_non_numeric_col = -1    
	meta_dose_header = []
        for field in header:
            try:
    	        num = float(field)
    	    # it doesn't convert, must be ascii
    	    except ValueError:
    	        meta_dose_header.append(field)
    	    else:
    	        self.group_header.append(num)

        if len(meta_dose_header) + len(self.group_header) != allrows.shape[1]:
            raise Exception('The header did not split nicely into non-numeric and numeric sections')
        
	return allrows, meta_dose_header
    
    def process_dose_header(self):
        """ Process the dose_header attribute to get information about
	    particle names and number

        Attributes
        -------
	self.num_species i--  gotten from the dose_header
	self.particle_names -- list of names that matches the input data

	Raises
	------
	Exception if the calculated number of species does not equal the number of names
	"""
        num = int((len(self.dose_header) - 2) / 2)

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
            self.tr_meta[i] = tuple(data[i,:self.num_meta])
        return 
    
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
    

    def get_default_structs
        """
        HZETRN Default output on sphere and surface 
        1 TotalDose, 
        1 TotalDoseEq, 
        N Dose(per particle), 
        N DoseEq(per particle), 
        N TotalFlux(per particle), 
        1 TotalFluxBackNeutrons 
        This is 2+N*3+1 plots.  For SPE N is 6 so we have 21, for GCR N is 59 so we have 180.
        """
        # Construct the by-column record array
        ple_names = get_particle_names(dose_header)
        self.num_species = len(ple_names)

        ######################### Dose  #################################
        dose_header = self.dose_meta_header[self.num_meta:]
        num_dose_response = len(dose_header) 

        # Construct the dose struct
        self.tr_dose = np.zeros(num_dose_response, dtype=self.response_type)
	i = 0
	self.tr_all[i] = ('Dose', 'All', -1, self.tr_all[i]['VALUES'])
	i = 1
	self.tr_all[i] = ('Doseq', 'All', -1, self.tr_all[i]['VALUES'])
	start = i+1
        for i in range(start,start + self.num_species):
	    self.tr_all[i] = ('Dose', ple_names[i-start], -1, self.tr_all[i]['VALUES'])
	start = i+1
        for i in range(start,start + self.num_species):
	    self.tr_all[i] = ('Doseq', ple_names[i-start], -1, self.tr_all[i]['VALUES'])
	# may not need
        end_dose = i+1
	for j in range(0,end_dose):
            self.tr_dose[j] = self.tr_all[j]
	return

def get_data_structs(dose_meta_header, group_header, rows):
    """ Translate all the data into a set of record arrays.

    Parameters
    ----------
    dose_meta_header : The header array (1xn) consisting of 
        meta_data      u,v,w,density,depth, depthq, not used in this function
	dose header    ascii names of total and particle dose and doseq

    group_header : The rest of the header, i.e. the energy boundaries (float)
        matching each column of data.  The energy boundary value groups are
	repeated for particles and responses

    rows : The raw data from the data file, both meta- and response

    Returns
    -------
    A set of record arrays, all with the same dtype:
    [('RESPONSE', 'a12'), ('PARTICLE', 'a12'),  ('GROUP', 'i4'), ('VALUES', 'f4', nrows)])
    return tr_dose, tr_let, tr_flux, tr_intflux, tr_neutflux, tr_flux_part, tr_neutflux_part

    Note that the values section of the dtype is an array (the column)

    These records contain all the response data:
    tr_dose*     RESPONSE in 'dose', 'doseq' PARTICLE in ['all',particle_names]         GROUP = -1
    tr_let       RESPONSE = 'let'  	     PARTICLE in 'int', 'dif'                   GROUP in 1-700
    tr_flux*     RESPONSE = 'flux'           PARTICLE in particle_names                 GROUP in 1-100 foreach pname
    tr_intflux   RESPONSE = 'intflux'        PARTICLE in particle_names                 GROUP in 1-100 foreach pname
    tr_neutflux  RESPONSE = 'neutflux'       PARTICLE in 'forward', 'backward', 'total' GROUP in 1-100 foreach pname
    Two more records are sums over flux groups:
    tr_flux_part*      RESPONSE = 'flux_part'      PARTICLE in particle_names                 GROUP = -1
    tr_neutflux_part*  RESPONSE = 'neutflux_part'  PARTICLE in 'forward', 'backward', 'total' GROUP = -1

    """ 
    global response_type
    # remove meta_columns from arrays
    data        = rows[:,num_meta:]
    dose_header = dose_meta_header[num_meta:]

    nrows = data.shape[0]
    ncols = data.shape[1]

    # Construct the by-column record array
    ple_names = get_particle_names(dose_header)
    num_species = len(ple_names)

    ######################### Dose  #################################
    num_dose_response = len(dose_header) 

    # Construct the dose struct
    tr_dose = np.zeros(num_dose_response, dtype=response_type)
    for col in range(len(dose_header)):
        if col == 0:
            tr_dose[col] = ('Dose',  'All', -1, data[:, col])
	elif col == 1:
            tr_dose[col] = ('Doseq', 'All', -1, data[:, col])
	elif col < num_species + 2:
            tr_dose[col] = ('Dose',  ple_names[col-2], -1, data[:, col])
        else:
            tr_dose[col] = ('Doseq', ple_names[col-2-num_species], -1, data[:, col])

    ######################### LET  #################################
    num_let_groups = 700
    tr_let   = np.zeros(num_let_groups*2, dtype=response_type)
    let_data = data[:,num_dose_response:]
    for i in range(num_let_groups):
       tr_let[i] = ('dif', 'All', i, let_data[:,i]) 
       tr_let[i+num_let_groups] = ('int', 'All', i, \
                                    let_data[:,i+num_let_groups]) 
    ######################### Flux #################################
    num_flux_groups = 100
    num_flux_cols = num_flux_groups*num_species

    # Chop off the let columns... 
    flux_data     = let_data[:,num_let_groups*2:]
    # ...then the flux-by-particle column
    neutflux_data = flux_data[:,num_flux_cols*2:]

    # Create separate record arrays because intflux is not in the default set
    num_neutflux_ples = 3
    tr_flux        = np.zeros(num_flux_cols, dtype=response_type)
    tr_flux_part   = np.zeros(num_species, dtype=response_type)
    tr_intflux     = np.zeros(num_flux_cols, dtype=response_type)
    tr_neutflux    = np.zeros(num_neutflux_ples*num_flux_groups, \
                              dtype=response_type)
    tr_neutflux_part = np.zeros(num_neutflux_ples, dtype=response_type)
    neutflux_names = ['forward', 'backward', 'total']
    for j in range(num_species):
	# The group columns for this particle start here
	pname_j = ple_names[j]
        part_j_start = j*num_flux_groups

	part_flux_col   = np.sum(flux_data[:,part_j_start:part_j_start + num_flux_groups])
        tr_flux_part[j] = ('flux_part', pname_j, -1, part_flux_col)

	if num_neutflux_ples > j:
            part_neutflux_col = np.sum(neutflux_data[:,part_j_start:part_j_start + num_flux_groups])
	    tr_neutflux_part[j] = ('neutflux_part', neutflux_names[j], \
	                            -1, part_neutflux_col) 
	
        for group_i in range(num_flux_groups):
	    flux_col_k = group_i + part_j_start

	    tr_flux[flux_col_k] = ('flux', pname_j, group_i, \
	                            flux_data[:,flux_col_k])

	    tr_intflux[flux_col_k] = ('intflux', pname_j, group_i, \
	                               flux_data[:,flux_col_k + num_flux_cols])
	    if num_neutflux_ples > j:
	         tr_neutflux[flux_col_k] = ('neutflux', neutflux_names[j], \
	                                  group_i, neutflux_data[:,flux_col_k])
	
    return tr_dose, tr_let, tr_flux, tr_intflux, tr_neutflux, tr_flux_part, tr_neutflux_part

def check_neutron_flux(tr_neut_flux):
    forward  = tr_neut_flux['VALUES'][tr_neut_flux['PARTICLE'] == 'forward']
    backward = tr_neut_flux['VALUES'][tr_neut_flux['PARTICLE'] == 'backward']
    total    = tr_neut_flux['VALUES'][tr_neut_flux['PARTICLE'] == 'total']
    check    = np.add(forward, backward)

    groups = total.shape[0]
    cols   = total.shape[1]
    resid  = np.zeros((groups,cols))
    presid = np.zeros((groups,cols))
    max_resid = 0.0
    for i in range(cols):
       for j in range(groups):
           if total[j,i] != check[j,i]:
	      resid[j,i] = (total[j,i] - check[j,i])/total[j,i]
	      presid[j,i] = 100*np.abs(resid[j,i])
	      if presid[j,i] > max_resid:
	          max_resid = presid[j,i]
		  print i,j, "max_resid", max_resid
	      if presid[j,i] > .1:
                  print i, j, ':', forward[j,i], backward[j,i], check[j,i], total[j,i], presid[j,i]
    print max_resid
    return presid

def parse_arguments():
    """ Argument parsing
    """
    
    min_col = num_meta
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


def get_group_struct(start, group_size, header):
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
    # global response_type 
    args = parse_arguments()

    db = resultrec(args.infile)

    ########################
    # Get the header with the meta-data and dose,
    # the complete header, and the floating point row data
    # dm_header, group_header, row_data = get_data(args.infile)
    # nrows = row_data.shape[0]
    # Format for all response data
    # set_response_type(nrows)
    ##########################################################
    ##########################################################
    tr_dose, tr_let, tr_flux, tr_intflux, tr_neut_flux, tr_flux_part, tr_neutflux_part \
         = get_data_structs(dm_header, group_header, row_data)

    tr_meta = get_meta_struct(row_data[:,:num_meta])
    #dose_header = dm_header[4:]
    #dose_data = row_data[:,4:len(dm_header)]
    #tr_dose = get_dose_struct(dose_header, dose_data)
    ##########################################################
    # tr_let_group, tr_flux_group = get_groups_from_file_data(group_header)
    ##########################################################
    group_data = row_data[:,len(dm_header):]
    # tr_get_let_structs(group_header, group_data)

    ####################################################################
    # Do all data cols at once
    tr_all_cols, tr_all_resp = get_data_struct(all_header[num_meta:], 
                                               row_data[:,num_meta:], 700, 100)
    with open('myfile.dat', 'w') as f:
	for row in tr_meta:
	    for el in row:
	        f.write(str(el) + ' ')
	    f.write('\n')

if __name__ == '__main__':
    main()
