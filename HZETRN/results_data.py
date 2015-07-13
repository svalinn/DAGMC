import argparse
import numpy as np

# Number of columns devoted to direction and other
# non-response data
num_meta = 6
# Format for all response data
# response_type = []

# Format for meta data
meta_type=[('u','f4'), ('v','f4'), ('w','f4'), 
           ('density','f4'), ('depth','f4'), ('depthq','f4')]

def get_data(infile):
    """ Read the lines of a database file, stripping 
    off the header and statistics lines.

    Parameters
    ----------
    infile : file with a header, 2 statistics rows, and otherwise
             row by column floats

    Returns
    -------
    The data is broken up into three  sections
    meta_dose_header -- list of meta-data, dose and doseq headers,
                        all of which are ascii
    group_header     -- list of header fields which are numeric,
                        e.g. energy boundaries for let data
			This list has many repeated sequences
			Note that meta_dose_header + group_header 
			a) is the entire header and
			b) has a length related to the row data
    allrows -- a 2D array of floats containing all the meta-data
               AND responses.
	       len(allrows) = len(meta_dose_header) + len(group_header)

    Notes
    -----
    File is the result of calling encode.py on an HZETRN dataset

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
	           
    dose_meta_header = []
    group_header = []
    # Provide a subset of the header that is just the ascii column headers.
    # The rest of the headers can be converted to energy group floats
    max_non_numeric_col = -1    
    for field in header:
        try:
	    num = float(field)
	# it doesn't convert, must be ascii
	except ValueError:
	    dose_meta_header.append(field)
	else:
	    group_header.append(num)
    if len(dose_meta_header) + len(group_header) != allrows.shape[1]:
        raise Exception('The header did not split nicely into non-numeric and numeric sections')
    return dose_meta_header, group_header, allrows
    
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
    Six record arrays, all with the same dtype:
    
    [('RESPONSE', 'a12'), ('PARTICLE', 'a12'),  
     ('GROUP', 'i4'), ('VALUES', 'f4', nrows)])

    Note that the values section of the dtype is an array (the column)

    These records contain all the response data:
    tr_dose      RESPONSE in 'dose', 'doseq' PARTICLE in ['all',particle_names]         GROUP = -1
    tr_let       RESPONSE = 'let'  	     PARTICLE in 'int', 'dif'                   GROUP in 1-700
    tr_flux      RESPONSE = 'flux','intflux' PARTICLE in particle_names                 GROUP in 1-100
    tr_neutflux  RESPONSE = 'neutflux'       PARTICLE in 'forward', 'backward', 'total' GROUP in 1-100

    Two more records are sums over flux groups
    tr_flux_part RESPONSE = 'flux_part'      PARTICLE in particle_names                 GROUP = -1
    tr_neutflux_part RESPONSE = 'neutflux_part'  PARTICLE in 'forward', 'backward', 'total' GROUP = -1

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
    tr_dose = np.zeros(ncols, dtype=response_type)
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

	part_flux_col        = np.sum(flux_data[:,part_j_start:part_j_start + num_flux_groups])
        tr_flux_part[j]  = ('flux_part', pname_j, -1, part_flux_col)

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

def get_particle_names(dose_header):
    # The dose header has two times the number of species, plus two total categories
    num_species = int((len(dose_header) - 2) / 2)
    # To fill in the PARTICLE field, extract name list depending
    # on the number of species, i.e. the radiation environment
    if 6 == num_species:
       particle_names = [x.split('_')[1] for x in dose_header[2:2+num_species]] 
    else:
       particle_names = ['particle_' + str(i) for i in range(num_species)]
    return particle_names

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

def get_meta_struct(data):
    """ Get a structured array consisting of the first num_meta columns
    of the passed in array.

    Parameters
    ----------
    data : numpy array with at least num_meta columns
        a two-dimensional array of floats where the 
	columns are over direction and the rows are the (u,v,w,depth) tuples
	for each direction

    Returns
    -------
    A structured array with the same number of tuples as the rows in the input
    data array

    """ 
    nrows = data.shape[0]
    tr_meta = np.zeros(nrows, dtype=meta_type) 
    for i in range(nrows):
        tr_meta[i] = tuple(data[i,:num_meta])
    return tr_meta

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

def set_response_type(nrows):
    global response_type
    response_type = [('RESPONSE', 'a12'), ('PARTICLE', 'a12'), \
	             ('GROUP', 'i4'), ('VALUES', 'f4', nrows)]
    return

def main():
    global response_type 
    args = parse_arguments()

    ########################
    # Get the header with the meta-data and dose,
    # the complete header, and the floating point row data
    dm_header, group_header, row_data = get_data(args.infile)
    nrows = row_data.shape[0]
    # Format for all response data
    set_response_type(nrows)
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
