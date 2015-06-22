import argparse
import numpy as np

# Number of columns devoted to direction and depth
num_meta = 4

"""
    Read the lines of a  database file, stripping 
    off the header and statistics lines
"""
def get_data(infile):

    rows = np.array([])
    with open(infile) as fp:
        for index, line in enumerate(fp):
	    i = index + 1
	    if i == 1:
	        header = line.split()
	    else:
	        try: 
	            row = map(float, line.split())
	            if len(rows) == 0:
	                rows = row
	            else:
	                rows = np.vstack((rows, row))
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
	except ValueError:
	    dose_meta_header.append(field)
	else:
	    group_header.append(num)
    if len(dose_meta_header) + len(group_header) != rows.shape[1]:
        raise Exception('The header did not split nicely into non-numeric and numeric sections')
    return dose_meta_header, group_header, rows
    
def get_data_structs(dose_meta_header, group_header, rows):
    
    nrows = rows.shape[0]
    # remove meta_columns from data array
    data = rows[nrows,num_meta:]
    ncols = data.shape[1]
    # ToDo Check
    # Likewise, remove meta columns from the dose_meta_header
    dose_header = dose_meta_header[num_meta:]
    num_dose_response = len(dose_header) 
    # The dose header has two times the number of species, plus two total categories
    num_species = int((len(dose_header) - 2) / 2)

    # Construct the by-column record array
    dtype= [('RESPONSE', 'a12'), ('PARTICLE', 'a12'), \
	    ('GROUP', 'i4'), ('VALUES', 'f4', (nrows,1))])

    tr_dose = np.zeros(ncols,  dtype=dtype)
    tr_meta = get_meta_struct(data)
    
    for col in range(len(dose_header)):
        if col == 0:
            tr_all_cols[col] = ('Dose', 'All',  -1, data[:, col])
	elif col == 1:
            tr_all_cols[col] = ('Doseq', 'All', -1, data[:, col])
	else if col < num_species - 2:
            tr_all_cols[col] = ('Dose',  dose_header[col], -1, data[:, col])
        else:
            tr_all_cols[col] = ('Doseq', dose_header[col], -1, data[:, col])

    num_let_groups = 700
    tr_dif_let = np.zeros(num_let_groups, dtype=dtype)
    tr_int_let = np.zeros(num_let_groups, dtype=dtype)
    let_data = data[:,num_dose_response:]
    for i in range(num_let_groups):
       col = i
       tr_dif_let = np.zeros('let', 'All', i, let_data[:,col])  
    for i in range(num_let_groups):
       col = i + num_let_groups
       tr_int_let = np.zeros('let', 'All', i, let_data[:,col])  


    particle_names = dose_header[2:2+num_species]
    if 6 == num_species:
       particle_names = [x.split('_')[1] for x in particle_names] 
    else:
       particle_names = ['particle_' + str(i) for i in range(num_species)]

    num_flux_groups = 100
    num_flux_cols = num_flux_groups*num_species
    flux_data = let_data[:,num_let_groups*2:]
    intflux_data = flux_data[:,num_flux_cols:]
    # check data.shape[1] == 4 + len(dose_header) + 2*num_let_groups + 
    # flux_data.shape[1] + intflux_data.shape[1] + neutflux_data.shape[1]
    tr_flux = np.zeros(num_flux_cols, dtype)
    tr_int_flux = np.zeros(num_flux_cols, dtype)
    for ple_num in range(num_species):
        for group in range(num_flux_groups):
	    flux_col = group + num_flux_groups*ple_num
	    tr_flux[flux_col]     = ('flux', particle_names[ple_num], group, flux_data[:,flux_col])
	    tr_int_flux[flux_col] = ('int_flux', particle_names[ple_num], group, intflux_data[:,flux_col])
    
    tr_neut_flux = np.zeros(num_flux_cols, dtype)
    neutflux_data = intflux_data[:,num_flux_cols:]
    neutflux_names = ['forward', 'backward', 'total']
    for neutflux_type in range(3):
        for group in range(num_flux_groups):
	    flux_col = group + num_flux_groups*neutflux_type
            tr_neutflux[flux_col] = ('neut_flux', neutflux_names[neutflux_type], group, flux_data[:,flux_col])
	    
    return tr_all_data

"""
Argument parsing
returns : args: -d for the run directory
ref     : DAGMC/tools/parse_materials/dagmc_get_materials.py
"""
def parsing():
    
    min_col = num_non_data
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', action='store', dest='infile',
	help='The name of the .dat file with one line per direction.')

    parser.add_argument(
        '-c', action='store', dest='data_column', 
	help='The column, column range, or "a" for all columns.  Must be {} or greater'.format(num_non_data))

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
    """ Get a structured array consisting of the first 4 columns
    of the passed in array.

    Parameters
    ----------
    data : numpy array with at least 4 columns
        a two-dimensional array of floats where the 
	columns are over direction and the rows are the (u,v,w,depth) tuples
	for each direction

    Returns
    -------
    A structured array with the same number of tuples as the rows in the input
    data array

    """ 
    nrows = data.shape[0]
    tr_meta = np.zeros(nrows, \
            dtype=[('id','i4'), ('u','f4'), ('v', 'f4'), \
	           ('w', 'f4'), ('depth', 'f4')])
    for i in range(nrows):
        tr_meta[i] = (i, data[i,0], data[i,1], data[i,2], data[i,3])
    return tr_meta

# ToDo: Remove
def get_dose_struct(header, data):
    """  Get a structured array of the named dose columns.
    
    Parameters
    ----------
    header : array of ascii words
        just the dose headers are expected

    data : 2d array of floats
        the columns of the array should match the header

    Returns
    -------
    A structured array with a tuple for each element in the header
    

    Notes
    -----
    To get the result column for particle "dose_1" you would
    type
    tp[tp['particle'] == 'dose_1']['result']
    or 
    result_col -f inputfile -p particle_name -o outfile
    To get the ids of the particle result
    tp[tp['particle'] == 'dose_1']['id']

    """
    nrows = data.shape[0]
    tr_dose = np.zeros(len(header), \
              dtype=[('id', 'i4', nrows), ('particle','a12'), \
	             ('result', 'f4', nrows)])

    for i in range(len(header)):
        particle = header[i]
        tr_dose[i] = (range(nrows), particle,data[0:nrows,i])
    return tr_dose

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

def get_data_struct(header, data, let_size, flux_size):
    nrows = data.shape[0]
    tr_all_cols = np.zeros(data.shape[1], \
                  dtype = [('particle', 'a12'), ('group_id', 'i4'),
	                   ('response_col', 'f4', nrows)])
    tr_all_resp = np.zeros(data.shape[0]*data.shape[1], \
                  dtype = [('id', 'i4'),
                           ('particle', 'a12'), ('group_id', 'i4'),
	                   ('response', 'f4')])
    
    # The header contains Dose_all and Doseq_all
    ple_header = header[2:]
    # print 'ple_header', ple_header[:150]

    num_dose = len(header)
    for i in range(data.shape[1]):
	if i < num_dose:
            tr_all_cols = (header[i], -1, data[0:nrows,i])
	    """
	    for j in range(nrows):
                index = j+j*i
	        tr_all_resp[index] = (j, header[i], -1, data[j,i])
	    """
        elif i < num_dose + let_size*2:
	    group_index = (i - num_dose)%let_size
	    tr_all_cols = ('TOTAL', group_index, data[0:nrows,i])
	    """
	    for j in range(1,nrows+1):
                index = j+j*i
	        tr_all_resp[index] = (j, 'TOTAL', group_index, data[j,i])
            """
	else:
	    group_index = (i - num_dose*2)%flux_size
	    ple = group_index/flux_size
	    particle = ple_header[ple]
	    tr_all_cols = (particle, group_index, data[0:nrows,i])
	    """
	    for j in range(1,nrows+1):
                index = j+j*i
	        tr_all_resp[index] = (j, particle, group_index, data[j,i])
            """
    return tr_all_cols 

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

def get_dose_from_filedata(dm_header, row_data):
    """ Get a structured array from data returned from file

    Parameters
    ----------
    dm_header : array of ascii column descriptors
        col 0-3 'u','v','w','depth'
	col 4-N(radiation-environmnt-used) 
	       dose, dose-equivalent for each species
	       and dose, dose-equivalent for all
	Total number of columns is 4 + 2(N) + 2

    row_data : 2d numpy array of response data
        Columns 0 to (2(N) + 5) are the data for the named header.
	The rest of the columns are not used in this function

    Notes
    -----
        Each row in the row_data Numpy array is the data for  
	radiation in the direction indicated by the first three
	columns
    """
    dose_header = dm_header[4:]
    dose_data = row_data[:,4:len(dm_header)]
    return get_dose_struct(dose_header, dose_data)

def get_dose_from_file(infile):
    """ Get a structured array from data in a file.

    Parameters
    ----------
    filename : file with a header, 2 statistics rows, and otherwise
               row by column floats

    Reference
    ---------
    See get_dose_struct
    """
    dm_header, group_header, row_data = get_data(infile)
    dose_header = dm_header[4:]
    dose_data = row_data[:,4:len(dm_header)]
    return get_dose_struct(dose_header, dose_data)


def get_let_structs(group_header, group_data):
    """ Get the structures for 
    a) let energy boundaries
    b) dif let value sets (columns) for each energy group
    c) int let value sets for each energy group

    Parameters
    ----------
    group_header : 1-d float array
       As read from the data file, truncated to start with the first
       let energy group

    Returns
    -------
    Numpy record arrays for the let energy boundaries and for the dif_let and
    int_let results
    """
    nrows = group_data.shape[0]
    dif_let_start = 0
    let_size = 700
    tr_let_group = get_group_struct(dif_let_start, \
                                    let_size, group_header)

    dif_let_end = dif_let_start + let_size + 1
    tr_dif_let  = np.zeros(let_size, \
        dtype=[('group_let_id', 'i4'), ('value', 'f4', nrows)])
    for i in range(let_size):
	# First member of tuple could be i
        tr_dif_let[i] = (tr_let_group['id'][i], group_data[:,i])

    int_let_start = dif_let_end
    int_let_end   = int_let_start + let_size + 1
    tr_int_let    = np.zeros(let_size, \
        dtype=[('group_let_id', 'i4'), ('value', 'f4', nrows)])
    for i in range(let_size):
	# First member of tuple could be i
        tr_int_let[i] = (tr_let_group['id'][i], \
	                 group_data[:,int_let_start+i])

    return tr_let_group, tr_dif_let, tr_int_let


def get_let_from_file(infile):
    """ Get structure arrays from a file  

    Parameters
    ----------
    filename : file with a header, 2 statistics rows, and otherwise
               row by column floats
        File is the result of calling results_encoding on a HZETRN dataset

    Returns
    -------
    Numpy record arrays for the let energy boundaries and for the dif_let and
    int_let results
    """
    dm_header, group_header, row_data = get_data(infile)
    group_data = row_data[:,len(dm_header):]
    return get_let_structs(group_header, group_data)

def write_let_group(dm_header, group_header, row_data, outfile):
    tr_let_group, tr_dif_let, tr_int_let = \
        get_let_structs(group_header, row_data[:,len(dm_header):])
    with open(outfile, 'w') as f:
        np.savetxt(f, zip(tr_let_group['id'],tr_let_group['energy']))
    return


# Todo: remove
def get_groups_from_filedata(group_header):
    """ Get a structured array from the data returned by
    get_data.

    Parameters
    ----------
    """
    dif_let_start = 0
    let_size = 700
    tr_let_group = get_group_struct(dif_let_start, let_size, group_header)

    flux_start = dif_let_start + 700 + 700
    flux_size = 100
    tr_flux_group = get_group_struct(flux_start, flux_size, group_header)

    return tr_let_group, tr_flux_group

# Don't use
def get_groups_from_file(infile):
    """ Get a structured array from the data file.

    Parameters
    ----------
    filename : file with a header, 2 statistics rows, and otherwise
               row by column floats

    Reference
    ---------
    See get_group_struct

    """
    dm_header, group_header, row_data = get_data(infile)
    return get_groups_from_filedata(group_header)


def main():
    args = parsing()

    ########################
    # Get the header with the meta-data and dose,
    # the complete header, and the floating point row data
    dm_header, group_header, row_data = get_data(args.infile)
    nrows = row_data.shape[0]
    ##########################################################
    tr_meta = get_meta_struct(row_data[:,:4])
    ##########################################################
    dose_header = dm_header[4:]
    dose_data = row_data[:,4:len(dm_header)]
    tr_dose = get_dose_struct(dose_header, dose_data)
    ##########################################################
    # tr_let_group, tr_flux_group = get_groups_from_file_data(group_header)
    ##########################################################
    group_data = row_data[:,len(dm_header):]
    # tr_get_let_structs(group_header, group_data)

    ####################################################################3
    # Do all data cols at once
    tr_all_cols, tr_all_resp = get_data_struct(all_header[4:], 
                                               row_data[:,4:], 700, 100)
    with open('myfile.dat', 'w') as f:
	for row in tr_meta:
	    for el in row:
	        f.write(str(el) + ' ')
	    f.write('\n')

if __name__ == '__main__':
    main()
