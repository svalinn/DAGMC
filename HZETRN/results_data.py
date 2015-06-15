import argparse
import numpy as np

# Number of columns devoted to direction and depth
num_non_data = 4

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
    tr_rows = np.zeros(nrows, \
            dtype=[('id','i4'), ('u','f4'), ('v', 'f4'), \
	           ('w', 'f4'), ('depth', 'f4')])
    for i in range(nrows):
        tr_rows[i] = (i, data[i,0], data[i,1], data[i,2], data[i,3])
    return tr_rows

def get_dose_struct(header, data):
    """  Get a structured array of the named dose columns.
    
    Parameters
    ----------
    header : array of ascii words
        just the dose headers

    data : 2d array of floats
        the columns of the array should match the header

    To get the result column for particle "dose_1" you would
    type
    tp[tp['particle'] == 'dose_1']['result']
    or 
    result_col -f inputfile -p particle_name -o outfile
    To get the ids of the particle result
    tp[tp['particle'] == 'dose_1']['id']
    """
    # Note: the 'id' field is redundant and should be remvoed
    nrows = data.shape[0]
    tr_cols = np.zeros(len(header), \
              dtype=[('id', 'i4', nrows), ('particle','a12'), \
	             ('result', 'f4', nrows)])

    for i in range(len(header)):
        particle = header[i]
        tr_cols[i] = (range(nrows), particle, data[0:nrows,i])
    return tr_cols

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
    dif_let_start = 0
    let_size = 700
    tr_let_group = get_group_struct(dif_let_start, let_size, group_header)

    flux_start = dif_let_start + 700 + 700
    flux_size = 100
    tr_flux_group = get_group_struct(flux_start, flux_size, group_header)
    ##########################################################
    # There are always 
    """
    dif_let_end = dif_let_start + let_size + 1
    dif_let_data = row_data[:,dif_let_start:dif_let_end]

    tr_dif_let = np.zeros(let_size, \
              dtype=[('group_let_id', 'i4'), ('result', 'f4', nrows)])

    for i in range(let_size):
        # tr_let_group[tr_let_group['id'] == i]['energy'][0]
        tr_cols[i] = (i, data[:,i])

    int_let_start = dif_let_end
    int_let_end = int_let_start + let_size
    int_let_data = row_data[:,int_let_start:int_let_end]

    tr_int_let = np.zeros(let_size, \
                 dtype=[('group_let_id', 'i4'), ('result', 'f4', nrows)])
    """
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
