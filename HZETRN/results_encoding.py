import subprocess
import argparse
import os
import sys
import glob
import numpy as np
import read_utils as ru
import logging

# Global names to be read from a separate text file
dose_filename  = ''
doseq_filename = ''
dose_ple_filename  = ''
doseq_ple_filename = ''
diflet_filename = '' 
intlet_filename = ''
flux_filename = ''
intflux_filename = ''
neutron_flux_filename = ''
slabs_to_ref_filename = 'slabs_to_ref.txt'

###############################
# globals
config_log  = 'config.log'

class Results:
    """
    Base class for reader classes
    3 + depth_header(env) + let700_headers + flux100_headers(env) = {2917, 13623}
    """

    directory=''

    def __init__(self, dir):
        self.directory = dir

    def words_last_line(self, name, n):
        result_line = Results.last_lines(self, name, 1)
	return result_line[0].split()

    def last_lines(self, name, n):
	filepath = self.directory + name
        return ru.check_last_n_lines(filepath, n)


class depth_by_m(Results):
    """
    Reader for files with lines at depth; 
    Old version: last depth is of interest.
    Backscatter: 
        line of interest = HEADER_LINES + slabs_to_ref + 1
	HEADER_LINES = 1
    Four files are read like this: dose and dose-equivalant, overall
        and per particle
    'm' varies based on per particle or not; if per particle file,
        it is environment-dependent
    
    The header column is repeated once for each additional column.
    
    Contribution to depth header = dose + doseq + dose(env) + doseq(env):
    1 + 1 + 6 + 6    = 14
    - or - 
    1 + 1 + 59 + 59  = 120
    """
    HEADER_LINES = 1
    def line_at_n(self, name, n):
        "Read the nth line as floating point values." 

	# Because of the enumerate, deth_line is 0-based, so don't add 1
	ref_line = self.HEADER_LINES + n
	words = []
        fp = open(name)
	for i, line in enumerate(fp):
	    if i == ref_line:
	       words = line.split()
	       break 
	fp.close()     
	row = np.array([float(x) for x in words])

        # return the first value separately
	return round(row[0],6), np.array([row[1:len(row)]])

    def last_lines(self, name):
        "Read all but the first value in the last line."
	words = Results.words_last_line(self, name, 1)
	# We are skipping the first value, which is the depth
	# return  np.array([float(x) for x in words[1:len(words)]])
	row = np.array([float(x) for x in words[1:len(words)]])
	return np.array([row])

class let700x2(Results):
    """
    let700x2
    Reader for files with blocks of 700 energy boundaries at depth, 
    first column is header, second column is data.  
    Contribution to header is NOT env-dependent.  Two files have 
    this characteristic, diflet and intlet:

    2 x first_col_of_last_block = 1400

    Old version: read last block always
    Backscatter: 
        read nth block starting at ref_line, where n = slabs_to_ref:
	HEADER_LINES + slabs_to_ref*(700+SEPARATOR_LINES) + 1

	HEADER_LINES = 2
	SEPARATOR_LINES = 1

    """
    HEADER_LINES    = 2
    BLOCK_SIZE      = 700
    SEPARATOR_LINES = 1

    def block_at_n(self, name, n):
        "Read a block of 700x2 values."

	ref_line = self.HEADER_LINES + n*(self.BLOCK_SIZE + self.SEPARATOR_LINES)  
	words = []
	rows = np.array([])
	keyword = 'Depth'
        fp = open(name)
	for index, line in enumerate(fp):
	    i = index + 1
	    if i == ref_line:
	        words = line.split()
		if words[0] == keyword:
		    depth = float(words[len(words) - 1])
	        else:
		    msg = "First word of {} at line {} of {} is not {}".format(line, i,name, keyword)
		    print msg
		    raise ValueError(msg)
	    if i > ref_line and i <= ref_line + self.BLOCK_SIZE:
	        words = line.split()
		# print i, name, words
	        # Convert the strings on the line into a float array
	        row = [float(x) for x in words]
	        if len(rows) == 0:
	           rows = row
	        else:
	           rows = np.vstack((rows, row))
	fp.close()     
        return round(depth,6), np.array([rows])
        
	
    def last_lines(self, name):

	lines = Results.last_lines(self, name, 700)

        rows = np.array([])
	for line in lines:
	    words = line.split()
	    # Convert the strings on the line into a float array
	    row = [float(x) for x in words]
	    if len(rows) == 0:
	        rows = row
	    else:
	        rows = np.vstack((rows, row))

        # index and data columns are extracted by caller
        return np.array([rows])


class flux100xm(Results):
    """
    flux100xm
    Reader for files with blocks of 100 flux values at depth, 
    first column is header, remaining columns are data.  
    Contribution to header IS env-dependent for two files and 
    NOT for a third.  The one that is not env-dependent has
    3 columns besides the header column.  
    
    The header column is repeated once for each data column.

    2 x first_col_of_last_block X ple(env) + 100 * 3 = 
    2 * 100 * 6 + 300 =  1500
    -or-
    2 * 100 * 59 + 300 = 12100

    Backtrace Modifications
    _______________________
    Instead of the last block we must read the nth block, where
    n = 1 + slabs-to-ref

    and

    slabs-to-ref is read out of the file slabs_to_ref.txt in each
    directional subdirectory.

        read block starting at 
	HEADER_LINES + slabs_to_ref*(BLOCK_LINES + SEPARATOR_LINES) + 1

	BLOCK_SIZE = 100
	HEADER_LINES = 2
	SEPARATOR_LINES = 1

    """
    keyword         = 'Depth'
    HEADER_LINES    = 2
    BLOCK_SIZE      = 100
    SEPARATOR_LINES = 1


    """
    def block_at_n(self, name, n, key):
        self.keyword = key
	return block_at_n(self, name, n)
    """
    def block_at_n(self, name, n, key):
        "Read a block of 100xm values."

	self.keyword = key
	ref_line = self.HEADER_LINES + n*(self.BLOCK_SIZE + self.SEPARATOR_LINES)  
	words = []
	rows = np.array([])
        fp = open(name)
	for index, line in enumerate(fp):
	    i = index + 1
	    if i == ref_line:
	        words = line.split()
		if words[0] == self.keyword:
		    depth = float(words[len(words) - 1])
	        else:
		    msg = "First word in {} at line {} is not {}".format(name, i,keyword)
		    raise ValueError(msg)
	    if i > ref_line and i <= ref_line + self.BLOCK_SIZE:
	        words = line.split()
	        # Convert the strings on the line into a float array
	        row = [float(x) for x in words]
	        if len(rows) == 0:
	           rows = row
	        else:
	           rows = np.vstack((rows, row))
	fp.close()     

        index = rows[:,0]
	num_reps = rows.shape[1] - 1
	# Concatenate the indices once for each remaining column
	full_index_header = np.array([list(index) * num_reps])

	all_nums = np.array([np.swapaxes(rows,0,1).ravel()])
	return (round(depth,6), full_index_header, all_nums)
	
        
    # Set once
    num_species=0

    def index_and_vals(self, name):
        "Read the last block of 100 lines of data, with m values per line."
	index = []
        rows = np.array([])
	lines = Results.last_lines(self, name, 100)

	for line in lines:
	    words = line.split()
	    # Grab the first column for the header
	    index.append(float(words[0]))
	    # We are skipping the first value, which is a flux bin
            # use list comprehension for generality and flexibility
	    row = np.array([float(x) for x in words[1:len(words)]])
	    if len(rows) == 0:
	        rows = row
	    else:
	        rows = np.vstack((rows, row))
        # Here the array of data taken from the last 100 lines of the file, 
	#     less the leftmost column, is transposed so that each 100 values 
	#     that was in a column is in a row.  Then it is unraveled,
	#     so that all the rows are end-to-end.
	all_particle_vals_as_row = np.swapaxes(rows,0,1).ravel()
	
	full_index_hdr = index_hdr = np.array(index)
        # Concatenate the first col to be repeated for other raveled cols.
	# rows.shape[1] = dimension 1 = num_cols
	for i in range(1,rows.shape[1]):
	   full_index_hdr = np.hstack((full_index_hdr, index_hdr))

        # This is an aside for later header calculation
	if self.num_species == 0:
            self.num_species = rows.shape[1]

	# Change the shape to [1,something] rather than [something,]
	all_nums = np.array([all_particle_vals_as_row])
	return (np.array([full_index_hdr]), all_nums)
      
def get_global_names(file):
    "Fill in the global filenames from a text file."
    global dose_filename, doseq_filename, dose_ple_filename, doseq_ple_filename
    global diflet_filename, intlet_filename
    global flux_filename, intflux_filename
    global neutron_flux_filename 

    lines = []
    with open(file,'r') as f:
	# Get read of the newline
        for line in f:
            lines.append(line.split('\n')[0])
    dose_filename         = lines[0]
    doseq_filename        = lines[1]
    dose_ple_filename     = lines[2]
    doseq_ple_filename    = lines[3]
    diflet_filename       = lines[4]
    intlet_filename       = lines[5]
    flux_filename         = lines[6]
    intflux_filename      = lines[7]
    neutron_flux_filename = lines[8]

    return 

def depth_header_ary(species):
    "Return a header for the depth data."
    header_list = ['{0: >12}'.format('X'), 
              '{0: >12}'.format('Y'),
              '{0: >12}'.format('Z'),
	      #'{0: >12}'.format('Depth'),
              '{0: >12}'.format('Dose_All'),
              '{0: >12}'.format('Doseq_All') ]
    if species == 6:
        header_list +=  ['{0: >12}'.format('dose-neutron'),
              '{0: >12}'.format('dose_proton'),
              '{0: >12}'.format('dose_deut.'),
              '{0: >12}'.format('dose_trit.'),
              '{0: >12}'.format('dose_He3'),
              '{0: >12}'.format('dose_He4'),
              '{0: >12}'.format('doseq_neut'),
              '{0: >12}'.format('doseq_proton'),
              '{0: >12}'.format('doseq_deut'),
              '{0: >12}'.format('doseq_trit'),
              '{0: >12}'.format('doseq_He3'),
              '{0: >12}'.format('doseq_He4')]
    elif species == 59:
        for i in range(1,60):
	    id = "{0: >12}".format('dose_%s' % i)
	    header_list.append(id)
        for i in range(1,60):
	    iq = "{0: >12}".format('doseq_%s' % i)
	    header_list.append(iq)
    header = np.array([header_list])
    return header

def get_filepaths(r):
    "Convenience function to construct filepaths and check for existence"
    global dose_filename, doseq_filename, dose_ple_filename, doseq_ple_filename
    global diflet_filename, intlet_filename
    global flux_filename, intflux_filename
    global neutron_flux_filename 
    global slabs_to_ref_filename
    
    ret = {}
    ret['dose'] = r + '/' + dose_filename 
    if not os.path.exists(ret['dose']):
        return {}
    ret['doseq'] = r + '/' + doseq_filename
    if not os.path.exists(ret['doseq']):
        return {}
    ret['dose_part'] = r + '/' + dose_ple_filename 
    if not os.path.exists(ret['dose_part']):
        return {}
    ret['doseq_part'] = r + '/' + doseq_ple_filename
    if not os.path.exists(ret['doseq_part']):
        return {}

    ret['dif_let'] = r + '/' + diflet_filename
    if not os.path.exists(ret['dif_let']):
        return {}
    ret['int_let'] = r + '/' + intlet_filename
    if not os.path.exists(ret['int_let']):
        return {}

    ret['flux'] = r + '/' + flux_filename
    if not os.path.exists(ret['flux']):
        return {}
    ret['int_flux'] = r + '/' + intflux_filename
    if not os.path.exists(ret['int_flux']):
        return {}
    ret['neutron'] = r + '/' + neutron_flux_filename
    if not os.path.exists(ret['neutron']):
        return {}

    ret['slabs_to_ref'] = r + '/' + slabs_to_ref_filename
    if not os.path.exists(ret['slabs_to_ref']):
        return {}
    
    return ret

"""
    slabs_to_ref = get_slabs_to_ref(datapaths['slabs_to_ref'])
"""
def get_slabs_to_ref(filename):
    with open(filename) as f:
        return int(f.readline())


"""
    Check that all elements are the same
    Ref: 
    http://stackoverflow.com/questions/3844801/check-if-all-elements-in-a-list-are-identical
"""
def checkEqual(lst):
    return lst[1:] == lst[:-1]
    
"""
Argument parsing
returns : args: -d for the run directory
ref     : DAGMC/tools/parse_materials/dagmc_get_materials.py
"""
def parsing():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-d', action='store', dest='run_dir',
	help='The name of the holding directory for all hzetrn runs')

    parser.add_argument(
        '-o', action='store', dest='data_file', 
	help='The relative path to the output file')
    parser.set_defaults(data_file='tsta.dat')

    args = parser.parse_args()
    if not args.run_dir:
        raise Exception('Run directory not specified. [-d] not set')
    
    if not args.data_file:
        raise Exception('Data file (output) not specified. [-o] not set')
    return args

"""
    results_encoding.py
    - Log this in the run_directory config.log
"""
def main():
    # Encode the transport results into a single file by direction.

    global config_log

    args = parsing()

    # Set up the names of the data files 
    get_global_names("names.txt")
    
    # Original Working Directory - will be used for the results file
    owd = os.path.dirname(os.path.abspath(__file__)) + '/'
    dat_outfile = owd + args.data_file
    print 'Writing to', dat_outfile

    ################################
    # Logging
    cur_path = os.path.dirname(os.path.abspath(__file__)) + '/'
    run_path = cur_path + args.run_dir + '/'
    config_path = run_path + config_log
    logging.basicConfig(filename=config_path, level=logging.INFO, format='%(asctime)s %(message)s')
    message = 'python results_encoding.py -d ' + args.run_dir + \
              ' -o ' + args.data_file 
    logging.info(message)
    ################################

    # ToDo: it may or may not be ok to hardcode this
    # Directory containing the ray subdirectories
    run_path = args.run_dir + '/data/'

    os.chdir(run_path)
    depth_particle_reader = depth_by_m('')
    dif_int_block_reader  = let700x2('')
    flux_block_reader     = flux100xm('')

    # Numpy Arrays
    all_values  = np.array([])
    full_header = np.array([])
    empty_ray_vecs = np.array([])
    
    rad_env = ''
    rad_env_file = owd + run_path + 'rad_env.txt'
    with open(rad_env_file, 'r') as f:
        rad_env = f.readline()

    print rad_env
    species = 0
    zero_fields = np.array([])
    if rad_env == 'spe':
       print 'species is 6'
       species = 6
       zero_fields  = np.zeros(2914)
    else:
       print 'species is 59'
       species = 59
       zero_fields = np.zeros(13620)

    # This picks up files that have at least one digit
    ray_subdirs = glob.glob('*[0-9]*')
    i = 0

    for r in ray_subdirs:
	# String list: this is what we send to the data line
        ray = r.split('_') 
	ray_vec = np.array([float(x) for x in ray])

	# If any files are missing the dictionary will be empty
        datapaths = get_filepaths(r)
	print datapaths

        # Save 0-data if there isn't any
        if len(datapaths) == 0:
	    print i, ray_vec, "no data"
	    # Data, including direction vector
	    # Header has already been written
	    if all_values.size != 0:
                with open(dat_outfile, 'a') as f:
	           np.savetxt(f, ray_vec, fmt='%12.6f', newline=' ')
	           np.savetxt(f, zero_fields, fmt='%12.6E', newline=' ')
	           f.write('\n')
	    # Hold off til the header is written
            else:
		if empty_ray_vecs.size == 0:
		    empty_ray_vecs = ray_vec
		    print 'empty_ray_vecs size', empty_ray_vecs.size
		else:
	            empty_ray_vecs = np.vstack((empty_ray_vecs, ray_vec))
		    print 'empty_ray_vecs size', empty_ray_vecs.size
	# this ray direction has data  
        else:
	    print i, ray_vec
	    slabs_to_ref = get_slabs_to_ref(datapaths['slabs_to_ref'])
	    print 'slabs_to_ref for ', r, slabs_to_ref

	    #######################################################
	    # Depth data
            # Add dose and doseq at depth for environment-determined particle set
	    """
	    dose      = depth_particle_reader.last_lines(datapaths['dose'])
	    doseq     = depth_particle_reader.last_lines(datapaths['doseq'])
	    dose_ple  = depth_particle_reader.last_lines(datapaths['dose_part'])
	    doseq_ple = depth_particle_reader.last_lines(datapaths['doseq_part'])
	    """

	    dose_depth, dose = depth_particle_reader.line_at_n(datapaths['dose'], slabs_to_ref)
	    doseq_depth, doseq = depth_particle_reader.line_at_n(datapaths['doseq'], slabs_to_ref)
	    dose_ple_depth, dose_ple = depth_particle_reader.line_at_n(datapaths['dose_part'], slabs_to_ref)
	    doseq_ple_depth, doseq_ple = depth_particle_reader.line_at_n(datapaths['doseq_part'], slabs_to_ref)

	    print 'Dose depth', dose_depth
	    if not checkEqual([dose_depth, doseq_depth, dose_ple_depth, doseq_ple_depth]):
		sys.exit("Problem with dose files for direction {}".format(r))
	        
	    # ray_vec.append(dose_depth)
	    dose_values = np.hstack((dose, doseq, dose_ple, doseq_ple))
	    #######################################################
	    # LET - 700 values per ray for dif and int each
	    # col 0 is header, col 1 is data
	    # dif_let = dif_int_block_reader.last_lines(datapaths['dif_let'])
	    # int_let = dif_int_block_reader.last_lines(datapaths['int_let'])

	    depth_let, dif_let = dif_int_block_reader.block_at_n(datapaths['dif_let'], slabs_to_ref)
	    depth_int, int_let = dif_int_block_reader.block_at_n(datapaths['int_let'], slabs_to_ref)

	    print 'LET depth', depth_let
	    if not checkEqual([dose_depth, depth_let, depth_int]):
	        sys.exit("Problem with let files for direction {}".format(r))
	    
	    # Pick of col 
	    let_values = np.hstack((dif_let[...,1], int_let[...,1]))
	    #######################################################
	    # FLUX - 100 rows by 60(gcr)/7(spe) or 4 columns
	    # col 0, repeated, is header for other columns
	    # flux_hdr, flux_nums         = flux_block_reader.index_and_vals(datapaths['flux'])
	    # intflux_hdr, intflux_nums   = flux_block_reader.index_and_vals(datapaths['int_flux'])
	    # neutflux_hdr, neutflux_nums = flux_block_reader.index_and_vals(datapaths['neutron'])

	    depth_flux1, flux_hdr, flux_nums         = flux_block_reader.block_at_n(datapaths['flux'], slabs_to_ref, 'Slab')
	    depth_flux2, intflux_hdr, intflux_nums   = flux_block_reader.block_at_n(datapaths['int_flux'], slabs_to_ref, 'Depth')
	    depth_flux3, neutflux_hdr, neutflux_nums = flux_block_reader.block_at_n(datapaths['neutron'], slabs_to_ref, 'Slab')
	    if not checkEqual([dose_depth, depth_flux1, depth_flux2, depth_flux3]):
	        sys.exit("Problem with flux files for direction {}".format(r))

	    flux_values = np.hstack((flux_nums, intflux_nums, neutflux_nums))

	    # Put the dose, let and flux values together into a row
 	    all_values_by_ray = np.hstack((dose_values, let_values, flux_values))

	    # For the the first valid direction write the header and start the overall data stack
            if all_values.size == 0:
		
                # Make a header array of words for the depth particles
		# Must be done after number of species is known
                hdr_depth = depth_header_ary(species)

		# Extract the let portion of the header
		let_hdr_nums = np.hstack((dif_let[...,0],int_let[...,0]))

		print let_hdr_nums.shape, flux_hdr.shape, intflux_hdr.shape, neutflux_hdr.shape
		hdr_nums = np.hstack((let_hdr_nums, flux_hdr, intflux_hdr, neutflux_hdr))

		# Header words, previously formatted
		with open(dat_outfile, 'w') as f:
		    np.savetxt(f, hdr_depth, fmt='%s', newline=' ')
		    
		# Header Numbers, format in place
		with open(dat_outfile, 'a') as f:
		    np.savetxt(f, hdr_nums, fmt='%12.6E', newline='')
	            f.write('\n')
		    # Write out any of the ray vectors with no data, 
		    #    now that the header has been written
	            if empty_ray_vecs.size != 0:
		        for ray_vec in empty_ray_vecs:
	                    np.savetxt(f, ray_vec,     fmt='%12.6f', newline=' ')
	                    np.savetxt(f, zero_fields, fmt='%12.6E', newline=' ')
	                    f.write('\n')

	 	all_values = all_values_by_ray
            # end if all_values.size == 0:

	    else:
		# The pure-numeric array for statistics
	        all_values = np.vstack((all_values, all_values_by_ray))

	    # Data, including direction vector
            with open(dat_outfile, 'a') as f:
	       np.savetxt(f, ray_vec, fmt='%12.6f', newline=' ')
	       np.savetxt(f, all_values_by_ray, fmt='%12.6E', newline='')
	       f.write('\n')

        i = i + 1

    # NOTES
    # .....
    ###############################################################

    # Averages, Std. Dev.
    ave_all_values = np.average(all_values, axis=0)
    std_all_values = np.std(all_values, axis=0)

    with open(dat_outfile, 'a') as f:
        np.savetxt(f, [['Average', 'over', 'direction']], fmt=['%12s', '%12s', '%12s'], newline=' ')
        np.savetxt(f, [ave_all_values], fmt='%12.6E', newline='\n')
        np.savetxt(f, [['Std.Dev.', 'over', 'direction']], fmt=['%12s', '%12s', '%12s'], newline=' ')
        np.savetxt(f, [std_all_values], fmt='%12.6E', newline='\n')

    ###############################################################

if __name__ == '__main__':
    main()
