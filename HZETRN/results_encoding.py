import subprocess
import argparse
import os
from itertools import islice
import glob
import numpy as np

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

class Results:

    directory=''

    def __init__(self, dir):
        self.directory = dir

    def index_from(self, name, index, numlines):
	vals = np.empty((1,numlines), float)
	lines = Results.last_lines(self, name, numlines)
	col = 0
	for line in lines:
	    words = line.split()
	    num = float(words[index])
	    vals[0,col] = num
	    col = col + 1
	return vals

    def last_lines(self, name, n):
	filepath = self.directory + name
        return check_last_n_lines(filepath, n)

class nx2(Results):

    def last_lines(self, name):
        result_line = Results.last_lines(self, name, 1)
	words = result_line[0].split()
	# We are skipping the first value, which is the depth
        vals = np.array(float(words[1]))
	return vals

class nx7(Results):
    def last_lines(self, name):
        "Read the last line containing 7 values."
        result_line = Results.last_lines(self, name, 1)
	words = result_line[0].split()
	# We are skipping the first value, which is the depth
	nums = np.array([float(x) for x in words[1:len(words)]])
	return nums

class n700x2(Results):

    def index_and_lines(self, name):

	lines = Results.last_lines(self, name, 700)

        rows = np.array([])
	for line in lines:
	    words = line.split()
	    # Turn the strings on the line into a float array
	    row = [float(x) for x in words]
	    if len(rows) == 0:
	        rows = row
	    else:
	        rows = np.vstack((rows, row))
        return rows

class n100xm(Results):
    "Read the last block of 100 lines of data, with m values per line."

    def index_and_vals(self, name):
	index = []
	vals = np.array([])
	lines = Results.last_lines(self, name, 100)
	for line in lines:
	    words = line.split()
	    # Grab the first column for the header
	    index.append(float(words[0]))
	    # We are skipping the first value, which is a flux bin
            # use list comprehension for generality and flexibility
	    nums = np.array([float(x) for x in words[1:len(words)]])
	    
	    if vals.size == 0:
	        vals = nums
	    else:
	        vals = np.vstack((vals, nums))

	all_vals_as_row = np.swapaxes(vals,0,1).ravel()
	full_index_hdr = index_hdr = np.array(index)
        # Concatenate the first col to be repeated for other raveled cols.
	for i in range(1,vals.shape[1]):
	   full_index_hdr = np.hstack((full_index_hdr, index_hdr))

	return (full_index_hdr, all_vals_as_row)
      
def reversed_lines(file):
    "Generate the lines of file in reverse order."
    part = ''
    for block in reversed_blocks(file):
	# Go from the end of the block to the beginning
	# c is a '\n' prefixed line
        for c in reversed(block):
            if c == '\n' and part:
	        # Since the block is reversed, the line is also reversed
                yield part[::-1]
                part = ''
            part += c
    # if part: yield part[::-1]
    if part: yield part


def reversed_blocks(file, blocksize=4096):
    "Generate blocks of file's contents in reverse order."
    file.seek(0, os.SEEK_END)
    here = file.tell()
    while 0 < here:
        delta = min(blocksize, here)
        here -= delta
        file.seek(here, os.SEEK_SET)
        yield file.read(delta)


def check_last_n_lines(file, n):
    "Return the last n lines of the file."
    lines=[]
    with open(file) as f:
        for line in islice(reversed_lines(f), n):
	    lines.append(line)
    # lines has oldest (last) line first; this returns
    # the last n lines in the same order as they are in the file
    return lines[::-1]

def get_names(file):
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
    dose_filename      = lines[0]
    doseq_filename     = lines[1]
    dose_ple_filename  = lines[2]
    doseq_ple_filename = lines[3]
    diflet_filename    = lines[4]
    intlet_filename    = lines[5]
    flux_filename      = lines[6]
    intflux_filename   = lines[7]
    neutron_flux_filename = lines[8]

    return 

def depth_header_ary():
    "Return a header for the depth data as a 1 x 17 array"
    dh = []
    dh = [ ['{0: >13}'.format('X'), 
              '{0: >12}'.format('Y'),
              '{0: >12}'.format('Z'),
              '{0: >12}'.format('Dose_All'),
              '{0: >12}'.format('Doseq_All'),
              '{0: >12}'.format('dose-neutron'),
              '{0: >12}'.format('dose-proton'),
              '{0: >12}'.format('dose-deut.'),
              '{0: >12}'.format('dose-trit.'),
              '{0: >12}'.format('dose-He3'),
              '{0: >12}'.format('dose-He4'),
              '{0: >12}'.format('doseq-neut'),
              '{0: >12}'.format('doseq-proton'),
              '{0: >12}'.format('doseq-deut'),
              '{0: >12}'.format('doseq-trit'),
              '{0: >12}'.format('doseq-He3'),
              '{0: >12}'.format('doseq-He4')]]
     
    # return dh
    return np.array(dh)

def get_filepaths(r):
    "Convenience function to construct filepaths and check for existence"
    global dose_filename, doseq_filename, dose_ple_filename, doseq_ple_filename
    global diflet_filename, intlet_filename
    global flux_filename, intflux_filename
    global neutron_flux_filename 
    
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
    
    return ret

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
    
    return args

def main():

    # Setup: parse the the command line parameters
    args = parsing()

    # Set up the names of the data files without them being in this file
    get_names("names.txt")
    
    # Original Working Directory - will be used for the results file
    owd = os.path.dirname(os.path.abspath(__file__)) + '/'
    outfile = owd + args.data_file
    print 'Writing to', outfile

    # Directory containing the ray subdirectories
    run_path = args.run_dir + '/data/'

    os.chdir(run_path)
    depth_data_reader = nx2('')
    depth_particle_reader = nx7('')
    dif_int_block_reader  = n700x2('')
    flux_block_reader     = n100xm('')

    # Numpy Arrays
    all_values  = np.array([])
    full_header = np.array([])
    
    # Make a header array of words for the depth particles
    hdr_depth = depth_header_ary()
    
    ray_subdirs = glob.glob('*')
    for r in ray_subdirs:

	# String list: this is what we send to the data line
        ray = r.split('_') 
	# If any files are missing the dictionary will be empty
        datapaths = get_filepaths(r)

        if len(datapaths) != 0:
	    ray_vec = np.array([float(x) for x in ray])
	    print ray_vec

	    #######################################################
	    # Depth data
            # Add dose and doseq at depth for all or 6 particles
	    dose      = depth_data_reader.last_lines(datapaths['dose'])
	    doseq     = depth_data_reader.last_lines(datapaths['doseq'])
	    dose_ple  = depth_particle_reader.last_lines(datapaths['dose_part'])
	    doseq_ple = depth_particle_reader.last_lines(datapaths['doseq_part'])
	   
	    all_values_by_ray = np.hstack((dose, doseq, dose_ple, doseq_ple))
	    #######################################################
	    # LET - 700 values per ray for dif and int each
	    # col 0 is header, col 1 is data
	    dif_let = dif_int_block_reader.index_and_lines(datapaths['dif_let'])
	    int_let = dif_int_block_reader.index_and_lines(datapaths['int_let'])
	    
	    all_values_by_ray = np.hstack((all_values_by_ray, dif_let[...,1], int_let[...,1]))

	    # 100 rows by 7 or 4 columns
	    flux_hdr, flux_nums         = flux_block_reader.index_and_vals(datapaths['flux'])
	    intflux_hdr, intflux_nums   = flux_block_reader.index_and_vals(datapaths['int_flux'])
	    neutflux_hdr, neutflux_nums = flux_block_reader.index_and_vals(datapaths['neutron'])

	    all_values_by_ray = np.hstack((all_values_by_ray, flux_nums, intflux_nums, neutflux_nums))

	    # The first time through write the header and start the overall data stack
            if all_values.size == 0:
		# 700 numbers, from 1st column, append to header 2x, for dif and int
		let_hdr_nums = np.hstack((dif_let[...,0],int_let[...,0]))
		hdr_nums = np.hstack((let_hdr_nums, flux_hdr, intflux_hdr, neutflux_hdr))

		with open(outfile, 'w') as f:
		    np.savetxt(f, hdr_depth, fmt='%s', newline='')
		    print 'hdr_depth shape', hdr_depth.shape
		    
		with open(outfile, 'a') as f:
		    np.savetxt(f, hdr_nums, fmt='%13.6E', newline='')
	            f.write('\n')
		    print 'hdr_nums shape', hdr_nums.shape

	 	all_values = all_values_by_ray

	    else:
		# The pure-numeric array for statistics
	        all_values = np.vstack((all_values, all_values_by_ray))

            with open(outfile, 'a') as f:
	       np.savetxt(f, ray_vec, fmt='%13.6f', newline='')
	       np.savetxt(f, all_values_by_ray, fmt='%13.6E', newline='')
	       f.write('\n')

    # NOTES
    # All data files must be present
    # print 'Saving ray_vec', ray_vec, 'and all_value array, shaped', all_values.shape
    ###############################################################

    # Averages
    ave_all_values = np.average(all_values, axis=0)
    with open(outfile, 'a') as f:
        np.savetxt(f, ['Average', 'over', 'direction'], fmt='%13s', newline='')
        np.savetxt(f, ave_all_values, fmt='%13.6E', newline='')
	f.write('\n')

    # Std. Dev.
    std_all_values = np.std(all_values, axis=0)
    with open(outfile, 'a') as f:
        np.savetxt(f, ['Std.Dev.', 'over', 'direction'], fmt='%13s', newline='')
        np.savetxt(f, std_all_values, fmt='%13.6E', newline='')
        f.write('\n')
        
    ###############################################################

if __name__ == '__main__':
    main()
