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
        result_line = Results.last_lines(self, name, 1)
	words = result_line[0].split()
	# We are skipping the first value, which is the depth
	nums = np.array([float(x) for x in words[1:len(words)]])
	return nums

class n700x2(Results):

    def index_from(self, name, index):
        return Results.index_from(self, name, index, 700)

    def last_lines(self, name):
        vals = self.index_from(name, 1)
	return vals


class n100xm(Results):
    "Read the last block of 100 lines of data, with m values per line"

    def index_from(self, name, index):
	return Results.index_from(self, name, index,  100)

    def last_lines(self, name):
	vals = np.array([])
	nums = np.array([])
	lines = Results.last_lines(self, name, 100)
	for line in lines:
	    words = line.split()
	    # We are skipping the first value, which is a flux bin
            # use list comprehension for generality and flexibility
	    nums = np.array([float(x) for x in words[1:len(words)]])
	    if vals.size == 0:
	        vals = nums
	    else:
	        vals = np.vstack((vals, nums))
	return vals
      
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
    dh = [ ['{0: >12}'.format('X'), 
              '{0: >15}'.format('Y'),
              '{0: >15}'.format('Z'),
              '{0: >15}'.format('Dose_All'),
              '{0: >15}'.format('Doseq_All'),
              '{0: >15}'.format('dose-neutron'),
              '{0: >15}'.format('dose-proton'),
              '{0: >15}'.format('dose-deut.'),
              '{0: >15}'.format('dose-trit.'),
              '{0: >15}'.format('dose-He3'),
              '{0: >15}'.format('dose-He4'),
              '{0: >15}'.format('doseq-neutron'),
              '{0: >15}'.format('doseq-proton'),
              '{0: >15}'.format('doseq-deut.'),
              '{0: >15}'.format('doseq-trit.'),
              '{0: >15}'.format('doseq-He3'),
              '{0: >15}'.format('doseq-He4')]]
     
    return np.array(dh)


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
    
    print "Saving processed data to file", args.data_file
    return args

def main():
    global dose_filename, doseq_filename, dose_ple_filename, doseq_ple_filename
    global diflet_filename, intlet_filename
    global flux_filename, intflux_filename
    global neutron_flux_filename 

    # Setup: parse the the command line parameters
    args = parsing()

    # Set up the names of the data files without them being in this file
    get_names("names.txt")
    
    # Original Working Directory - will be used for the results file
    owd = os.path.dirname(os.path.abspath(__file__)) + '/'
    outfile = owd + args.data_file
    print 'writing to', outfile

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
    # flux_nums_for_all = np.array([])
    
    # Make a header array of words for the depth particles
    dhl = depth_header_ary()
    # dh = depth_header()
    
    ray_subdirs = glob.glob('*')
    for r in ray_subdirs:
	# Depth files
	dose_filepath       = r + '/' + dose_filename 
	doseq_filepath      = r + '/' + doseq_filename
	dose_part_filepath  = r + '/' + dose_ple_filename 
	doseq_part_filepath = r + '/' + doseq_ple_filename
	# LET files
	dif_filepath        = r + '/' + diflet_filename
	int_filepath        = r + '/' + intlet_filename
	# Flux files
	flux_filepath       = r + '/' + flux_filename
	intflux_filepath    = r + '/' + intflux_filename
        neutron_flux_filepath = r + '/' + neutron_flux_filename

	# String list: this is what we send to the data line
        ray = r.split('_') 
	print ray

	# ray_vec = np.array([float(x) for x in ray])

	# ToDo: perform this check earlier
	if os.path.exists(dose_filepath) and os.path.exists(doseq_filepath) and \
	   os.path.exists(dose_part_filepath) and os.path.exists(doseq_part_filepath) and \
	   os.path.exists(dif_filepath) and os.path.exists(int_filepath) and \
	   os.path.exists(flux_filepath) and os.path.exists(intflux_filepath) \
	                                 and os.path.exists(neutron_flux_filepath):
	   
            # Add dose and doseq at depth for all or 6 particles
	    h_dose      = depth_data_reader.last_lines(dose_filepath)
	    i_doseq     = depth_data_reader.last_lines(doseq_filepath)
	    j_dose_ple  = depth_particle_reader.last_lines(dose_part_filepath)
	    k_doseq_ple = depth_particle_reader.last_lines(doseq_part_filepath)
	   
	    #######################################################
	    # LET - 700 values per ray for dif and int each
	    l_dif_let = dif_int_block_reader.last_lines(dif_filepath)
	    m_int_let = dif_int_block_reader.last_lines(int_filepath)

	    # 100 rows by 6 columns
            nnums = flux_block_reader.last_lines(flux_filepath)
	    onums = flux_block_reader.last_lines(intflux_filepath)
            pnums = flux_block_reader.last_lines(neutron_flux_filepath)

	    # Flip the rows and columns: shape (100, n) -> (n, 100)
	    rnnums = np.swapaxes(nnums,0,1)
	    ronums = np.swapaxes(onums,0,1)
	    rpnums = np.swapaxes(pnums,0,1)
	    
	    all_values_by_ray = np.hstack((ray, h_dose, i_doseq, j_dose_ple, k_doseq_ple))
	    all_values_by_ray = np.hstack((all_values_by_ray,l_dif_let.ravel(), m_int_let.ravel()))
	    all_values_by_ray = np.hstack((all_values_by_ray, rnnums.ravel(), ronums.ravel(), rpnums.ravel()))

	    # The first time through write the header and start the overall data stack
            if all_values.size == 0:
		# 700 numbers, from 1st column, append to header 2x, for dif and int
	        let_hdr_nums = dif_int_block_reader.index_from(dif_filepath, 0)
	        flux_depth_nums = flux_block_reader.index_from(flux_filepath, 0)
		# 700 + 700 + 100
		flux_nums_for_all = flux_depth_nums
		for i in range(0,14):
		    flux_nums_for_all = np.hstack((flux_nums_for_all, flux_depth_nums))
		# Create one very long line: first write, then append, no newlines
		hdr_nums = np.hstack((let_hdr_nums, let_hdr_nums, flux_nums_for_all))
		print 'hdr_nums shape', hdr_nums.shape
		# print 'hdr_nums', np.swapaxes(hdr_nums,0,1)

		with open(outfile, 'w') as f:
		    np.savetxt(f, dhl, fmt='%s')
		    
		with open(outfile, 'a') as f:
		    np.savetxt(f, hdr_nums, fmt='%13.6E')

		all_values = all_values_by_ray

	    else:
	        all_values = np.vstack((all_values, all_values_by_ray))

	    # Write one ray at a time; could also store in memory 
	    # with open(args.data_file, 'a') as f:
	    #   np.savetxt(f, all_values_by_ray, fmt=='%14.6E', newline='\n'))
                # np.savetxt(args.data_file, all_depth, '%14.6E', delimiter=' ', newline='\n', header=depth_hdr, footer=footer) 

    # NOTES
    # All data is stored in memory 	
    print 'Saving all_value array, shaped', all_values.shape, 'to', args.data_file
    # np.savetxt(f_handle, all_values, delimiter='', fmt='%s')
    ###############################################################
    print 'all_values shape', all_values.shape 
    for_statistics = all_values[...,3:]
    print 'for_statistics shape', for_statistics.shape
    ave_all_values = np.average(for_statistics.shape,0)

    # f_handle.close()
    ###############################################################

if __name__ == '__main__':
    main()
