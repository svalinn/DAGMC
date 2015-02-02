import subprocess
import argparse
import os
from itertools import islice
import glob
import numpy as np
import csv

# Global names to be read from a separate text file
dose_filename  = ''
doseq_filename = ''
dose_ple_filename  = ''
doseq_ple_filename = ''
diflet_filename = '' 
intlet_filename = ''
flux_filename = ''
intflux_filename = ''

class Results:

    directory=''

    def __init__(self, dir):
        self.directory = dir

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
        vals = np.array([float(words[1]), float(words[2]), float(words[3]), float(words[4]), float(words[5]), float(words[6])])
	return vals

class n700x2(Results):

    def index_from(self, name, index):
        # vals = np.array([[]])
	vals = np.empty((1,700), float)
	lines = Results.last_lines(self, name, 700)
	col = 0
	for line in lines:
	    words = line.split()
	    num = float(words[index])
	    vals[0,col] = num
	    col = col + 1
	return vals
      
    def last_lines(self, name):
        vals = self.index_from(name, 1)
	return vals

    """
    def last_lines(self, name):
	vals = np.array([])
        result_lines = Results.last_lines(self, name, 700)
	for line in result_lines:
	    words = line.split()
	    vals = np.hstack((vals, float(words[1])))
	return vals
   """

class n100x7(Results):

    def last_lines(self, name):
	vals = np.array([])
	lines = Results.last_lines(self, name, 100)
	for line in lines:
	    words = line.split()
	    # We are skipping the first value, which is a flux bin
            nums = np.array([float(words[1]), float(words[2]), float(words[3]), float(words[4]), float(words[5]), float(words[6])])
	    vals = np.vstack((vals, nums))
	return vals
      
    def index_from(self, name, index):
	vals = np.empty((1,100), float)
	lines = Results.last_lines(self, name, 100)
	col = 0
	for line in lines:
	    words = line.split()
	    num = float(words[index])
	    vals[0,col] = num
	    col = col + 1
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

    return 

def depth_header():
    "Return a header for the depth data"
    header =  '\n' + \
              '{0: >12}'.format('X') + \
              '{0: >15}'.format('Y') + \
              '{0: >15}'.format('Z') + \
              '{0: >15}'.format('Dose_All') + \
              '{0: >15}'.format('Doseq_All') + \
              '{0: >15}'.format('dose-neutron') + \
              '{0: >15}'.format('dose-proton') + \
              '{0: >15}'.format('dose-deut.') + \
              '{0: >15}'.format('dose-trit.') + \
              '{0: >15}'.format('dose-H3') + \
              '{0: >15}'.format('dose-H4') + \
              '{0: >15}'.format('doseq-neutron') + \
              '{0: >15}'.format('doseq-proton') + \
              '{0: >15}'.format('doseq-deut.') + \
              '{0: >15}'.format('doseq-trit.') + \
              '{0: >15}'.format('doseq-H3') + \
              '{0: >15}'.format('doseq-H4') + '\n'
    return header

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
    args = parser.parse_args()
    if not args.run_dir:
        raise Exception('Run directory not specified. [-d] not set')

    return args

def main():
    global dose_filename, doseq_filename, dose_ple_filename, doseq_ple_filename
    global diflet_filename, intlet_filename
    global flux_filename, intflux_filename

    # Setup: parse the the command line parameters
    args = parsing()

    # Set up the names of the data files without them being in this file
    get_names("names.txt")
    # print dose_filename, doseq_filename, dose_ple_filename, doseq_ple_filename
    
    # Original Working Directory - will be used for the results file
    owd = os.path.dirname(os.path.abspath(__file__)) + '/'

    # Directory containing the ray subdirectories
    run_path = args.run_dir + '/data/'

    os.chdir(run_path)
    depth_data_reader = nx2('')
    depth_particle_reader = nx7('')

    dif_int_block_reader  = n700x2('')
    flux_block_reader     = n100x7('')

    all_depth = np.array([])
    # dif_block700 = np.array([1,703])
    dif_block700 = np.array([])
    int_block700 = np.array([])
    block700_header = np.array([])
    block100_header = np.array([])
    btest           = np.array([])
    # This could be up to 10,000 subdirectories
    ray_subdirs = glob.glob('*')
    for r in ray_subdirs:
	dose_filepath  = r + '/' + dose_filename 
	doseq_filepath = r + '/' + doseq_filename
	dose_part_filepath  = r + '/' + dose_ple_filename 
	doseq_part_filepath = r + '/' + doseq_ple_filename

	dif_filepath     = r + '/' + diflet_filename
	int_filepath     = r + '/' + intlet_filename
	    
	flux_filepath     = r + '/' + flux_filename
	intflux_filepath     = r + '/' + intflux_filename

        ray = r.split('_') 
	ray_vec = np.array([float(ray[0]), float(ray[1]), float(ray[2])])

	if os.path.exists(dose_filepath) and os.path.exists(doseq_filepath) and \
	   os.path.exists(dose_part_filepath) and os.path.exists(doseq_part_filepath):
            # Add dose and doseq at depth for 1 or 6 particles
	    hnums = depth_data_reader.last_lines(dose_filepath)
	    inums = depth_data_reader.last_lines(doseq_filepath)
	    jnums = depth_particle_reader.last_lines(dose_part_filepath)
	    knums = depth_particle_reader.last_lines(doseq_part_filepath)
	   
	    one_line = np.hstack((ray_vec,hnums,inums,jnums,knums))
	    # First time through...
	    if all_depth.size == 0:
	        all_depth = one_line
	    # ...otherwise append ('vertically stack')
	    else:
	        all_depth = np.vstack((all_depth,one_line)) 
    
	if os.path.exists(dif_filepath) and os.path.exists(int_filepath):
	    # Only need to read this once
	    ray_vec.shape = (1,3)
	    if block700_header.size == 0:
	        fnums = dif_int_block_reader.index_from(dif_filepath, 0)
		pad = np.array([[0.0, 0.0, 0.0]])
		block700_header = np.hstack((pad,fnums))

	    hnums = dif_int_block_reader.last_lines(dif_filepath)
	    dif_line = np.hstack((ray_vec,hnums))

	    inums = dif_int_block_reader.last_lines(int_filepath)
	    int_line = np.hstack((ray_vec,inums))

	    if dif_block700.size == 0:
	        dif_block700 = dif_line
		int_block700 = int_line
	    else:
	        dif_block700 = np.vstack((dif_block700,dif_line))
	        int_block700 = np.vstack((int_block700,int_line))

	if os.path.exists(flux_filepath) and os.path.exists(intflux_filepath):
	    ray_vec.shape = (1,3)
	    if block100_header.size == 0:

	        enums = flux_block_reader.index_from(flux_filepath, 0)
	        fnums = flux_block_reader.index_from(intflux_filepath, 0)
		pad = np.array([[0.0, 0.0, 0.0]])
		block100_header = np.hstack((pad,fnums))
		btest           = np.hstack((pad,enums))
	    
    ###############################################################
    ave_at_depth = np.average(all_depth,0)
    std_at_depth = np.std(all_depth,0)

    dashes15  = '{0: >15}'.format('--') + '{0: >15}'.format('--')
    ave_words = '{0: >12}'.format('--') + dashes15
    std_words = '{0: >12}'.format('--') + dashes15

    for i in range(3,17):
        ave_words += "{:15.6E}".format(ave_at_depth[i])
        std_words += "{:15.6E}".format(std_at_depth[i])

    depth_hdr = depth_header()
    footer = depth_hdr + '\n' + \
            '{0: <12}'.format('Averages')  + '\n' + ave_words + '\n' + \
	    '{0: <12}'.format('Std.-Dev.') + '\n' + std_words + '\n' 

    ###############################################################
    os.chdir(owd)
    np.savetxt('tmp1.dat', all_depth, '%14.6E', delimiter=' ', newline='\n', header=depth_hdr, footer=footer) 

    ################################################################
    ave_diflet = np.average(dif_block700, axis=0)
    std_diflet = np.std(dif_block700,axis=0)
    stat_diflet = np.vstack((ave_diflet, std_diflet))
    # stat_diflet[0:2, 0:3] = 0

    ave_intlet   = np.average(int_block700,0)
    std_intlet   = np.std(int_block700,0)
    stat_intlet = np.vstack((ave_intlet, std_intlet))
    # stat_intlet[0:2, 0:3] = 0
    ################################################################
    # Go into append mode
    f_handle = file('tmp1.dat', 'a')
    np.savetxt(f_handle, block700_header, '%14.6E', delimiter=' ', newline='\n', header='DIF_LET\n')
    np.savetxt(f_handle, dif_block700,    '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, stat_diflet,     '%14.6E', delimiter=' ', newline='\n', header='DIF-AVERAGES-STD.DEV\n', footer='\n')
    # get rows 0 and 1
    # np.savetxt(f_handle, stat_diflet[0,:],  '%14.6E', delimiter=' ', newline='\n', header='DIF-AVERAGES-STD.DEV\n', footer='\n')
    # np.savetxt(f_handle, stat_diflet[1,:],  '%14.6E', delimiter=' ', newline='\n', header='DIF-AVERAGES-STD.DEV\n', footer='\n')
    np.savetxt(f_handle, block700_header, '%14.6E', delimiter=' ', newline='\n', header='INT_LET\n')
    np.savetxt(f_handle, int_block700,    '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, stat_intlet,     '%14.6E', delimiter=' ', newline='\n', header='INT-AVERAGES-STD.DEV\n', footer='\n')


if __name__ == '__main__':
    main()
