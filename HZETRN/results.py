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
neutron_flux_filename = ''

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

class n100x7(Results):

    def last_lines(self, name):
	vals = np.array([])
	nums = np.array([])
	lines = Results.last_lines(self, name, 100)
	for line in lines:
	    words = line.split()
	    # We are skipping the first value, which is a flux bin
            nums = np.array([float(words[1]), float(words[2]), float(words[3]), float(words[4]), float(words[5]), float(words[6])])
	    if vals.size == 0:
	        vals = nums
	    else:
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

class n100x4(Results):

    def last_lines(self, name):
	vals = np.array([])
	nums = np.array([])
	lines = Results.last_lines(self, name, 100)
	for line in lines:
	    words = line.split()
	    # We are skipping the first value, which is a flux bin
            nums = np.array([float(words[1]), float(words[2]), float(words[3])] )
	    if vals.size == 0:
	        vals = nums
	    else:
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
              '{0: >15}'.format('dose-He3') + \
              '{0: >15}'.format('dose-He4') + \
              '{0: >15}'.format('doseq-neutron') + \
              '{0: >15}'.format('doseq-proton') + \
              '{0: >15}'.format('doseq-deut.') + \
              '{0: >15}'.format('doseq-trit.') + \
              '{0: >15}'.format('doseq-He3') + \
              '{0: >15}'.format('doseq-He4') + '\n'
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
    global neutron_flux_filename 

    # Setup: parse the the command line parameters
    args = parsing()

    # Set up the names of the data files without them being in this file
    get_names("names.txt")
    
    # Original Working Directory - will be used for the results file
    owd = os.path.dirname(os.path.abspath(__file__)) + '/'

    # Directory containing the ray subdirectories
    run_path = args.run_dir + '/data/'

    os.chdir(run_path)
    depth_data_reader = nx2('')
    depth_particle_reader = nx7('')

    dif_int_block_reader  = n700x2('')
    flux_block_reader     = n100x7('')
    nflux_block_reader    = n100x4('')

    # Numpy Arrays
    all_depth = np.array([])
   
    dif_block700 = np.array([])
    int_block700 = np.array([])
    block700_header = np.array([])
    block100_header = np.array([])
    nblock100_header = np.array([])

    flux_block100        = np.array([])
    ray_flux_block100    = np.array([])
    ray_intflux_block100 = np.array([])
    ray_flux_stack       = np.array([])
    ray_intflux_stack    = np.array([])

    nflux_block100     = np.array([])
    ray_nflux_block100 = np.array([])
    ray_nflux_stack    = np.array([])
    
    ray_subdirs = glob.glob('*')
    for r in ray_subdirs:
	dose_filepath  = r + '/' + dose_filename 
	doseq_filepath = r + '/' + doseq_filename
	dose_part_filepath  = r + '/' + dose_ple_filename 
	doseq_part_filepath = r + '/' + doseq_ple_filename

	dif_filepath     = r + '/' + diflet_filename
	int_filepath     = r + '/' + intlet_filename
	    
	flux_filepath     = r + '/' + flux_filename
	intflux_filepath  = r + '/' + intflux_filename

        neutron_flux_filepath = r + '/' + neutron_flux_filename

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
	    if block100_header.size == 0:
	        fnums = flux_block_reader.index_from(flux_filepath, 0)
		pad = np.array([[0.0, 0.0, 0.0]])
		block100_header = np.hstack((pad,fnums))

            hnums = flux_block_reader.last_lines(flux_filepath)
	    inums = flux_block_reader.last_lines(intflux_filepath)

	    # Flip the rows and columns: shape (100, 6) -> (6, 100)
	    rhnums = np.swapaxes(hnums,0,1)
	    rinums = np.swapaxes(inums,0,1)
	    
	    # Force the ray vector shape to be such that it can be used for padding
	    ray_vec.shape = (3)

	    # Pad the first row to initialize 
	    ray_flux_block100    = np.hstack((ray_vec,(rhnums[0,:])))
	    ray_intflux_block100 = np.hstack((ray_vec,(rinums[0,:])))
    	    
	    for i in range(1,6):
	        # pad each row with the ray vector values
	        tmp_row = np.hstack((ray_vec,rhnums[i,:]))
		ray_flux_block100 = np.vstack((ray_flux_block100, tmp_row))

	        tmp_row = np.hstack((ray_vec,rinums[i,:]))
		ray_intflux_block100 = np.vstack((ray_intflux_block100, tmp_row))
	    
	    # Now stack up the blocks like they are planes in a 3d image
            if ray_flux_stack.size == 0:
    	        print 'ray_flux_block100 shape after vstack', ray_flux_block100.shape
	        ray_flux_stack    = ray_flux_block100
	        ray_intflux_stack = ray_intflux_block100
            else:
		print 'ray_flux_stack shape', ray_flux_stack.shape
	        # Stack the block of data along the third dimension (= 2)
	        ray_flux_stack    = np.dstack((ray_flux_stack,   ray_flux_block100))
	        ray_intflux_stack = np.dstack((ray_intflux_stack,ray_intflux_block100))
	# Neutron Flux: could possibly be joined with rest of fluxes
	if os.path.exists(neutron_flux_filepath):
	    if nblock100_header.size == 0:
	        nfnums = nflux_block_reader.index_from(neutron_flux_filepath, 0)
		pad = np.array([[0.0, 0.0, 0.0]])
		nblock100_header = np.hstack((pad,nfnums))

            jnums = nflux_block_reader.last_lines(neutron_flux_filepath)

	    # Flip the rows and columns: shape (100, 6) -> (6, 100)
	    rjnums = np.swapaxes(jnums,0,1)
	    
	    # Pad the first row to initialize 
	    ray_nflux_block100 = np.hstack((ray_vec,(rjnums[0,:])))
    	    
	    for i in range(1,3):
	        # pad each row with the ray vector values
	        tmp_row = np.hstack((ray_vec,rjnums[i,:]))
		ray_nflux_block100 = np.vstack((ray_nflux_block100, tmp_row))

	    # Now stack up the blocks like they are planes in a 3d image
            if ray_nflux_stack.size == 0:
	        ray_nflux_stack = ray_nflux_block100
            else:
	        # Stack the block of data along the third dimension (= 2)
	        ray_nflux_stack = np.dstack((ray_nflux_stack,  ray_nflux_block100))
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
    np.savetxt('tmp2.csv', all_depth, '%14.6E', delimiter=' ', newline='\n', header=depth_hdr, footer=footer) 

    ################################################################
    ave_diflet  = np.average(dif_block700, axis=0)
    std_diflet  = np.std(dif_block700,axis=0)
    stat_diflet = np.vstack((ave_diflet, std_diflet))

    ave_intlet  = np.average(int_block700,0)
    std_intlet  = np.std(int_block700,0)
    stat_intlet = np.vstack((ave_intlet, std_intlet))

    # Fluxes
    ave_flux  = np.average(ray_flux_stack, axis=2)
    std_flux  = np.std(ray_flux_stack, axis=2)
    stat_flux = np.dstack((ave_flux, std_flux))
    
    ave_intflux  = np.average(ray_intflux_stack, axis=2)
    std_intflux  = np.std(ray_intflux_stack, axis=2)
    stat_intflux = np.dstack((ave_intflux, std_intflux))

    # Neutron Flux
    ave_nflux  = np.average(ray_nflux_stack, axis=2)
    std_nflux  = np.std(ray_nflux_stack, axis=2)
    stat_nflux = np.dstack((ave_nflux, std_nflux))
    ################################################################
    # Go into append mode
    f_handle = file('tmp2.csv', 'a')
    np.savetxt(f_handle, block700_header, '%14.6E', delimiter=' ', newline='\n', header='DIF_LET\n')
    np.savetxt(f_handle, dif_block700,    '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, stat_diflet,     '%14.6E', delimiter=' ', newline='\n', header='DIF-AVERAGES-STD.DEV\n', footer='\n')
    # get rows 0 and 1
    np.savetxt(f_handle, block700_header, '%14.6E', delimiter=' ', newline='\n', header='INT_LET\n')
    np.savetxt(f_handle, int_block700,    '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, stat_intlet,     '%14.6E', delimiter=' ', newline='\n', header='INT-AVERAGES-STD.DEV\n', footer='\n')

    ####################################################################
    # FLUX
    # write out the blocks and block averages for all the flux particles
    ####################################################################
    np.savetxt(f_handle, block100_header, '%14.6E', delimiter=' ', newline='\n', header='PROTON-FLUX' )
    np.savetxt(f_handle, np.swapaxes(ray_flux_stack[0,:,:],0,1), '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, np.swapaxes(stat_flux[0,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='Proton-Average-Std.Dev\n', footer='\n')

    np.savetxt(f_handle, block100_header, '%14.6E', delimiter=' ', newline='\n', header='NEUTRON-FLUX')
    np.savetxt(f_handle, np.swapaxes(ray_flux_stack[1,:,:],0,1), '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, np.swapaxes(stat_flux[1,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='Neutron-Average-Std.Dev\n', footer='\n')
    np.savetxt(f_handle, block100_header, '%14.6E', delimiter=' ', newline='\n', header='FLUX')
    np.savetxt(f_handle, np.swapaxes(ray_flux_stack[2,:,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='DEUTERON-FLUX')
    np.savetxt(f_handle, np.swapaxes(stat_flux[2,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='Deuteron-Average-Std.Dev\n', footer='\n')

    np.savetxt(f_handle, block100_header, '%14.6E', delimiter=' ', newline='\n', header='TRITON-FLUX')
    np.savetxt(f_handle, np.swapaxes(ray_flux_stack[3,:,:],0,1), '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, np.swapaxes(stat_flux[3,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='Triton-Average-Std.Dev\n', footer='\n')
    
    np.savetxt(f_handle, block100_header, '%14.6E', delimiter=' ', newline='\n', header='He3-FLUX')
    np.savetxt(f_handle, np.swapaxes(ray_flux_stack[4,:,:],0,1), '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, np.swapaxes(stat_flux[4,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='He3-Average-Std.Dev\n', footer='\n')
    
    np.savetxt(f_handle, block100_header, '%14.6E', delimiter=' ', newline='\n', header='He4-FLUX')
    np.savetxt(f_handle, np.swapaxes(ray_flux_stack[5,:,:],0,1), '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, np.swapaxes(stat_flux[5,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='He4-Average-Std.Dev\n', footer='\n')
    ####################################################################
    # INTFLUX
    # write out the blocks and block averages for all the flux particles
    ####################################################################
    np.savetxt(f_handle, block100_header, '%14.6E', delimiter=' ', newline='\n', header='PROTON-INTFLUX' )
    np.savetxt(f_handle, np.swapaxes(ray_intflux_stack[0,:,:],0,1), '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, np.swapaxes(stat_intflux[0,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='Proton-Average-Std.Dev\n', footer='\n')

    np.savetxt(f_handle, block100_header, '%14.6E', delimiter=' ', newline='\n', header='NEUTRON-INTFLUX')
    np.savetxt(f_handle, np.swapaxes(ray_intflux_stack[1,:,:],0,1), '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, np.swapaxes(stat_intflux[1,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='Neutron-Average-Std.Dev\n', footer='\n')
    np.savetxt(f_handle, block100_header, '%14.6E', delimiter=' ', newline='\n', header='INTFLUX')
    np.savetxt(f_handle, np.swapaxes(ray_intflux_stack[2,:,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='DEUTERON-INTFLUX')
    np.savetxt(f_handle, np.swapaxes(stat_intflux[2,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='Deuteron-Average-Std.Dev\n', footer='\n')

    np.savetxt(f_handle, block100_header, '%14.6E', delimiter=' ', newline='\n', header='TRITON-INTFLUX')
    np.savetxt(f_handle, np.swapaxes(ray_intflux_stack[3,:,:],0,1), '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, np.swapaxes(stat_intflux[3,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='Triton-Average-Std.Dev\n', footer='\n')
    
    np.savetxt(f_handle, block100_header, '%14.6E', delimiter=' ', newline='\n', header='He3-INTFLUX')
    np.savetxt(f_handle, np.swapaxes(ray_intflux_stack[4,:,:],0,1), '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, np.swapaxes(stat_intflux[4,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='He3-Average-Std.Dev\n', footer='\n')
    
    np.savetxt(f_handle, block100_header, '%14.6E', delimiter=' ', newline='\n', header='He4-INTFLUX')
    np.savetxt(f_handle, np.swapaxes(ray_intflux_stack[5,:,:],0,1), '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, np.swapaxes(stat_intflux[5,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='He4-Average-Std.Dev\n', footer='\n')
    
    np.savetxt(f_handle, nblock100_header, '%14.6E', delimiter=' ', newline='\n', header='FORWARD-NEUTRON-FLUX')
    np.savetxt(f_handle, np.swapaxes(ray_nflux_stack[0,:,:],0,1), '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, np.swapaxes(stat_nflux[0,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='Forward-Average-Std.Dev\n', footer='\n')

    np.savetxt(f_handle, nblock100_header, '%14.6E', delimiter=' ', newline='\n', header='BACKWARD-NEUTRON-FLUX')
    np.savetxt(f_handle, np.swapaxes(ray_nflux_stack[1,:,:],0,1), '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, np.swapaxes(stat_nflux[1,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='Backward-Average-Std.Dev\n', footer='\n')

    np.savetxt(f_handle, nblock100_header, '%14.6E', delimiter=' ', newline='\n', header='TOTAL-NEUTRON-FLUX')
    np.savetxt(f_handle, np.swapaxes(ray_nflux_stack[2,:,:],0,1), '%14.6E', delimiter=' ', newline='\n', footer='\n')
    np.savetxt(f_handle, np.swapaxes(stat_nflux[2,:],0,1), '%14.6E', delimiter=' ', newline='\n', header='Total-Average-Std.Dev\n', footer='\n')

if __name__ == '__main__':
    main()
