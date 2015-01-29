import subprocess
import argparse
import os
from itertools import islice
import glob
import numpy as np
import csv

dose_filebase  = 'ascii_dose_table'

doseq_filebase = 'ascii_doseq_table'

class Results:

    directory=''

    def __init__(self, dir):
        self.directory = dir

    def process(self, name, n):
	filepath = self.directory + name
        return check_last_n_lines(filepath, n)

class nx2(Results):

    def process(self, name, n):
        result_line = Results.process(self, name, n)
	cols = result_line[0].split()
	print cols
        return cols


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

    # Setup: parse the the command line parameters
    args = parsing()
    curdir = os.path.dirname(os.path.abspath(__file__)) + '/'
    run_path = args.run_dir + '/data/'

    os.chdir(run_path)
    g = nx2('')

    ray_subdirs = glob.glob('*')
    doses = []
    doseqs = []
    header = ",X Y Z Dose-Depth Dose Doseq-Depth Doseq"
    # Open as write, then append, in order to clear previous
    with open('../results_data.csv', 'w') as s:
	s.write(header + '\n')
        s.close
    with open('../results_data.csv', 'a') as s:
        for r in ray_subdirs:
            ray = r.split('_') 
	    dosefile  = r + '/' + dose_filebase + '.dat'
	    doseqfile = r + '/' + doseq_filebase + '.dat'
	    # print 'ray = ', ray
	    if os.path.exists(dosefile) and os.path.exists(doseqfile):
		line = "," + ray[0] + "," +  ray[1] + "," + ray[2] + ","
	        cols = g.process(dosefile, 1)
		line += cols[0] + "," + cols[1] + ","
	        doses.append(float(cols[len(cols)-1]))
	        # print ray, cols

	        cols = g.process(doseqfile,1)
	        doseqs.append(float(cols[len(cols)-1]))
		line += cols[0] + "," + cols[1] + "\n"
		s.write(line)

        print 'doses', doses
	dose_ave = np.average(doses)
	dose_std = np.std(doses)
        print 'Average dose', "{:1.6E}".format(dose_ave)
        print 'Std. Dev. dose', "{:1.6E}".format(dose_std)
    

        print 'doseqs', doseqs
	doseq_ave = np.average(doseqs)
	doseq_std = np.std(doseqs)
        ave_doseq_str =  "{:1.6E}".format(doseq_ave)
        std_doseq_str = "{:1.6E}".format(doseq_std)
	s.write("Averages,,,,," +  "{:1.6E}".format(dose_ave) + ",," +  ave_doseq_str + "\n") 
	#s.write('Std Dev', '', '', '', dose_std, '', doseq_std) 


if __name__ == '__main__':
    main()
