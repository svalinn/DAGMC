import argparse
import os
import sys
import numpy as np
import fileinput

"""
    Convenience function to create list of column header names
"""
def construct_depth_header(species):
    # Number of values in the header that are not data results
    "Return a header for the depth data."
    header_list = ['{0: >12}'.format('X'), 
              '{0: >12}'.format('Y'),
              '{0: >12}'.format('Z'),
	      ###########################
	      # Add Depth Column
	      '{0: >12}'.format('Depth'),
	      ###########################
              '{0: >12}'.format('Dose_All'),
              '{0: >12}'.format('Doseq_All') ]
    if species == 6:
        header_list +=  ['{0: >12}'.format('dose_neutron'),
              '{0: >12}'.format('dose_proton'),
              '{0: >12}'.format('dose_deut'),
              '{0: >12}'.format('dose_trit'),
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

def get_slabs_to_ref(filename):
    """ Read an integer from a text file.
    """
    with open(filename) as f:
        return int(f.readline())

"""
Argument parsing
returns : args: -d for the run directory
ref     : DAGMC/tools/parse_materials/dagmc_get_materials.py
"""
def parse_command_line_arguments():
    global data_subdir_name

    parser = argparse.ArgumentParser()

    parser.add_argument( '-d', '--run_dir', 
        help='The name of the holding directory for all hzetrn runs')
    parser.set_defaults(run_dir='rundir')

    parser.add_argument(
        '-s', '--data_dir', 
	help='The name of the data subdirectory holding the direction subdirs')
    parser.set_defaults(data_dir='data')

    parser.add_argument( '-o', '--out_file', 
	help='The relative path to the output file')
    parser.set_defaults(out_file='out.dat')

    args = parser.parse_args()
    
    return args

def main():
    args = parse_command_line_arguments()

    run_path  = os.path.join(os.getcwd(), args.run_dir)
    data_path = os.path.join(run_path, args.data_dir)
    out_path  = os.path.join(os.getcwd(), args.out_file)

    # Start header section; depends on species, data_path
    rad_env = ''
    rad_env_path = os.path.join(data_path, 'rad_env.txt')
    with open(rad_env_path) as f:
        rad_env = f.readline()
    species = 0
    if rad_env == 'spe':
       species = 6
    elif rad_env == 'gcr':
       species = 59
    else:
        raise Exception('Unknown radiation environment {} in {}'.format(rad_env, rad_env_filepath))
    depth_header = construct_depth_header(species)
    index_path = os.path.join(data_path, 'index_header')
    # No need to convert these formatted values to floats, they're just getting written to a file
    with open(index_path) as f:
        index_header_line = f.readline()
    # The index header was already formatted so append it after a space
    with open(out_path, 'w') as f:
        np.savetxt(f, depth_header, fmt='%s', newline=' ' + index_header_line)
    # End of header section: header has been written    
    ###############################################

    
    # Start of data section
    index = 1
    row_path = os.path.join(data_path, str(index) + '/row')
    all_values = np.array([])
    with open(out_path, 'a') as fout:
        while(os.path.isfile(row_path)):
	    with open(row_path) as fin:
	        line = fin.readline()
	        fout.write(line)

	        # Will need the numbers for stats
	        row = np.array([map(float, line.split())])
	        print 'np row shape', row.shape
		# Stack up the rows, but drop the first four columns
                if all_values.size == 0:
	            all_values = row[:,4:]
		    print 'all_values shape', all_values.shape
	        else:
	            all_values = np.vstack((all_values, row[:,4:]))
		    print index, row[0,0:4]
	    index = index + 1
            row_path = os.path.join(data_path, str(index) + '/row')
    # End of data section
    #######################
   
    # Start of statistics section
    print 'all_values shape', all_values.shape
    
    # Averages, Std. Dev.
    ave_all_values = np.average(all_values, axis=0)
    std_all_values = np.std(all_values, axis=0)

    with open(out_path, 'a') as f:
        np.savetxt(f, [['Average', 'over', 'direction', 'at_Depth']], fmt=['%12s', '%12s', '%12s', '%12s'], newline=' ')
        np.savetxt(f, [ave_all_values], fmt='%12.6E', newline='\n')
        np.savetxt(f, [['Std.Dev.', 'over', 'direction', 'at_Depth']], fmt=['%12s', '%12s', '%12s', '%12s'], newline=' ')
        np.savetxt(f, [std_all_values], fmt='%12.6E', newline='\n')
      
if __name__ == '__main__':
    main()
