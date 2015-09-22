import argparse
import os
import numpy as np

# Number of columns not containing results data
meta = 6

def construct_depth_header(species):
    """ Return a header for the depth data.  This portion of the
    header is ascii titles.  The rest of the header is numeric.

    Parameters
    ----------
    species : integer
        Number of species.  Dependent on a command line argument.

    Returns
    -------
    A numpy row array of ascii words, formatted to sent to the database
    
    """
    header_list = ['{0: >12}'.format('X'), 
              '{0: >12}'.format('Y'),
              '{0: >12}'.format('Z'),
	      ###########################
	      # Add meta columns
	      '{0: >12}'.format('Density'),
	      '{0: >12}'.format('Tot_Depth'),
	      '{0: >12}'.format('Tot_g/cm^2'),
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
    Parameters
    ----------
    filename : string

    Returns
    -------
    An integer read from the file
    """
    with open(filename) as f:
        return int(f.readline())

def write_header(data_path, out_path):
    """ Construct a header for the data and write it to the database file (out_path)

    Parameters
    ----------
    data_path : string
        Path to the data directory

    species : integer
        The number of particles irradiating the geometry

    out_path : string
        Path to the file to write the header to.
    """
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
    # No need to convert these formatted values to floats, since they are being written to a file
    with open(index_path) as f:
        index_header_line = f.readline()
    # The index header was already formatted so append it after a space
    with open(out_path, 'w') as f:
        np.savetxt(f, depth_header, fmt='%s', newline=' ' + index_header_line)
    return

def write_lines_to_db(data_path, out_path):
    """ Collect all the row data previously stored and write it to a 
    database file

    Parameters
    ----------
    data_path : string
        Path to the data directory

    out_path : string
        Path to the file to write the header to.
    
    Returns
    -------
    A two-dimensional numpy array of rows containing all the data (at the depth
    of interest) for each direction.  This array is to be used for statistics,
    so the meta data (direction and depth) is not included.
    """
    index = 1
    row_path = os.path.join(data_path, str(index) + '/row')
    all_values = np.array([])
    with open(out_path, 'a') as fout:
        while(os.path.isfile(row_path)):
	    with open(row_path) as fin:
	        line = fin.readline()
	        fout.write(line)

	        # Done writing, but will need the numbers for stats
	        row = np.array([map(float, line.split())])
		# Stack up the rows, but drop the meta columns
                if all_values.size == 0:
	            all_values = row[:,meta:]
	        else:
	            all_values = np.vstack((all_values, row[:,meta:]))
		print index, row[0,0:meta]
	    index = index + 1
            row_path = os.path.join(data_path, str(index) + '/row')
    return all_values

def parse_command_line_arguments():
    """Perform command line argument parsing
    Notes
    -----
    - All arguments have a default
    - A run directory must have been previously set up with 
      HZETRN directories and whatever cross-sections are needed
    - action = 'store' is default, so removed

    Returns
    -------
    An object whose attributes are the parameter names
    """

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

    write_header(data_path, out_path)
    all_values = write_lines_to_db(data_path, out_path)
    
    # Averages, Std. Dev.
    ave_all_values = np.average(all_values, axis=0)
    std_all_values = np.std(all_values, axis=0)

    # stat_words = ['over', 'direction', 'at', 'Reference',  'Depth']
    # stat_words_fmt = [['%12s']*meta]
    # all_stat_words = stat_words.insert(0,'Average')
    with open(out_path, 'a') as f:
        np.savetxt(f, [['Average', 'over', 'direction', 'at', 'Reference', 'Depth']], \
	               fmt=['%12s', '%12s', '%12s', '%12s', '%12s', '%12s'], newline=' ')
        # np.savetxt(f, [['Average'] + stat_words], fmt=stat_words_fmt, newline=' ')
        np.savetxt(f, [ave_all_values], fmt='%12.6E', newline='\n')
        np.savetxt(f, [['Std.Dev.', 'over', 'direction', 'at', 'Reference', 'Depth']], \
	               fmt=['%12s', '%12s', '%12s', '%12s', '%12s', '%12s'], newline=' ')
        # np.savetxt(f, [['Std.Dev.'] + stat_words], fmt=stat_words_fmt, newline=' ')
        np.savetxt(f, [std_all_values], fmt='%12.6E', newline='\n')
      
if __name__ == '__main__':
    main()
