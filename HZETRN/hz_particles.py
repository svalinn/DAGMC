import results_data as rd
import argparse

def parse_arguments():
    """ Argument parsing
    """
    
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', action='store', dest='infile',
	help='The name of the .dat file with one line per direction.')

    args = parser.parse_args()
    if not args.infile:
        raise Exception('Input file not specified. [-f] not set.')
    return args

def main():
    args = parse_arguments()
    names = rd.list_particle_names(args.infile)
    print "Particle names are:", names

if __name__ == '__main__':
    main()
