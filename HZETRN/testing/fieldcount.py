import os
import argparse
import glob
import read_utils as ru


def parsing():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', action='store', dest='filename',
        help='Name of file to read.')
	
    parser.add_argument(
        '-n', action='store', dest='last', type=int,
        help='Lines relative to LAST. "-1" means first line for --all == False')
    parser.set_defaults(last=1)

    parser.add_argument(
        '--all', action='store_true', dest='get_all',
	help='Perform count for all files starting with "ascii" in named directory.')
    parser.set_defaults(get_all=False)

    args = parser.parse_args()
    
    if not args.filename:
       raise Exception('File not specified. [-f] not set.')

    return args


def main():
    
    args = parsing()

    location = args.filename
    n = args.last
    print 'filename', location
    print 'from end', n

    owd = os.path.dirname(os.path.abspath(__file__)) + '/'
    print 'OWD: ', owd
    
    if args.get_all:
        os.chdir(location)
	print 'Line from end: ', n
	listing = glob.glob('ascii*')
	print 'Files are', listing
        for file in listing: 
	    line = ru.line_from_end(file, n)
	    words = line.split()
            print '{0: <30}'.format(file), 'fields =', len(words)
            # print 'file', file, 'all or first 10', words[0:9] if len(words) >= 10 else words
    else:
	if n == -1:
	    line = ru.first_line(location)
	else:
            line = ru.line_from_end(location, n)
        words = line.split()

        print 'Number of fields on ', n, 'th last line is', len(words)
        print words[0:9]
    return

if __name__ == '__main__':
    main()
