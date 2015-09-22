#!/bin/bash

# Simple bash script with one parameter, the name of the file.  
# It tests the dat file by doing the following
# a) gets the number of lines in the file
# b) for each line in the file, prints out the line number and the number of fields in each line
# usage() { echo "Usage: $0 [-f <infile>] [-r <response>] [-p <particle>] [-g <group>] [-d] [-o <outfile>]" 1>&2; exit 1; }

# http://stackoverflow.com/questions/6022384/bash-tool-to-get-nth-line-from-a-file
# Use double quotes to parse the command line.  Need the extra quotes
lines=`wc -l $1 | cut -f1 -d' '`

# words=`wc -w $2 | cut -f1 -d' '`
for i in `seq 1 "$lines"`;
do
    num=`sed ""$i"q;d" $1 | wc -L`
    echo $i $num 
done
