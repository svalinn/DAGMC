#!/bin/bash
#
# References
# http://stackoverflow.com/questions/16483119/example-of-how-to-use-getopts-in-bash
# http://www.theunixschool.com/2012/08/getopts-how-to-pass-command-line-options-shell-script-Linux.html
# http://www.unix.com/man-page/all/1/getopts/
# http://stackoverflow.com/questions/7529856/bash-getopts-retrieving-multiple-variables-from-one-flag
# 
usage() { echo "Usage: $0 [-f <infile>] [-r <response>] [-p <particle>] [-g <group>] [-d] [-o <outfile>]" 1>&2; exit 1; }

# d is optional
while getopts ":f:r:p:g:do:" clp; do
    case "${clp}" in 
        f)
	   f=${OPTARG}
	   ;;
	r)
	   resp+=("${OPTARG}")
	   ;;
	p)
	   ple+=("${OPTARG}")
	   ;;
	g)
	   g=${OPTARG}
	   ;;
	d)
	   d=1
	   ;;
	o)
	   o=${OPTARG}
	   ;;
	*)
	   usage
	   ;;
    esac
done

shift $((OPTIND-1))

if [ -z "${f}" ] || [ -z "${resp}" ] || [ -z "${ple}" ] || [ -z "${o}" ]; then
    usage
fi

echo "Input data file f = ${f}"
echo "Response r = ${resp[@]}"
echo "Or:"
for val in "${resp[@]}"; do 
    echo "  - $val"
done
echo "Particle p = ${ple[@]}"
if [ ! -z "$d" ]; then
   echo "Write direction d on"
fi
echo "Output file o = ${o}"
