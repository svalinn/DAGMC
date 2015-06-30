#!/bin/bash

usage() { echo "Usage: $0 [-f <infile>] [-r <response>] [-p <particle>] [-g <group>] [-d] [-o <outfile>]" 1>&2; exit 1; }

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

echo "f = ${f}"
echo "r = ${resp[@]}"
echo "Or:"
for val in "${resp[@]}"; do 
    echo "  - $val"
done
echo "p = ${ple[@]}"
if [ ! -z "$d" ]; then
   echo "d on"
fi
echo "o = ${o}"
	   
