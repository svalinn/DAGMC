#!/bin/bash

usage() { echo "Usage: $0 [-f <infile>] [-r <response>] [-p <particle>] [-g <group>] [-d] [-o <outfile>]" 1>&2; exit 1; }

while getopts ":f:r:p:g:o:" clp; do
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
echo "p = ${ple[@]}"
if [ ! -z "${g}" ]; then
   OIFS=$IFS
   IFS='-'
   arr="${g}"
   echo "g ="
   for x in $arr
   do
      echo "  $x"
   done
   IFS=$OIFS
fi

echo "o = ${o}"
	   
