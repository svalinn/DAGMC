# Makefile for hzetrn_paths
#
# Needs to be called while setting the MOAB_BASE_DIR variable.
# MOAB_BASE_DIR should be the directory in which there is an
# 'include' and 'lib' directory with the MOAB include files and
# libraries.
#
# prompt%> make MOAB_BASE_DIR=/path/to/installed/cgm
#
# MOAB_BASE_DIR = /path/to/installed/cgm

include ${MOAB_BASE_DIR}/lib/moab.make


CXXSOURCES = hzetrn_paths.cc 
CXXOBJS = hzetrn_paths.o 

hzetrn_paths: ${CXXOBJS} Makefile
	${CXX} ${CXXFLAGS} -g -o $@ ${CXXOBJS} ${MOAB_LIBS_LINK}

.cc.o:
	${CXX} ${CXXFLAGS} ${MOAB_CXXFLAGS} ${MOAB_INCLUDES} -g -o $@ -c $<

clean:
	rm -rf hzetrn_paths *.o
