# jcz fludag.make after MOAB/lib/moab.make
# jcz to be included in the fludag/src makefile


# MOAB_CXXFLAGS =  -Wall -pipe -pedantic -Wno-long-long -O2 -DNDEBUG 
# MOAB_CFLAGS =  -Wall -pipe -pedantic -Wno-long-long -O2 -DNDEBUG 
# MOAB_LDFLAGS =        -L/filespace/people/z/zachman/dagmc_bld/HDF5/lib -L/filespace/people/z/zachman/dagmc_bld/CGM/lib  -L/filespace/people/z/zachman/dagmc_bld/bin



# MOAB_LIBS_LINK = ${MOAB_LDFLAGS} -L${MOAB_LIBDIR} -lMOAB  -lm  -lnetcdf -lhdf5    -lcgm  -L/filespace/people/z/zachman/dagmc_bld/bin -lcubiti19

DAG_LIBDIR     =   /filespace/people/z/zachman/dagmc_bld/MOAB/lib
MOAB_LDFLAGS   = -L/filespace/people/z/zachman/dagmc_bld/bin
MOAB_LIBS_LINK = -L${DAG_LIBDIR} -lMOAB  -lm  -lnetcdf  -ldagmc -L/filespace/people/z/zachman/dagmc_bld/bin
FLUKA_LIBS_LINK = -L${FLUPRO} -lflukahp
# jcz added for fortran compilation
FLDLIBS        = -lgfortran -lstdc++

# MOAB_CXX = g++
# MOAB_CC  = gcc
# jcz added for fortran compilation
FXX = gfortran

# Override MOAB_LIBDIR and MOAB_INCLUDES from above with the correct
# values for the installed MOAB.

MOAB_LIBDIR=/filespace/people/z/zachman/dagmc_bld/MOAB/lib
MOAB_INCLUDES=-I/filespace/people/z/zachman/dagmc_bld/MOAB/include
