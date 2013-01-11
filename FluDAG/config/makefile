########################################################################
#
# Julie Zachman
# December 14, 2012
# Rudimentary makefile to test compiling and linking of heterogenous sources 
# 

MOAB_INCLUDES  = -I/filespace/people/z/zachman/dagmc_bld/MOAB/include
MOAB_LIBDIR    =   /filespace/people/z/zachman/dagmc_bld/MOAB/lib
MOAB_LIBS_LINK = -L${MOAB_LIBDIR} -lMOAB  -lm  -lnetcdf  -ldagmc -L/filespace/people/z/zachman/dagmc_bld/bin

FSRCS = 

CCSRCS = 

OBJS = $(FSRCS:.f=.o) $(CCSRCS:.c=.o)
PRGS = $(FSRCS:.f=) $(CCSRCS:.c=)

.SUFFIXES:
.SUFFIXES: .f .c .o

.f.o:
	gfortran -c $<

.c.o: 
	g++ -c $<

CPPFLAGS = ${MOAB_INCLUDES}
CXXFLAGS = -g

default: all

all: mainDag mainf

mainDag: mainDag.o obj/fludagW.o
	$(CXX) -o $@ $< $(MOAB_LIBS_LINK) $(FLUKA_LIBS_LINK) -lgfortran
        
mainDag.o : mainDag.cpp
	g++ -c $< $(CPPFLAGS) 

fludagW.o : cpp/fludagW.cpp
	cd cpp
	g++ -c $(CPPFLAGS) $<
	cd ..
	mv -f fludagW.o obj
 
# Create a fortrain main linking object files compiled from fortran, C and C++
mainf: obj/main_f.o obj/simple.o obj/test.o obj/fludagW.o
	gfortran -o  $@ obj/main_f.o obj/simple.o obj/test.o obj/fludagW.o $(MOAB_LIBS_LINK) $(FLUKA_LIBS_LINK) -lstdc++

main_f.o: for/main_f.f
	cd for
	gfortran -c $<
	cd ..
	mv -f $@ obj

simple.o: for/simple.f 
	cd for
	gfortran -c $<
	cd ..
	mv -f simple.o obj

test.o: c/test.c
	cd c
	gcc -c $<
	cd ..; ls
	mv -f $@ obj

clean:
	cp *.o object
	rm *.o

print-%:
	@echo $* = $($*)
