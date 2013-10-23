#include /home/patrick/moab_4.6.0/MOAB/lib/moab.make
include /filespace/people/s/shriwise/dagmc_bld/MOAB/lib/moab.make

#INC = -I/home/patrick/scratch/moab_tools/make_watertight
INC = -I/filespace/people/s/shriwise/scratch/moab_tools/make_watertight
MOAB_CXXFLAGS =  -Wall -pipe -pedantic -Wno-long-long ${INC}
MOAB_CFLAGS = -Wall -pipe -pedantic -Wno-long-long ${INC}
CXXFLAGS += ${MOAB_CXXFLAGS} -g 
CC = g++
LD_FLAGS = -g
CPPFLAGS += ${MOAB_INCLUDES} -g
CFLAGS   += ${MOAB_CFLAGS} -g
# add -g -pg to both CXX and LD flags to profile

all: test

gen.o: ../gen.cpp ../gen.hpp
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c ../gen.cpp

arc.o: ../arc.cpp ../arc.hpp
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c ../arc.cpp

zip.o: ../zip.cpp ../zip.hpp
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c ../zip.cpp

cleanup.o: ../cleanup.cpp ../cleanup.hpp
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c ../cleanup.cpp

cw_func.o: ../cw_func.cpp ../cw_func.hpp
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c ../cw_func.cpp

mw_func.o: ../mw_func.cpp ../mw_func.hpp
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c ../mw_func.cpp

test: test_cyl.o gen.o arc.o zip.o cleanup.o cw_func.o mw_func.o
	$(CC) $(LD_FLAGS) -o test_cyl test_cyl.o gen.o arc.o zip.o cleanup.o cw_func.o  \
	mw_func.o ${MOAB_LIBS_LINK} -ldagmc


clean:
	rm -f make_watertight.o make_watertight gen.o arc.o zip.o \
	cleanup.o post_process.o post_process cw_func.o mw_fix mw_fix.o test_cyl test_cyl.o mw_func.o
