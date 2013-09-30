#include /home/patrick/moab_4.6.0/MOAB/lib/moab.make
include /filespace/people/s/shriwise/make_watertight_test_suite/moabs/wo_cgm/4.6.0/moab/lib/moab.make

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

check_watertight_func.o: ../check_watertight_func.cpp ../check_watertight_func.hpp
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c ../check_watertight_func.cpp

test: test_cyl.o gen.o arc.o zip.o cleanup.o check_watertight_func.o
	$(CC) $(LD_FLAGS) -o test_cyl test_cyl.o gen.o arc.o zip.o cleanup.o check_watertight_func.o  \
	${MOAB_LIBS_LINK} -ldagmc


clean:
	rm -f make_watertight.o make_watertight gen.o arc.o zip.o \
	cleanup.o post_process.o post_process check_watertight_func.o mw_fix mw_fix.o test_cyl test_cyl.o
