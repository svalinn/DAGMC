include /home/patrick/moab_4.6.0/MOAB/lib/moab.make
#include /filespace/people/s/shriwise/make_watertight_test_suite/moabs/wo_cgm/4.6.0/moab/lib/moab.make

MOAB_CXXFLAGS =  -Wall -pipe -pedantic -Wno-long-long 
MOAB_CFLAGS = -Wall -pipe -pedantic -Wno-long-long
CXXFLAGS += ${MOAB_CXXFLAGS} -g 
CC = g++
LD_FLAGS = -g
CPPFLAGS += ${MOAB_INCLUDES} -g
CFLAGS   += ${MOAB_CFLAGS} -g
# add -g -pg to both CXX and LD flags to profile

all: make_watertight post_process check_watertight fix

gen.o: gen.cpp gen.hpp
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c gen.cpp

arc.o: arc.cpp arc.hpp
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c arc.cpp

zip.o: zip.cpp zip.hpp
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c zip.cpp

cleanup.o: cleanup.cpp cleanup.hpp
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c cleanup.cpp

make_watertight: make_watertight.o gen.o arc.o zip.o cleanup.o
	$(CC) $(LD_FLAGS) -o make_watertight make_watertight.o gen.o \
	arc.o zip.o cleanup.o ${MOAB_LIBS_LINK} 

post_process: post_process.o gen.o arc.o zip.o cleanup.o
	$(CC) $(LD_FLAGS) -o post_process post_process.o gen.o \
	arc.o zip.o cleanup.o ${MOAB_LIBS_LINK} 

check_watertight: check_watertight.o gen.o cleanup.o
	$(CC) $(LD_FLAGS) -o check_watertight check_watertight.o gen.o \
	arc.o zip.o cleanup.o ${MOAB_LIBS_LINK} 

fix: mw_fix.o gen.o arc.o zip.o cleanup.o
	$(CC) $(LD_FLAGS) -o mw_fix mw_fix.o gen.o arc.o zip.o cleanup.o  \
	${MOAB_LIBS_LINK} -ldagmc

clean:
	rm -f make_watertight.o make_watertight gen.o arc.o zip.o \
	cleanup.o post_process.o post_process check_watertight.o check_watertight mw_fix mw_fix.o
