include /opt/dagmc/moab/lib/moab.make
MOAB_CXXFLAGS =  -Wall -pipe -pedantic -Wno-long-long 
MOAB_CFLAGS = -Wall -pipe -pedantic -Wno-long-long
CXXFLAGS += ${MOAB_CXXFLAGS} -g 
CC = g++
LD_FLAGS = -g
CPPFLAGS += ${MOAB_INCLUDES} -g
CFLAGS   += ${MOAB_CFLAGS} -g
# add -g -pg to both CXX and LD flags to profile

all: uhu

uhu.o: uhu.cpp 
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c uhu.cpp
gen.o: gen.cpp 
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c gen.cpp
geometry.o: geometry.cpp 
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c geometry.cpp

uhu: uhu.o gen.o geometry.o
	$(CC) $(LD_FLAGS) -o uhu uhu.o gen.o geometry.o \
	${MOAB_LIBS_LINK} -ldagmc -lMOAB

clean:
	rm -f uhu *.o 
