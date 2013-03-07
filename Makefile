include /opt/dagmc/moab/lib/moab.make
CXXFLAGS += ${MOAB_CXXFLAGS} -g -pg
CC = g++
LD_FLAGS = -pg 
CPPFLAGS += ${MOAB_INCLUDES} -pg
CFLAGS   += ${MOAB_CFLAGS} -pg
# add -g -pg to both CXX and LD flags to profile

all: print_vols

print_vols: print_vols.o
	$(CC) $(LD_FLAGS) -o print_vols print_vols.o ${MOAB_LIBS_LINK} -ldagmc

print_vols.o: print_vols.cpp
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c print_vols.cpp

clean:
	rm -f print_vols.o
