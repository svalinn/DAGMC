include /data/opt/dagmc/moab/lib/moab.make
CXXFLAGS += ${MOAB_CXXFLAGS} -g -pg
CC = mpicxx
LD_FLAGS = -pg -std=c++11 
CPPFLAGS += ${MOAB_INCLUDES} -pg
CFLAGS   += ${MOAB_CFLAGS} -pg
# add -g -pg to both CXX and LD flags to profile

all: tartarus

tartarus: tartarus.o 
	$(CC) $(LD_FLAGS) -o tartarus tartarus.o ${MOAB_LIBS_LINK} -ldagmc

tartarus.o: tartarus.cpp 
	$(CC) $(CXXFLAGS) ${MOAB_INCLUDES} -c tartarus.cpp $(LD_FLAGS) -ldagmc

clean:
	rm -f print_vols.o
