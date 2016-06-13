include ${MOAB_MAKE}

INC = -I../
MOAB_CXXFLAGS =  -Wall -pipe -pedantic -Wno-long-long ${INC}
MOAB_CFLAGS = -Wall -pipe -pedantic -Wno-long-long ${INC}
CXXFLAGS += ${MOAB_CXXFLAGS} -g 
CC = g++
LD_FLAGS = -g
CPPFLAGS += ${MOAB_INCLUDES} -g
CFLAGS   += ${MOAB_CFLAGS} -g
# add -g -pg to both CXX and LD flags to profile

all: test_cyl 

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

test_funcs.o: gen.o
	$(CC) $(CXXFLAGS) $(MOAB_INCLUDES) -c test_funcs.cpp

test_cyl: test_cyl.o test_funcs.o gen.o arc.o zip.o cleanup.o cw_func.o mw_func.o
	$(CC) $(LD_FLAGS) -o test_cyl test_funcs.o test_cyl.o gen.o arc.o zip.o cleanup.o cw_func.o  \
	mw_func.o ${MOAB_LIBS_LINK} -ldagmc
check: test_cyl
	./cone_tests && echo && ./test_cyl

clean:
	rm -f test_cyl *.o *_mod*
