# $Id: GNUmakefile,v 1.1 1999-01-07 16:05:40 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exampleN01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

DAGSOLID_HOME  := /data/opt/andy_dag/DAGMC/Geant4/dagsolid
DAGSOLID_LIB   := ${DAGSOLID_HOME}/lib
DAGSOLID_INC   := -I${DAGSOLID_HOME}/include

HDF5_DIR       := /data/opt/dagmc/hdf5

MOAB_DIR       := /data/opt/moab_dev/
MOAB_LDFLAGS   := -L/data/opt/dagmc/hdf5/lib 
MOAB_LIBDIR    := ${MOAB_DIR}/lib
MOAB_INCLUDES  :=-I${MOAB_DIR}/include
MOAB_LIBS_LINK := ${MOAB_LDFLAGS} -L${MOAB_LIBDIR} -lMOAB -lm  -lhdf5 -g -pg

PYNE_DIR       := /home/davisa/.local/
PYNE_LDFLAGS   := -L${PYNE_DIR}/lib ${MOAB_LDFLAGS}
PYNE_LIBDIR    := ${PYNE_DIR}/lib
PYNE_INCLUDES  :=-I${PYNE_DIR}/include/pyne -I${HDF5_DIR}/include
PYNE_LIBS_LINK :=-L${PYNE_LIBDIR} -lpyne -lhdf5 


CPPFLAGS  += $(MOAB_INCLUDES) $(DAGSOLID_INC) $(PYNE_INCLUDES) -g -pg
CXXFLAGS  = $(CPPFLAGS) -fPIC
LDFLAGS   += $(MOAB_LIBS_LINK) -g -pg 
LDLIBS += ${MOAB_DIR}/lib/libdagmc.so
LDLIBS += ${MOAB_LIBDIR}/libMOAB.so
LDLIBS += ${DAGSOLID_LIB}/libdagsolid.so
LDLIBS += ${PYNE_LIBDIR}/libpyne.so
