# --- DAGMC option.
CXXFLAGS = $(CFLAGS)

MOAB_DIR = $(HOME)/.local
CUBIT_DIR = 

include $(MOAB_DIR)/lib/moab.make   
CUBIT_LINK_PATH=$(CUBIT_DIR)

ifneq (,$(CUBIT_LINK_PATH))
  # Cubit-based MOAB build specified; ensure library paths work
  DAGMC_CFLAGS += -DCUBIT_LIBS_PRESENT
  MOAB_LDFLAGS += -Wl,-rpath=$(CUBIT_LINK_PATH)
endif

CPP_FLAGS += $(MOAB_CPPFLAGS)
CXXFLAGS += $(MOAB_CXXFLAGS) $(DAGMC_CFLAGS) 
INCLUDES += $(MOAB_INCLUDES)
LDFLAGS = $(MOAB_LDFLAGS) $(CXX_FORTRAN_LDFLAGS) 

DAGMC_LIBS += $(MOAB_LIBS_LINK) -ldagmc -lMOAB -L$(HOME)/UW/research/software/DAGMC/bld -ldagmciface -ldagtally -lstdc++

DAGMC_MOD=  dagmc_mod$(OBJF)
