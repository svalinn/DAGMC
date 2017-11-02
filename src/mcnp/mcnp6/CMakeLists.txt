# Get date and time
execute_process(COMMAND date +%m/%d/%y OUTPUT_VARIABLE ENV_DATE OUTPUT_STRIP_TRAILING_WHITESPACE)

set(LINK_LIBS)

# Source files
file(GLOB SRC_FILES
  "Source/dedx/*.F90"
  "Source/fluka89/*.F90"
  "Source/hexs/*.F90"
  "Source/import/*.cpp"
  "Source/import/*.f"
  "Source/import/*.F"
  "Source/import/*.F90"
  "Source/import/cgm/*.c"
  "Source/import/cgm/*.cpp"
  "Source/lcs/*.F90"
  "Source/meshtal/*.F90"
  "Source/regl/*.F90"
  "Source/spabi/*.F90"
  "Source/src/*.c"
  "Source/src/*.cpp"
  "Source/src/*.F90"
  "Source/src/partisn/*.F90"
  "Source/src/utils/*.c"
  "Source/src/utils/*.F90"
)

# Figure out which version of MCNP we're compiling based on which source files are present
set(MCNP_VERSION Unknown)
if ("/${SRC_FILES};" MATCHES "/main.F90;")
  set(MCNP_VERSION 6_beta2)
endif ()
# These files are present in 6.1 but not in 6_beta2
if ("/${SRC_FILES};" MATCHES "/cosmic_src.F90;" AND
    "/${SRC_FILES};" MATCHES "/cosmic_uod.F90;" AND
    "/${SRC_FILES};" MATCHES "/inelastic_xs_grid_mod.F90;" AND
    "/${SRC_FILES};" MATCHES "/italft.F90;" AND
    "/${SRC_FILES};" MATCHES "/kadjoint_data_mod.F90;" AND
    "/${SRC_FILES};" MATCHES "/ksen_interface_mod.F90;" AND
    "/${SRC_FILES};" MATCHES "/ksen_xs_mod.F90;" AND
    "/${SRC_FILES};" MATCHES "/limited_interaction.F90;" AND
    "/${SRC_FILES};" MATCHES "/model_cross_sections.F90;" AND
    "/${SRC_FILES};" MATCHES "/otf_mod.F90;" AND
    "/${SRC_FILES};" MATCHES "/read_background.F90;" AND
    "/${SRC_FILES};" MATCHES "/talfft.F90;")
  set(MCNP_VERSION 6.1)
endif ()
# These files are present in 6.1.1beta but not in 6.1
if ("/${SRC_FILES};" MATCHES "/cerenkov_mod.F90;" AND
    "/${SRC_FILES};" MATCHES "/cgm_interface.F90;" AND
    "/${SRC_FILES};" MATCHES "/parse_ichoic.F90;" AND
    "/${SRC_FILES};" MATCHES "/reflect_light.F90;")
  set(MCNP_VERSION 6.1.1beta)
endif ()
if (MCNP_VERSION STREQUAL "Unknown")
  message(FATAL_ERROR "MCNP6 source files not found in mcnp6/Source/src")
endif ()
if (NOT "/${SRC_FILES};" MATCHES "/dagmc_mod.F90;")
  message(FATAL_ERROR "MCNP6 source version ${MCNP_VERSION} files were found, but dagmc_mod.F90 was not. Perhaps the patch was not applied?")
endif ()
message(STATUS "DAG-MCNP6 version ${MCNP_VERSION} source files found")

# Optimization flags
set(CMAKE_C_FLAGS_RELEASE              "-O1"   )
set(CMAKE_C_FLAGS_RELWITHDEBINFO       "-O1 -g")
set(CMAKE_CXX_FLAGS_RELEASE            "-O1"   )
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO     "-O1 -g")
set(CMAKE_Fortran_FLAGS_RELEASE        "-O1"   )
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O1 -g")

# C compiler flags
if (CMAKE_C_COMPILER_ID STREQUAL "Intel")
  #set(CMAKE_C_FLAGS "${CMAKE_Fortran_FLAGS} -mcmodel=medium")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pc64")
else ()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")
endif ()

# Fortran compiler flags
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-vec")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -heap-arrays 1024")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -traceback")
  #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mcmodel=medium")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pc64")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -diag-disable 8290,8291,8577")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
else ()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mieee-fp")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdollar-ok")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-sign-zero")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-double-8")
endif ()

# Include directories
include_directories(${CMAKE_CURRENT_LIST_DIR}/Source/import)
if (MCNP_VERSION STREQUAL "6.1.1beta")
  include_directories(${CMAKE_CURRENT_LIST_DIR}/Source/import/cgm)
endif ()

# Variables in mcnp_env.F90
add_definitions(-DKODE="mcnp")
if (MCNP_VERSION STREQUAL "6_beta2")
  add_definitions(-DTHREAD="MCNP6")
  add_definitions(-DTVERS="BETA2")
elseif (MCNP_VERSION STREQUAL "6.1")
  add_definitions(-DTHREAD="MCNP6")
  add_definitions(-DTVERS="1.0")
elseif (MCNP_VERSION STREQUAL "6.1.1beta")
  add_definitions(-DTHREAD="MCNP")
  add_definitions(-DTVERS="6.1.1b")
endif ()
add_definitions(-DLODDAT="${ENV_DATE}")

# Preprocessor directives
if (MCNP_VERSION STREQUAL "6_beta2")
  add_definitions(-DALPHA_VERSION=ON)
  add_definitions(-DLANL=ON)
endif ()
add_definitions(-DLINUX=ON)
add_definitions(-DUNIX=ON)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  add_definitions(-DINTEL=ON)
  add_definitions(-DDIRACCESS_RECL_WORDS=ON)
else ()
  add_definitions(-DGFORTRAN=ON)
endif ()
add_definitions(-DLLNL_FISS=ON)
add_definitions(-DF90_INT_4BYTE=ON)
add_definitions(-DF90_REAL_8BYTE=ON)
if (MCNP_VERSION STREQUAL "6.1.1beta")
  add_definitions(-DCGM=ON)
endif ()
add_definitions(-DCEM=ON)
add_definitions(-DCEM_M6=ON)
add_definitions(-DCINDER=ON)
add_definitions(-DLAQGSM=ON)
if (MCNP_VERSION STREQUAL "6.1" OR MCNP_VERSION STREQUAL "6.1.1beta")
  add_definitions(-DUOD=ON)
endif ()
if (MESHTAL_DEBUG)
  add_definitions(-DMESHTAL_DEBUG=ON)
endif ()

# Make sure DATAPATH is set
if (MCNP6_DATAPATH)
  add_definitions(-DDPATH="${MCNP6_DATAPATH}")
else ()
  if (DEFINED ENV{DATAPATH})
    add_definitions(-DDPATH="$ENV{DATAPATH}")
  else ()
    message(FATAL_ERROR "MCNP DATAPATH not set. Either the environment variable DATAPATH or the CMake variable DMCNP6_DATAPATH must be set.")
  endif ()
endif ()

# MCNP Plotter
if (MCNP6_PLOT)
  find_package(X11 REQUIRED)
  add_definitions(-DPLOT=ON)
  if (MCNP_VERSION STREQUAL "6_beta2")
    add_definitions(-DMCPLOT=ON)
  endif ()
  add_definitions(-DXLIB=ON)
  list(APPEND LINK_LIBS ${X11_LIBRARIES})
  if (BUILD_STATIC_EXE)
    list(APPEND LINK_LIBS xcb pthread Xau Xdmcp)
  endif ()
endif ()

# OpenMP
if (OPENMP_BUILD)
  find_package(OpenMP REQUIRED)
  add_definitions(-DOMP=ON)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif ()

# MPI
if (MPI_BUILD)
  find_package(MPI REQUIRED)
  file(GLOB SRC_FILES ${SRC_FILES} "Source/dotcomm/*.c")
  add_definitions(-DVERS="6.mpi")
  add_definitions(-DMPI=ON)
  add_definitions(-D_DOTCOMM_MPI=ON)
  list(APPEND LINK_LIBS ${MPI_Fortran_LIBRARIES})
  include_directories(${MPI_INCLUDE_PATH})
  if (MCNP_VERSION STREQUAL "6_beta2")
    include_directories(${CMAKE_CURRENT_LIST_DIR}/Source/dotcomm)
  endif ()
else ()
  add_definitions(-DVERS="6")
endif ()

# DAGMC
list(APPEND SRC_FILES
  $<TARGET_OBJECTS:mcnp_funcs>
  $<TARGET_OBJECTS:meshtal_funcs>
)
list(APPEND LINK_LIBS dagmc pyne_dagmc uwuw dagtally)

# Executable
if (MPI_BUILD)
  set(MCNP_EXEC "mcnp6.mpi")
else ()
  set(MCNP_EXEC "mcnp6")
endif ()

get_directory_property(COMPDEFS COMPILE_DEFINITIONS)
message(STATUS "DAG-MCNP6 C_FLAGS: ${CMAKE_C_FLAGS}")
message(STATUS "DAG-MCNP6 CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "DAG-MCNP6 Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
message(STATUS "DAG-MCNP6 C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}: ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
message(STATUS "DAG-MCNP6 CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
message(STATUS "DAG-MCNP6 Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}: ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
message(STATUS "DAG-MCNP6 compile definitions: ${COMPDEFS}")

# Install
dagmc_install_exe(${MCNP_EXEC})
set_target_properties(${MCNP_EXEC} PROPERTIES LINKER_LANGUAGE Fortran)
