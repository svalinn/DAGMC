message("")

set(LINK_LIBS)

# Source files
file(GLOB SRC_FILES
  "Source/src/*.c"
  "Source/src/*.cpp"
  "Source/src/*.F90"
)

# PyNE mesh source files
if (BUILD_MCNP_PYNE_SOURCE)
  message(STATUS "Building DAG-MCNP5 with PyNE mesh source routine support")
  # Remove MCNP's source.F90 file
  list(FILTER SRC_FILES EXCLUDE REGEX "/source.F90")
  list(APPEND SRC_FILES
    "pyne_source/source.F90"
    $<TARGET_OBJECTS:measure>
    $<TARGET_OBJECTS:source_sampling>
  )
endif ()

# Figure out which version of MCNP we're compiling based on which source files are present
set(MCNP_VERSION Unknown)
if ("/${SRC_FILES};" MATCHES "/main.F90;")
  set(MCNP_VERSION 514)
endif ()
# These files are present in 515 but not in 514
if ("/${SRC_FILES};" MATCHES "/annihilation_gammas.F90;" AND
    "/${SRC_FILES};" MATCHES "/phtvr_mod.F90;" AND
    "/${SRC_FILES};" MATCHES "/read_line.F90;")
  set(MCNP_VERSION 515)
endif ()
# These files are present in 516 but not in 515
if ("/${SRC_FILES};" MATCHES "/isotopic_mod.F90;" AND
    "/${SRC_FILES};" MATCHES "/kadjoint_mod.F90;" AND
    "/${SRC_FILES};" MATCHES "/kinetics_mod.F90;")
  set(MCNP_VERSION 516)
endif ()
if (MCNP_VERSION STREQUAL "Unknown")
  message(FATAL_ERROR "MCNP5 source files not found in mcnp5/Source/src")
elseif (NOT "/${SRC_FILES};" MATCHES "/dagmc_mod.F90;")
  message(FATAL_ERROR "MCNP${MCNP_VERSION} source files were found, but dagmc_mod.F90 was not. Perhaps the patch was not applied?")
else ()
  message(STATUS "DAG-MCNP${MCNP_VERSION} source files found")
endif ()

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
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
else ()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mieee-fp")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdollar-ok")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-sign-zero")
  if (BUILD_MCNP_PLOT)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
  endif ()
endif ()

# Preprocessor directives
add_definitions(-DKODE="mcnp")
if (MCNP_VERSION STREQUAL "514")
  add_definitions(-DTHREAD="MCNP5_RSICC")
  add_definitions(-DTVERS="1.40")
elseif (MCNP_VERSION STREQUAL "515")
  add_definitions(-DTHREAD="MCNP5_RSICC")
  add_definitions(-DTVERS="1.51")
elseif (MCNP_VERSION STREQUAL "516")
  add_definitions(-DTHREAD="MCNP5")
  add_definitions(-DTVERS="1.60")
endif ()
add_definitions(-DLODDAT="${ENV_DATE}")
add_definitions(-DLINUX=ON)
add_definitions(-DUNIX=ON)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  add_definitions(-DINTEL=ON)
  add_definitions(-DDIRACCESS_RECL_WORDS=ON)
else ()
  add_definitions(-DGFORTRAN=ON)
endif ()
add_definitions(-DF90_INT_4BYTE=ON)
add_definitions(-DF90_REAL_8BYTE=ON)
if (MESHTAL_DEBUG)
  add_definitions(-DMESHTAL_DEBUG=ON)
endif ()

# Make sure DATAPATH is set
if (MCNP5_DATAPATH)
  add_definitions(-DDPATH="${MCNP5_DATAPATH}")
else ()
  if (DEFINED ENV{DATAPATH})
    add_definitions(-DDPATH="$ENV{DATAPATH}")
  else ()
    message(FATAL_ERROR "MCNP DATAPATH not set. Either the environment variable DATAPATH or the CMake variable DMCNP5_DATAPATH must be set.")
  endif ()
endif ()

# MCNP Plotter
if (BUILD_MCNP_PLOT)
  find_package(X11 REQUIRED)
  add_definitions(-DPLOT=ON)
  add_definitions(-DMCPLOT=ON)
  add_definitions(-DXLIB=ON)
  list(APPEND LINK_LIBS ${X11_LIBRARIES})
  if (BUILD_STATIC_EXE)
    list(APPEND LINK_LIBS xcb pthread Xau Xdmcp)
  endif ()
endif ()

# OpenMP
if (BUILD_MCNP_OPENMP)
  find_package(OpenMP REQUIRED)
  add_definitions(-DOMP=ON)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif ()

# MPI
if (BUILD_MCNP_MPI)
  find_package(MPI REQUIRED)
  file(GLOB SRC_FILES ${SRC_FILES} "Source/dotcomm/src/*.F90")
  file(GLOB SRC_FILES ${SRC_FILES} "Source/dotcomm/src/internals/mpi/*.c")
  add_definitions(-DVERS="5.mpi")
  add_definitions(-DMULTP=ON)
  add_definitions(-DMPI=ON)
  add_definitions(-DDMMP_NAME="mcnp5.mpi")
  add_definitions(-D_DOTCOMM_MPI=ON)
  list(APPEND LINK_LIBS ${MPI_Fortran_LIBRARIES})
  include_directories(${MPI_INCLUDE_PATH})
  include_directories(${CMAKE_CURRENT_LIST_DIR}/Source/dotcomm/include)
  include_directories(${CMAKE_CURRENT_LIST_DIR}/Source/dotcomm/src/internals/mpi)
else ()
  add_definitions(-DVERS="5")
endif ()

# DAGMC
list(APPEND SRC_FILES
  $<TARGET_OBJECTS:mcnp_funcs>
  $<TARGET_OBJECTS:meshtal_funcs>
)
list(APPEND LINK_LIBS dagmc pyne_dagmc uwuw dagtally)

# Executable
if (BUILD_MCNP_MPI)
  set(MCNP_EXEC "mcnp5.mpi")
else ()
  set(MCNP_EXEC "mcnp5")
endif ()

get_directory_property(COMPDEFS COMPILE_DEFINITIONS)
message(STATUS "DAG-MCNP5 C_FLAGS: ${CMAKE_C_FLAGS}")
message(STATUS "DAG-MCNP5 CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "DAG-MCNP5 Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
message(STATUS "DAG-MCNP5 C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}: ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
message(STATUS "DAG-MCNP5 CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
message(STATUS "DAG-MCNP5 Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}: ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
message(STATUS "DAG-MCNP5 compile definitions: ${COMPDEFS}")

# Install
dagmc_install_exe(${MCNP_EXEC})
set_target_properties(${MCNP_EXEC} PROPERTIES LINKER_LANGUAGE Fortran)
