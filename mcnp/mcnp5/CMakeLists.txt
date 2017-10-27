message(STATUS "Building DAG-MCNP5")

# Get date and time
execute_process(COMMAND date +%m/%d/%y OUTPUT_VARIABLE ENV_DATE OUTPUT_STRIP_TRAILING_WHITESPACE)

set(LINK_LIBS)

# Source files
set(MCNP_SRC_FILES)
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/src/*.F90")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/src/*.c")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/src/*.cpp")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "${CMAKE_SOURCE_DIR}/mcnp/*.cpp")
if (MPI_BUILD)
  file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/dotcomm/src/*.F90")
  file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/dotcomm/src/internals/mpi/*.c")
endif ()

# Figure out which version of MCNP we're compiling based on which source files are present
set(MCNP_VERSION Unknown)
if ("/${MCNP_SRC_FILES};" MATCHES "/main.F90;")
  set(MCNP_VERSION 1.40)
endif ()
# These files are present in 1.51 but not in 1.40
if ("/${MCNP_SRC_FILES};" MATCHES "/annihilation_gammas.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/phtvr_mod.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/read_line.F90;")
  set(MCNP_VERSION 1.51)
endif ()
# These files are present in 1.60 but not in 1.51
if ("/${MCNP_SRC_FILES};" MATCHES "/isotopic_mod.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/kadjoint_mod.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/kinetics_mod.F90;")
  set(MCNP_VERSION 1.60)
endif ()
if (${MCNP_VERSION} STREQUAL Unknown)
  message(FATAL_ERROR "MCNP5 source files not found in mcnp5/Source/src")
endif ()
if (NOT "/${MCNP_SRC_FILES};" MATCHES "/dagmc_mod.F90;")
  message(FATAL_ERROR "MCNP5 source version ${MCNP_VERSION} files were found, but dagmc_mod.F90 was not. Perhaps the patch was not applied?")
endif ()
message(STATUS "DAG-MCNP5 version ${MCNP_VERSION} source files found")

# Optimization flags
set(CMAKE_C_FLAGS_RELEASE              "-O1"   )
set(CMAKE_C_FLAGS_RELWITHDEBINFO       "-O1 -g")
set(CMAKE_CXX_FLAGS_RELEASE            "-O1"   )
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO     "-O1 -g")
set(CMAKE_Fortran_FLAGS_RELEASE        "-O1"   )
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O1 -g")

# C compiler flags
if (${CMAKE_C_COMPILER_ID} STREQUAL Intel)
  #set(CMAKE_C_FLAGS "${CMAKE_Fortran_FLAGS} -mcmodel=medium")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pc64")
else ()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")
endif ()

# Fortran compiler flags
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
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
  if (MCNP5_PLOT)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
  endif ()
endif ()

# Preprocessor directives
add_definitions(-DKODE="mcnp")
if (MPI_BUILD)
  add_definitions(-DVERS="5.mpi")
else ()
  add_definitions(-DVERS="5")
endif ()
if (${MCNP_VERSION} STREQUAL 1.60)
  add_definitions(-DTHREAD="MCNP5")
else ()
  add_definitions(-DTHREAD="MCNP5_RSICC")
endif ()
add_definitions(-DTVERS="${MCNP_VERSION}")
add_definitions(-DLODDAT="${ENV_DATE}")
add_definitions(-DLINUX=ON)
add_definitions(-DUNIX=ON)
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
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
if (MCNP5_PLOT)
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
  add_definitions(-DMULTP=ON)
  add_definitions(-DMPI=ON)
  add_definitions(-DDMMP_NAME="mcnp5.mpi")
  add_definitions(-D_DOTCOMM_MPI=ON)
  list(APPEND LINK_LIBS ${MPI_Fortran_LIBRARIES})
  include_directories(${MPI_INCLUDE_PATH})
  include_directories(${CMAKE_CURRENT_LIST_DIR}/Source/dotcomm/include)
  include_directories(${CMAKE_CURRENT_LIST_DIR}/Source/dotcomm/src/internals/mpi)
endif ()

# DAGMC
include_directories(${CMAKE_SOURCE_DIR}/dagmc)
include_directories(${CMAKE_SOURCE_DIR}/tally)
if (BUILD_STATIC_EXE)
  list(APPEND LINK_LIBS dagmc-static dagtally-static)
else ()
  list(APPEND LINK_LIBS dagmc dagtally)
endif ()

# Executable
if (MPI_BUILD)
  set(MCNP_EXEC "mcnp5.mpi")
else ()
  set(MCNP_EXEC "mcnp5")
endif ()
add_executable(${MCNP_EXEC} ${MCNP_SRC_FILES})
set_target_properties(${MCNP_EXEC} PROPERTIES LINKER_LANGUAGE Fortran)

# Install DAG-MCNP5
target_link_libraries(${MCNP_EXEC} ${LINK_LIBS})
install(TARGETS ${MCNP_EXEC} DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
