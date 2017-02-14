message(STATUS "Building DAG-MCNP6")

enable_language(Fortran)

# Get date and time
# modifed from http://www.cmake.org/pipermail/cmake/2010-November/041034.html
macro (DATETIME RESULT)
  execute_process(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE ${RESULT})
  string(REGEX REPLACE "(..)/(..)/..(..).*" "\"\\2/\\1/\\3\"" ${RESULT} ${${RESULT}})
endmacro (DATETIME)
DATETIME(DATE)

# Source files
set(MCNP_SRC_FILES)
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/dedx/*.F90")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/fluka89/*.F90")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/hexs/*.F90")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/import/*.cpp")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/import/*.f")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/import/*.F")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/import/*.F90")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/import/cgm/*.c")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/import/cgm/*.cpp")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/lcs/*.F90")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/meshtal/*.F90")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/regl/*.F90")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/spabi/*.F90")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/src/*.c")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/src/*.F90")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/src/partisn/*.F90")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/src/utils/*.c")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/src/utils/*.F90")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "../*.cpp")
if (MPI_BUILD)
  file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/dotcomm/*.c")
endif (MPI_BUILD)

# Figure out which version of MCNP we're compiling based on which source files are present
set(MCNP_VERSION Unknown)
if ("/${MCNP_SRC_FILES};" MATCHES "/main.F90;")
  set(MCNP_VERSION 6_beta2)
endif ()
# These files are present in 6.1 but not in 6_beta2
if ("/${MCNP_SRC_FILES};" MATCHES "/cosmic_src.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/cosmic_uod.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/inelastic_xs_grid_mod.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/italft.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/kadjoint_data_mod.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/ksen_interface_mod.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/ksen_xs_mod.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/limited_interaction.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/model_cross_sections.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/otf_mod.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/read_background.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/talfft.F90;")
  set(MCNP_VERSION 6.1)
endif ()
# These files are present in 6.1.1beta but not in 6.1
if ("/${MCNP_SRC_FILES};" MATCHES "/cerenkov_mod.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/cgm_interface.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/parse_ichoic.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/reflect_light.F90;")
  set(MCNP_VERSION 6.1.1beta)
endif ()
if (${MCNP_VERSION} STREQUAL Unknown)
  message(FATAL_ERROR "MCNP6 source files not found in mcnp6/Source/src")
endif (${MCNP_VERSION} STREQUAL Unknown)
if (NOT "/${MCNP_SRC_FILES};" MATCHES "/dagmc_mod.F90;")
  message(FATAL_ERROR "MCNP6 source version ${MCNP_VERSION} files were found, but dagmc_mod.F90 was not. Perhaps the patch was not applied?")
endif (NOT "/${MCNP_SRC_FILES};" MATCHES "/dagmc_mod.F90;")
message(STATUS "DAG-MCNP6 version ${MCNP_VERSION} source files found")

# C compiler flags
if (${CMAKE_C_COMPILER_ID} STREQUAL Intel)
  #set(CMAKE_C_FLAGS "${CMAKE_Fortran_FLAGS} -mcmodel=medium")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pc64")
else ()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")
endif (${CMAKE_C_COMPILER_ID} STREQUAL Intel)

# Fortran compiler flags
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
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
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)

# Linker flags
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-intel")
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)

# Include directories
include_directories(${CMAKE_CURRENT_LIST_DIR}/Source/import)
if (${MCNP_VERSION} STREQUAL 6.1.1beta)
  include_directories(${CMAKE_CURRENT_LIST_DIR}/Source/import/cgm)
endif (${MCNP_VERSION} STREQUAL 6.1.1beta)

# Variables in mcnp_env.F90
add_definitions(-DKODE="mcnp")
if (MPI_BUILD)
  add_definitions(-DVERS="6.mpi")
else ()
  add_definitions(-DVERS="6")
endif (MPI_BUILD)
if (${MCNP_VERSION} STREQUAL 6_beta2)
  add_definitions(-DTHREAD="MCNP6")
  add_definitions(-DTVERS="BETA2")
elseif (${MCNP_VERSION} STREQUAL 6.1)
  add_definitions(-DTHREAD="MCNP6")
  add_definitions(-DTVERS="1.0")
elseif (${MCNP_VERSION} STREQUAL 6.1.1beta)
  add_definitions(-DTHREAD="MCNP")
  add_definitions(-DTVERS="6.1.1b")
endif (${MCNP_VERSION} STREQUAL 6_beta2)
add_definitions(-DLODDAT=${DATE})

# Preprocessor directives
if (${MCNP_VERSION} STREQUAL 6_beta2)
  add_definitions(-DALPHA_VERSION=ON)
  add_definitions(-DLANL=ON)
endif (${MCNP_VERSION} STREQUAL 6_beta2)
add_definitions(-DLINUX=ON)
add_definitions(-DUNIX=ON)
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  add_definitions(-DINTEL=ON)
  add_definitions(-DDIRACCESS_RECL_WORDS=ON)
else ()
  add_definitions(-DGFORTRAN=ON)
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
add_definitions(-DLLNL_FISS=ON)
add_definitions(-DF90_INT_4BYTE=ON)
add_definitions(-DF90_REAL_8BYTE=ON)
if (${MCNP_VERSION} STREQUAL 6.1.1beta)
  add_definitions(-DCGM=ON)
endif (${MCNP_VERSION} STREQUAL 6.1.1beta)
add_definitions(-DCEM=ON)
add_definitions(-DCEM_M6=ON)
add_definitions(-DCINDER=ON)
add_definitions(-DLAQGSM=ON)
if (${MCNP_VERSION} STREQUAL 6.1 OR ${MCNP_VERSION} STREQUAL 6.1.1beta)
  add_definitions(-DUOD=ON)
endif (${MCNP_VERSION} STREQUAL 6.1 OR ${MCNP_VERSION} STREQUAL 6.1.1beta)
if (MESHTAL_DEBUG)
  add_definitions(-DMESHTAL_DEBUG=ON)
endif (MESHTAL_DEBUG)

# Make sure DATAPATH is set
if (MCNP6_DATAPATH)
  add_definitions(-DDPATH="${MCNP6_DATAPATH}")
else (MCNP6_DATAPATH)
  if (DEFINED ENV{DATAPATH})
    add_definitions(-DDPATH="$ENV{DATAPATH}")
  else (DEFINED ENV{DATAPATH})
    message(FATAL_ERROR "MCNP DATAPATH not set. Either the environment variable DATAPATH or the CMake variable DMCNP6_DATAPATH must be set.")
  endif (DEFINED ENV{DATAPATH})
endif (MCNP6_DATAPATH)

# MCNP Plotter
if (MCNP6_PLOT)
  find_package(X11 REQUIRED)
  add_definitions(-DPLOT=ON)
  if (${MCNP_VERSION} STREQUAL 6_beta2)
    add_definitions(-DMCPLOT=ON)
  endif (${MCNP_VERSION} STREQUAL 6_beta2)
  add_definitions(-DXLIB=ON)
endif (MCNP6_PLOT)

# OpenMP
if (OPENMP_BUILD)
  find_package(OpenMP REQUIRED)
  add_definitions(-DOMP=ON)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")  # OpenMP_Fortran_FLAGS is undefined
endif (OPENMP_BUILD)

# MPI
if (MPI_BUILD)
  find_package(MPI REQUIRED)
  add_definitions(-DMPI=ON)
  add_definitions(-D_DOTCOMM_MPI=ON)
  include_directories(${MPI_INCLUDE_PATH})
  if (${MCNP_VERSION} STREQUAL 6_beta2)
    include_directories(${CMAKE_CURRENT_LIST_DIR}/Source/dotcomm)
  endif (${MCNP_VERSION} STREQUAL 6_beta2)
endif (MPI_BUILD)

# DAGMC
include_directories(${CMAKE_SOURCE_DIR}/common)
include_directories(${CMAKE_SOURCE_DIR}/tally)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Add the executable
if (MPI_BUILD)
  set(MCNP_EXEC "mcnp6.mpi")
else ()
  set(MCNP_EXEC "mcnp6")
endif (MPI_BUILD)
add_executable(${MCNP_EXEC} ${MCNP_SRC_FILES})

# Make sure the executable is compiled as Fortran
set_target_properties(${MCNP_EXEC} PROPERTIES LINKER_LANGUAGE Fortran)

# Link appropriate libraries
if (MCNP6_PLOT)
  target_link_libraries(${MCNP_EXEC} ${X11_LIBRARIES})
endif (MCNP6_PLOT)
if (MPI_BUILD)
  target_link_libraries(${MCNP_EXEC} ${MPI_Fortran_LIBRARIES})
endif (MPI_BUILD)
target_link_libraries(${MCNP_EXEC} ${HDF5_hdf5_hl_LIBRARY})
target_link_libraries(${MCNP_EXEC} ${HDF5_hdf5_LIBRARY})
target_link_libraries(${MCNP_EXEC} ${DAG_LINK_LIBS})
target_link_libraries(${MCNP_EXEC} dagmciface)
target_link_libraries(${MCNP_EXEC} dagtally)
target_link_libraries(${MCNP_EXEC} pyne_dagmc)
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "stdc++")
  set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "stdc++")
  set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "")
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)

# Install DAG-MCNP6
install(TARGETS ${MCNP_EXEC} DESTINATION ${INSTALL_BIN_DIR}
                             PERMISSIONS ${PERMS})
