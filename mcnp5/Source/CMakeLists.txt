cmake_minimum_required(VERSION 2.8)
project(MCNP5)

# modifed from http://www.cmake.org/pipermail/cmake/2010-November/041034.html
MACRO (DATETIME RESULT)
  EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE ${RESULT})
  STRING(REGEX REPLACE "(..)/(..)/..(..).*" "\"\\2/\\1/\\3\"" ${RESULT} ${${RESULT}})
ENDMACRO (DATETIME)
DATETIME(DATE)

# Build type
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif (NOT CMAKE_BUILD_TYPE)

# Make sure DATAPATH is set
if (MCNP5_DATAPATH)
  add_definitions(-DDPATH="${MCNP5_DATAPATH}")
else (MCNP5_DATAPATH)
  if(DEFINED ENV{DATAPATH})
    add_definitions(-DDPATH="$ENV{DATAPATH}")
  else(DEFINED ENV{DATAPATH})
    message(FATAL_ERROR "MCNP DATAPATH not set. Either the environment variable DATAPATH or the CMake variable DMCNP5_DATAPATH must be set.")
  endif(DEFINED ENV{DATAPATH})
endif (MCNP5_DATAPATH)

# Compilers
#if (MPI_BUILD)
#  execute_process(COMMAND which mpicc OUTPUT_VARIABLE CMAKE_C_COMPILER)
#  execute_process(COMMAND which mpic++ OUTPUT_VARIABLE CMAKE_CXX_COMPILER)
#  execute_process(COMMAND which mpif90 OUTPUT_VARIABLE CMAKE_Fortran_COMPILER)
#elseif (MPI_BUILD)
#  execute_process(COMMAND which gcc OUTPUT_VARIABLE CMAKE_C_COMPILER)
#  execute_process(COMMAND which g++ OUTPUT_VARIABLE CMAKE_CXX_COMPILER)
#  execute_process(COMMAND which gfortran OUTPUT_VARIABLE CMAKE_Fortran_COMPILER)
#endif (MPI_BUILD)

# Fortran compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mieee-fp")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdollar-ok")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-sign-zero")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-double-8")

# Preprocessor directives
add_definitions(-DKODE="mcnp")
add_definitions(-DTHREAD="MCNP5")
add_definitions(-DTVERS="1.60")
add_definitions(-DLODDAT=${DATE})
add_definitions(-DLINUX=ON)
add_definitions(-DUNIX=ON)
add_definitions(-DGFORTRAN=ON)
add_definitions(-DF90_INT_4BYTE=ON)
add_definitions(-DF90_REAL_8BYTE=ON)

# MCNP Plotter
if (MCNP_PLOT)
  find_package(X11)
  if (X11_FOUND)
    add_definitions(-DPLOT=ON)
    add_definitions(-DMCPLOT=ON)
    add_definitions(-DXLIB=ON)
  elseif (X11_FOUND)
    message(FATAL_ERROR "X11 not found.")
  endif (X11_FOUND)
endif (MCNP_PLOT)

# OpenMP
if (OMP_BUILD)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    add_definitions(-DOMP=ON)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")  # OpenMP_Fortran_FLAGS is undefined
  elseif (OPENMP_FOUND)
    message(FATAL_ERROR "OpenMP not found.")
  endif (OPENMP_FOUND)
endif (OMP_BUILD)

# MPI
if (MPI_BUILD)
  find_package(MPI)
  if (MPI_FOUND)
    add_definitions(-DVERS="5.mpi")
    add_definitions(-DMPI=ON)
    add_definitions(-DDMMP_NAME="mcnp5.mpi")
    add_definitions(-D_DOTCOMM_MPI=ON)
    add_definitions(-DDOTCOMM_INTERNAL="mpi")  # probably don't need this?
    include_directories(${MPI_INCLUDE_PATH})
    include_directories(${SRC_DIR}/dotcomm/include)
    include_directories(${SRC_DIR}/dotcomm/src/internals/mpi)
  elseif (MPI_FOUND)
    message(FATAL_ERROR "MPI not found.")
  endif (MPI_FOUND)
else (MPI_BUILD)
  add_definitions(-DVERS="5")
endif (MPI_BUILD)

# Source files
include(${SRC_DIR}/src/CMakeLists.txt)
if (MPI_BUILD)
  include(${SRC_DIR}/dotcomm/CMakeLists.txt)
endif (MPI_BUILD)

# Add the executable
if (MPI_BUILD)
  set(MCNP_EXEC "mcnp5.mpi")
else (MPI_BUILD)
  set(MCNP_EXEC "mcnp5")
endif (MPI_BUILD)
add_executable(${MCNP_EXEC} ${SRC_FILES})

# Make sure the executable is compiled as Fortran
set_target_properties(${MCNP_EXEC} PROPERTIES LINKER_LANGUAGE Fortran)

# Link appropriate libraries
if (MCNP_PLOT)
  target_link_libraries(${MCNP_EXEC} ${X11_LIBRARIES})
endif (MCNP_PLOT)
if (MPI_BUILD)
  target_link_libraries(${MCNP_EXEC} ${MPI_LIBRARIES})
endif (MPI_BUILD)

# Install the targets
set(CMAKE_SKIP_RPATH TRUE)
install(TARGETS ${MCNP_EXEC} RUNTIME DESTINATION ${INSTALL_BIN_DIR} PERMISSIONS ${PERMS})
