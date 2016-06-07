message(STATUS "Building DAG-MCNP5")

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
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "../Source/src/*.F90")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "../Source/src/*.c")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "../../*.cpp")
if (MPI_BUILD)
  file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "../Source/dotcomm/src/*.F90")
  file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "../Source/dotcomm/src/internals/mpi/*.c")
endif (MPI_BUILD)

# Figure out which version of MCNP we're compiling based on which source files are present
set(MCNP_VERSION Unknown)
if ("/${MCNP_SRC_FILES};" MATCHES "/ckchar.F90;")
  set(MCNP_VERSION 1.40)
endif ()
if ("/${MCNP_SRC_FILES};" MATCHES "/annihilation_gammas.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/phtvr_mod.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/read_line.F90;")
  set(MCNP_VERSION 1.51)
endif ()
if ("/${MCNP_SRC_FILES};" MATCHES "/isotopic_mod.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/kadjoint_mod.F90;" AND
    "/${MCNP_SRC_FILES};" MATCHES "/kinetics_mod.F90;")
  set(MCNP_VERSION 1.60)
endif ()
if (${MCNP_VERSION} STREQUAL Unknown)
  message(FATAL_ERROR "MCNP5 source files not found in mcnp5/Source/src")
endif (${MCNP_VERSION} STREQUAL Unknown)
if (NOT "/${MCNP_SRC_FILES};" MATCHES "/dagmc_mod.F90;")
  message(FATAL_ERROR "MCNP5 source version ${MCNP_VERSION} files were found, but dagmc_mod.F90 was not. Perhaps the patch was not applied?")
endif (NOT "/${MCNP_SRC_FILES};" MATCHES "/dagmc_mod.F90;")
message(STATUS "DAG-MCNP5 version ${MCNP_VERSION} source files found")

# C compiler flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")

# Fortran compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mieee-fp")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdollar-ok")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-sign-zero")
if (MCNP5_PLOT)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
endif (MCNP5_PLOT)
if (LIB_TYPE STREQUAL STATIC)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wl,--allow-multiple-definition")
endif (LIB_TYPE STREQUAL STATIC)

# Preprocessor directives
add_definitions(-DKODE="mcnp")
if (MPI_BUILD)
  add_definitions(-DVERS="5.mpi")
else ()
  add_definitions(-DVERS="5")
endif (MPI_BUILD)
if (${MCNP_VERSION} STREQUAL 1.60)
  add_definitions(-DTHREAD="MCNP5")
else ()
  add_definitions(-DTHREAD="MCNP5_RSICC")
endif (${MCNP_VERSION} STREQUAL 1.60)
add_definitions(-DTVERS="${MCNP_VERSION}")
add_definitions(-DLODDAT=${DATE})
add_definitions(-DLINUX=ON)
add_definitions(-DUNIX=ON)
add_definitions(-DGFORTRAN=ON)
add_definitions(-DF90_INT_4BYTE=ON)
add_definitions(-DF90_REAL_8BYTE=ON)

# Make sure DATAPATH is set
if (MCNP5_DATAPATH)
  add_definitions(-DDPATH="${MCNP5_DATAPATH}")
else ()
  if(DEFINED ENV{DATAPATH})
    add_definitions(-DDPATH="$ENV{DATAPATH}")
  else ()
    message(FATAL_ERROR "MCNP DATAPATH not set. Either the environment variable DATAPATH or the CMake variable DMCNP5_DATAPATH must be set.")
  endif (DEFINED ENV{DATAPATH})
endif (MCNP5_DATAPATH)

# MCNP Plotter
if (MCNP5_PLOT)
  find_package(X11 REQUIRED)
  add_definitions(-DPLOT=ON)
  add_definitions(-DMCPLOT=ON)
  add_definitions(-DXLIB=ON)
endif (MCNP5_PLOT)

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
  add_definitions(-DMULTP=ON)
  add_definitions(-DMPI=ON)
  add_definitions(-DDMMP_NAME="mcnp5.mpi")
  add_definitions(-D_DOTCOMM_MPI=ON)
  include_directories(${MPI_INCLUDE_PATH})
  include_directories(${CMAKE_CURRENT_LIST_DIR}/../Source/dotcomm/include)
  include_directories(${CMAKE_CURRENT_LIST_DIR}/../Source/dotcomm/src/internals/mpi)
endif (MPI_BUILD)

# DAGMC
include_directories(${PROJECT_SOURCE_DIR}/tally)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Add the executable
if (MPI_BUILD)
  set(MCNP_EXEC "mcnp5.mpi")
else ()
  set(MCNP_EXEC "mcnp5")
endif (MPI_BUILD)
add_executable(${MCNP_EXEC} ${MCNP_SRC_FILES})

# Make sure the executable is compiled as Fortran
set_target_properties(${MCNP_EXEC} PROPERTIES LINKER_LANGUAGE Fortran)

# Link appropriate libraries
if (MCNP5_PLOT)
  target_link_libraries(${MCNP_EXEC} ${X11_LIBRARIES})
endif (MCNP5_PLOT)
if (MPI_BUILD)
  target_link_libraries(${MCNP_EXEC} ${MPI_LIBRARIES})
endif (MPI_BUILD)
target_link_libraries(${MCNP_EXEC} ${HDF5_hdf5_hl_LIBRARY})
target_link_libraries(${MCNP_EXEC} ${HDF5_hdf5_LIBRARY})
target_link_libraries(${MCNP_EXEC} ${DAG_LINK_LIBS})
target_link_libraries(${MCNP_EXEC} dagmciface)
target_link_libraries(${MCNP_EXEC} dagtally)
target_link_libraries(${MCNP_EXEC} pyne_dagmc)

# Install DAG-MCNP5
install(TARGETS ${MCNP_EXEC} DESTINATION ${INSTALL_BIN_DIR}
                             PERMISSIONS ${PERMS})
