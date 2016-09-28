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
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/hexs/*.F90")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/import/*.cpp")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/import/*.f")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/import/*.F")
file(GLOB MCNP_SRC_FILES ${MCNP_SRC_FILES} "Source/import/*.F90")
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

# C compiler flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")

# Fortran compiler flags
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pc64")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -diag-disable 8290,8291")
else ()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mieee-fp")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdollar-ok")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-sign-zero")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-double-8")
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)

# Include directories
include_directories(${CMAKE_CURRENT_LIST_DIR}/Source/import)
include_directories(${PROJECT_SOURCE_DIR}/tally)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Preprocessor directives
add_definitions(-DKODE="mcnp")
if (MPI_BUILD)
  add_definitions(-DVERS="6.mpi")
else ()
  add_definitions(-DVERS="6")
endif (MPI_BUILD)
add_definitions(-DTHREAD="MCNP6")
add_definitions(-DTVERS="1.0")
add_definitions(-DLODDAT=${DATE})
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
add_definitions(-DCEM=ON)
add_definitions(-DCEM_M6=ON)
add_definitions(-DCINDER=ON)
add_definitions(-DLAQGSM=ON)
add_definitions(-DUOD=ON)

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
  include_directories(${MPI_INCLUDE_PATH})
endif (MPI_BUILD)

# Add the executable
if (MPI_BUILD)
  set(MCNP_EXEC "mcnp6.mpi")
else (MPI_BUILD)
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

# Install MCNP6
install(TARGETS ${MCNP_EXEC} DESTINATION ${INSTALL_BIN_DIR}
                             PERMISSIONS ${PERMS})
