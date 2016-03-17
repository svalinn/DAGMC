project(DAGMC)

enable_language(Fortran)

# modifed from http://www.cmake.org/pipermail/cmake/2010-November/041034.html
macro (DATETIME RESULT)
  execute_process(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE ${RESULT})
  string(REGEX REPLACE "(..)/(..)/..(..).*" "\"\\2/\\1/\\3\"" ${RESULT} ${${RESULT}})
endmacro (DATETIME)
DATETIME(DATE)

# Build type
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif (NOT CMAKE_BUILD_TYPE)

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

# Make sure DATAPATH is set
if (MCNP5_DATAPATH)
  add_definitions(-DDPATH="${MCNP5_DATAPATH}")
else ()
  if(DEFINED ENV{DATAPATH})
    add_definitions(-DDPATH="$ENV{DATAPATH}")
  else ()
    message(FATAL_ERROR "MCNP DATAPATH not set. Either the environment variable DATAPATH or the CMake variable DMCNP5_DATAPATH must be set.")
  endif(DEFINED ENV{DATAPATH})
endif (MCNP5_DATAPATH)

# MCNP Plotter
if (MCNP_PLOT)
  find_package(X11 REQUIRED)
  add_definitions(-DPLOT=ON)
  add_definitions(-DMCPLOT=ON)
  add_definitions(-DXLIB=ON)
endif (MCNP_PLOT)

# OpenMP
if (OMP_BUILD)
  find_package(OpenMP REQUIRED)
  add_definitions(-DOMP=ON)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")  # OpenMP_Fortran_FLAGS is undefined
endif (OMP_BUILD)

# MPI
if (MPI_BUILD)
  find_package(MPI REQUIRED)
  add_definitions(-DVERS="5.mpi")
  add_definitions(-DMPI=ON)
  add_definitions(-DDMMP_NAME="mcnp5.mpi")
  add_definitions(-D_DOTCOMM_MPI=ON)
  include_directories(${MPI_INCLUDE_PATH})
  include_directories(${CMAKE_CURRENT_LIST_DIR}/dotcomm/include)
  include_directories(${CMAKE_CURRENT_LIST_DIR}/dotcomm/src/internals/mpi)
else ()
  add_definitions(-DVERS="5")
endif (MPI_BUILD)

# DAGMC
if (DAGMC_BUILD)
  include_directories(${CMAKE_CURRENT_LIST_DIR}/../../tally)
endif (DAGMC_BUILD)

# Source files
include(${CMAKE_CURRENT_LIST_DIR}/src/CMakeLists.txt)
if (MPI_BUILD)
  include(${CMAKE_CURRENT_LIST_DIR}/dotcomm/CMakeLists.txt)
endif (MPI_BUILD)
if (DAGMC_BUILD)
  list(APPEND DAGMC_SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/../mcnp_funcs.cpp)
  list(APPEND DAGMC_SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/../meshtal_funcs.cpp)
  list(APPEND SRC_FILES ${DAGMC_SRC_FILES})
endif (DAGMC_BUILD)

# Add the executable
if (MPI_BUILD)
  set(MCNP_EXEC "mcnp5.mpi")
else ()
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
if (DAGMC_BUILD)
  target_link_libraries(${MCNP_EXEC} ${HDF5_hdf5_hl_LIBRARY})
  target_link_libraries(${MCNP_EXEC} ${HDF5_hdf5_LIBRARY})
  target_link_libraries(${MCNP_EXEC} ${DAG_LINK_LIBS})
  target_link_libraries(${MCNP_EXEC} dagmciface)
  target_link_libraries(${MCNP_EXEC} dagtally)
  target_link_libraries(${MCNP_EXEC} pyne_dagmc)
endif (DAGMC_BUILD)

# Install the targets
install(TARGETS ${MCNP_EXEC} RUNTIME DESTINATION ${INSTALL_BIN_DIR} PERMISSIONS ${PERMS})
