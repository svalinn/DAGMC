project(DAGMC)

cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0004 OLD)
if (CMAKE_MAJOR_VERSION GREATER 3 OR CMAKE_MAJOR_VERSION EQUAL 3)
  cmake_policy(SET CMP0040 OLD)
endif (CMAKE_MAJOR_VERSION GREATER 3 OR CMAKE_MAJOR_VERSION EQUAL 3)

# Macro to enable fortran and get rid of some variables
macro (enable_fortran)
  enable_language(Fortran)
  set(CMAKE_Fortran_FLAGS_DEBUG "-g")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g")
  set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "")
  set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "")
  if (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-intel")
  endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
endmacro (enable_fortran)

# Setup test macro
macro (setup_test test_name ext driver libs)
  set(drivers ${driver})
  # Convert driver from string to list
  string(STRIP "${drivers}" drivers)
  string(REGEX REPLACE "[ ]+" " " drivers "${drivers}")
  string(REPLACE " " ";" drivers "${drivers}")
  add_executable(${test_name} ${test_name}.${ext} ${drivers})
  target_link_libraries(${test_name} ${libs})
  install(TARGETS ${test_name} DESTINATION tests)
endmacro (setup_test)

# Installation directories
set(INSTALL_BIN_DIR     bin     CACHE PATH "Installation directory for executables")
set(INSTALL_LIB_DIR     lib     CACHE PATH "Installation directory for libraries")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
set(INSTALL_SHARE_DIR   share   CACHE PATH "Installation directory for other files")

# Default to a release build
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "CMAKE_BUILD_TYPE not specified, defaulting to Release")
  set(CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Release" AND
    NOT CMAKE_BUILD_TYPE STREQUAL "Debug" AND
    NOT CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  message(FATAL_ERROR "Specified CMAKE_BUILD_TYPE is invalid; valid options are Release, Debug, RelWithDebInfo")
endif ()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# Optimization flags: note that these are overwritten for MCNP source code to
# guarantee that nothing higher than -O1 is used
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

# Clear out some variables automatically set by CMake
set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "")
set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "")

if (${CMAKE_SYSTEM_NAME} STREQUAL Darwin AND ${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
  set(CXX_LIBRARY c++)
else ()
  set(CXX_LIBRARY stdc++)
endif (${CMAKE_SYSTEM_NAME} STREQUAL Darwin AND ${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)

# Determine if we're building static executables
option(BUILD_STATIC "Build static executables" OFF)
if (BUILD_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
else ()
  set(BUILD_SHARED ON)
  if (${CMAKE_C_COMPILER_ID} STREQUAL Clang)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  endif (${CMAKE_C_COMPILER_ID} STREQUAL Clang)
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
endif (BUILD_STATIC)

# Allow use of find_package(HDF5), find_package(MOAB), etc.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake)

# find the gfortran library name, this is required for OS/X which uses a variety
# of package managers to install gfortan and caould be in a number of places
if (${CMAKE_SYSTEM_NAME} STREQUAL Darwin)
  find_library(GFORTRAN_LIBRARY NAMES gfortran libgfortran
                                PATHS /usr/local/lib/gcc/* /lib /lib64 /usr/lib
                                      /usr/lib64 /usr/lib/gcc/*/*)
endif (${CMAKE_SYSTEM_NAME} STREQUAL Darwin)

# Make sure we have HDF5
find_package(HDF5 REQUIRED)

# Make sure we have MOAB
find_package(MOAB REQUIRED)
set(DAG_LINK_LIBS ${MOAB_LIBRARIES})
string(STRIP "${DAG_LINK_LIBS}" DAG_LINK_LIBS)
string(REGEX REPLACE ";" " " DAG_LINK_LIBS "${DAG_LINK_LIBS}")
string(REGEX REPLACE "[ ]+" " " DAG_LINK_LIBS "${DAG_LINK_LIBS}")
separate_arguments(DAG_LINK_LIBS)
list(REMOVE_DUPLICATES DAG_LINK_LIBS)
if (BUILD_STATIC)
  list(REMOVE_ITEM DAG_LINK_LIBS z)
  list(APPEND DAG_LINK_LIBS z)
  list(APPEND DAG_LINK_LIBS gfortran)
  list(APPEND DAG_LINK_LIBS quadmath)
endif (BUILD_STATIC)
message(STATUS "DAG_LINK_LIBS: ${DAG_LINK_LIBS}")

# Build Gtest
add_subdirectory(gtest)

# Figure out which physics codes to build
if (BUILD_ALL)
  set(BUILD_MCNP5 ON)
  set(BUILD_MCNP6 ON)
  set(BUILD_GEANT4 ON)
  set(BUILD_FLUKA ON)
endif (BUILD_ALL)
if (BUILD_MCNP5 OR BUILD_MCNP6)
  set(BUILD_TALLY ON)
endif (BUILD_MCNP5 OR BUILD_MCNP6)

# Build pyne_dagmc library
add_subdirectory(pyne)

# Build dagtally library
if (BUILD_TALLY)
  add_subdirectory(tally)
endif (BUILD_TALLY)

# Build MCNP
if (BUILD_MCNP5 OR BUILD_MCNP6)
  set(BUILD_ANY ON)
  add_subdirectory(mcnp)
endif (BUILD_MCNP5 OR BUILD_MCNP6)

# Build Geant4
if (BUILD_GEANT4)
  set(BUILD_ANY ON)
  add_subdirectory(geant4)
endif (BUILD_GEANT4)

# Build FLUKA
if (BUILD_FLUKA)
  set(BUILD_ANY ON)
  add_subdirectory(fluka)
endif (BUILD_FLUKA)

# Include UWUW in DAGMC library
add_subdirectory(uwuw)
add_subdirectory(dagmc) 

# Build tools
add_subdirectory(tools)

# Build tests that wouldn't normally be built for the CI
if (BUILD_CI_TESTS)
  if (NOT BUILD_TALLY)
    add_subdirectory(tally)
  endif (NOT BUILD_TALLY)
  if (NOT BUILD_MCNP5 AND NOT BUILD_MCNP6)
    add_subdirectory(mcnp)
  endif (NOT BUILD_MCNP5 AND NOT BUILD_MCNP6)
endif (BUILD_CI_TESTS)
