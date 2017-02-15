project(DAGMC)

cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0004 OLD)
if (CMAKE_MAJOR_VERSION GREATER 3 OR CMAKE_MAJOR_VERSION EQUAL 3)
  cmake_policy(SET CMP0040 OLD)
endif (CMAKE_MAJOR_VERSION GREATER 3 OR CMAKE_MAJOR_VERSION EQUAL 3)

# Setup test macro
macro (setup_test test_name ext driver libs)
  add_executable(${test_name} ${test_name}.${ext} ${driver})
  target_link_libraries(${test_name} ${libs})
  install(TARGETS ${test_name} DESTINATION tests
                               PERMISSIONS ${PERMS})
  add_test(${test_name} ${test_name})
endmacro (setup_test)

# Run test macro
macro(run_test test_target)
  add_custom_target(${test_target}_runtest COMMAND ${test_target}
                                           DEPENDS ${test_target}
                                           WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  add_dependencies(${test_target}_runtest gtest)
endmacro (run_test)

# Installation directories
set(INSTALL_BIN_DIR     bin     CACHE PATH "Installation directory for executables")
set(INSTALL_LIB_DIR     lib     CACHE PATH "Installation directory for libraries")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
set(INSTALL_SHARE_DIR   share   CACHE PATH "Installation directory for other files")

# Permissions of compiled binaries and libraries
set(PERMS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE)

# Determine if we're building static libraries
if (STATIC_LIB)
  set(LIB_TYPE "STATIC")
  set(BUILD_SHARED_LIBRARIES OFF)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
else ()
  set(LIB_TYPE "SHARED")
endif (STATIC_LIB)

# Allow use of find_package(HDF5), find_package(MOAB), etc.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake)


# Make sure we have HDF5 and MOAB
find_package(HDF5 REQUIRED)
find_package(MOAB REQUIRED)

# find the gfortran library name, this is required for OS/X which uses a variety
# of package managers to install gfortan and caould be in a number of places
find_library(GFORTRAN_LIBRARY NAMES gfortran libgfortran
  PATHS /usr/local/lib/gcc/*/ /lib /lib64 /usr/lib /usr/lib64 /usr/lib/gcc/*/*)

# DAGMC_LIBRARIES used instead of MOAB_LIBRARIES for MOAB 4.7+
if (DAGMC_LIBRARIES)
  set(DAG_LINK_LIBS "${DAGMC_LIBRARIES}")
else ()
  set(DAG_LINK_LIBS "dagmc ${MOAB_LIBRARIES}")
endif (DAGMC_LIBRARIES)

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

# Include UWUW in dagmciface library
add_subdirectory(uwuw)
add_subdirectory(common) 

# Build tools
add_subdirectory(tools)

# Build tests that wouldn't normally be built for the CI
if (BUILD_TESTS)
  if (NOT BUILD_TALLY)
    add_subdirectory(tally)
  endif (NOT BUILD_TALLY)
  if (NOT BUILD_MCNP5 AND NOT BUILD_MCNP6)
    add_subdirectory(mcnp/test)
  endif (NOT BUILD_MCNP5 AND NOT BUILD_MCNP6)
endif (BUILD_TESTS)
