project(DAGMC)

cmake_minimum_required(VERSION 2.8)

set(DAGMC_MAJOR_VERSION 3)
set(DAGMC_MINOR_VERSION 0)
set(DAGMC_PATCH_VERSION 0)
set(DAGMC_VERSION
  ${DAGMC_MAJOR_VERSION}.${DAGMC_MINOR_VERSION}.${DAGMC_PATCH_VERSION})
message(STATUS "Configuring DAGMC ${DAGMC_VERSION}")

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
  add_test(${test_name} ${test_name})
endmacro (setup_test)

# Run test macro
macro(run_test test_target)
  add_custom_target(${test_target}_runtest COMMAND ${test_target}
                                           DEPENDS ${test_target}
                                           WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  add_dependencies(${test_target}_runtest gtest)
endmacro (run_test)

# Installation directories
set(INSTALL_BIN_DIR     bin     CACHE PATH "Installation directory for executables")
set(INSTALL_LIB_DIR     lib     CACHE PATH "Installation directory for libraries")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
set(INSTALL_SHARE_DIR   share   CACHE PATH "Installation directory for other files")

set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "")
set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "")

# Determine if we're building static executables
option(BUILD_STATIC "Build static executables" OFF)
if (BUILD_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
  endif (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
else ()
  set(BUILD_SHARED ON)
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
# some cleanup of the MOAB libraries variable
string(STRIP "${MOAB_LIBRARIES}" MOAB_LIBRARIES)
string(REGEX REPLACE ";" " " MOAB_LIBRARIES "${MOAB_LIBRARIES}")
string(REGEX REPLACE "[ ]+" " " MOAB_LIBRARIES "${MOAB_LIBRARIES}")
separate_arguments(MOAB_LIBRARIES)
list(REMOVE_DUPLICATES MOAB_LIBRARIES)
if (BUILD_STATIC)
  list(REMOVE_ITEM MOAB_LIBRARIES z)
  list(APPEND MOAB_LIBRARIES z)
  list(APPEND MOAB_LIBRARIES gfortran)
  list(APPEND MOAB_LIBRARIES quadmath)
endif (BUILD_STATIC)
# separate the MOAB library directory from the DAG_LINK_LIBS variable
list(GET MOAB_LIBRARIES 0 MOAB_LIBRARY_DIRS) # ASSUMES THIS DIRECTORY COMES FIRST
string(REGEX REPLACE "-L" "" MOAB_LIBRARY_DIRS "${MOAB_LIBRARY_DIRS}")

set(DAG_LINK_LIBS ${MOAB_LIBRARIES})
message(STATUS "DAG_LINK_LIBS: ${DAG_LINK_LIBS}")

# Get some environment variables
set(ENV_USER "$ENV{USER}")
execute_process(COMMAND hostname       OUTPUT_VARIABLE ENV_HOST OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND uname -s       OUTPUT_VARIABLE ENV_OS   OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND date +%m/%d/%y OUTPUT_VARIABLE ENV_DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND date +%H:%M:%S OUTPUT_VARIABLE ENV_TIME OUTPUT_STRIP_TRAILING_WHITESPACE)

# Build things
add_subdirectory(src)

# Setup the configuration file and install
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/DAGMCConfig.cmake.in DAGMCConfig.cmake @ONLY)
install( FILES ${CMAKE_CURRENT_BINARY_DIR}/DAGMCConfig.cmake DESTINATION ${INSTALL_LIB_DIR}/cmake/)
