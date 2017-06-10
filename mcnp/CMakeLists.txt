message(STATUS "Building mcnpfuncs internal library")

enable_fortran()

# Source files
set(MCNP_FUNCS_SRC)
list(APPEND MCNP_FUNCS_SRC mcnp_funcs.cpp)
list(APPEND MCNP_FUNCS_SRC meshtal_funcs.cpp)

# Libraries
set(LINK_LIBS_STATIC)
list(APPEND LINK_LIBS_STATIC ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS_STATIC dagmc-static)
list(APPEND LINK_LIBS_STATIC dagtally-static)
if (${CMAKE_SYSTEM_NAME} STREQUAL Darwin AND ${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
  list(APPEND LINK_LIBS_STATIC c++)
else ()
  list(APPEND LINK_LIBS_STATIC stdc++)
endif (${CMAKE_SYSTEM_NAME} STREQUAL Darwin AND ${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)

# Include directories
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/dagmc)
include_directories(${CMAKE_SOURCE_DIR}/tally)

# Preprocessor directives
# Needed for correct name mangling to be applied
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  add_definitions(-DINTEL_FORTRAN=ON)
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  add_definitions(-DGFORTRAN=ON)
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)

# Build mcnpfuncs object library
add_library(mcnpfuncs STATIC ${MCNP_FUNCS_SRC})
target_link_libraries(mcnpfuncs ${LINK_LIBS_STATIC})

# Build MCNP5
if (BUILD_MCNP5)
  add_subdirectory(mcnp5)
endif (BUILD_MCNP5)

# Build MCNP6
if (BUILD_MCNP6)
  add_subdirectory(mcnp6)
endif (BUILD_MCNP6)

# Build tests
if (BUILD_MCNP5 OR BUILD_MCNP6 OR BUILD_CI_TESTS)
  add_subdirectory(test)
endif (BUILD_MCNP5 OR BUILD_MCNP6 OR BUILD_CI_TESTS)
