enable_fortran()

# Preprocessor directives
# Needed for correct name mangling to be applied
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  add_definitions(-DINTEL_FORTRAN=ON)
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  add_definitions(-DGFORTRAN=ON)
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)

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
