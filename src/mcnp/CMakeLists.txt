message("")

# Preprocessor directives needed for correct name mangling to be applied
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  add_definitions(-DINTEL_FORTRAN=ON)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_definitions(-DGFORTRAN=ON)
endif ()

include_directories(${CMAKE_SOURCE_DIR}/src/dagmc)
include_directories(${CMAKE_BINARY_DIR}/src/dagmc)
include_directories(${CMAKE_SOURCE_DIR}/src/tally)
include_directories(${CMAKE_SOURCE_DIR}/src/pyne)
include_directories(${CMAKE_SOURCE_DIR}/src/uwuw)

message(STATUS "Building object library: mcnp_funcs")
add_library(mcnp_funcs OBJECT mcnp_funcs.cpp)
message(STATUS "Building object library: meshtal_funcs")
add_library(meshtal_funcs OBJECT meshtal_funcs.cpp)

if (BUILD_MCNP5)
  add_subdirectory(mcnp5)
endif ()
if (BUILD_MCNP6)
  add_subdirectory(mcnp6)
endif ()
if (BUILD_TESTS)
  add_subdirectory(tests)
endif ()
