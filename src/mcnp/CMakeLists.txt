message("")

# Preprocessor directives needed for correct name mangling to be applied
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  add_definitions(-DINTEL_FORTRAN=ON)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_definitions(-DGFORTRAN=ON)
endif ()

enable_language(Fortran)

# Source files
set(MCNP_FUNCS_SRC)
list(APPEND MCNP_FUNCS_SRC mcnp_funcs.cpp)
list(APPEND MCNP_FUNCS_SRC meshtal_funcs.cpp)
list(APPEND MCNP_FUNCS_SRC source_sampling.cpp)
list(APPEND MCNP_FUNCS_SRC measure.cpp)

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

include_directories(${CMAKE_SOURCE_DIR}/src/dagmc)
include_directories(${CMAKE_BINARY_DIR}/src/dagmc)
include_directories(${CMAKE_SOURCE_DIR}/src/tally)

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
