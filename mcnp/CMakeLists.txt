# Source files
set(MCNP_FUNCS_SRC)
list(APPEND MCNP_FUNCS_SRC mcnp_funcs.cpp)
list(APPEND MCNP_FUNCS_SRC meshtal_funcs.cpp)

# Libraries
set(LINK_LIBS_STATIC)
list(APPEND LINK_LIBS_STATIC ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS_STATIC dagmc-static)
list(APPEND LINK_LIBS_STATIC dagtally-static)
if (${CMAKE_SYSTEM_NAME} STREQUAL Darwin)
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
    list(APPEND LINK_LIBS_STATIC c++)
  else ()
    list(APPEND LINK_LIBS_STATIC stdc++)
  endif (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
endif (${CMAKE_SYSTEM_NAME} STREQUAL Darwin)

# Include directories
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/dagmc)
include_directories(${CMAKE_SOURCE_DIR}/tally)

# Build mcnpfuncs internal library
add_library(mcnpfuncs-internal STATIC ${MCNP_FUNCS_SRC})
set_target_properties(mcnpfuncs-internal PROPERTIES OUTPUT_NAME mcnpfuncs)
target_link_libraries(mcnpfuncs-internal ${LINK_LIBS_STATIC})

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
