message(STATUS "Including MCNP in dagmciface library")

# dagmciface source files
set(DI_SRC_FILES ${DI_SRC_FILES} ${CMAKE_CURRENT_LIST_DIR}/mcnp_funcs.cpp PARENT_SCOPE)

# dagmciface public headers
set(DI_PUB_HEADERS ${DI_PUB_HEADERS} ${CMAKE_CURRENT_LIST_DIR}/mcnp_funcs.h PARENT_SCOPE)

# Build MCNP5
if (BUILD_MCNP5)
  add_subdirectory(mcnp5)
endif (BUILD_MCNP5)

# Build MCNP6
if (BUILD_MCNP6)
  add_subdirectory(mcnp6)
endif (BUILD_MCNP6)

# Build tests
if (BUILD_MCNP5 OR BUILD_MCNP6)
  add_subdirectory(test)
endif (BUILD_MCNP5 OR BUILD_MCNP6)
