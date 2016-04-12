message(STATUS "Including MCNP5 in dagmciface library")

# dagmciface source files
set(DI_SRC_FILES ${DI_SRC_FILES} ${CMAKE_CURRENT_LIST_DIR}/mcnp_funcs.cpp PARENT_SCOPE)

# dagmciface public headers
set(DI_PUB_HEADERS ${DI_PUB_HEADERS} ${CMAKE_CURRENT_LIST_DIR}/mcnp_funcs.h PARENT_SCOPE)

# Build DAG-MCNP5
add_subdirectory(build)
