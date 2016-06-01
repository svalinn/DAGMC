message(STATUS "Building fludag library")

enable_language(Fortran)

# Source files
set(SRC_FILES fluka_funcs.cpp)

# Public headers
set(PUB_HEADERS fluka_funcs.h)

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${FLUKA_DIR}/libflukahp.a)
list(APPEND LINK_LIBS pyne_dagmc)
list(APPEND LINK_LIBS dagmciface)
list(APPEND LINK_LIBS lapack)
list(APPEND LINK_LIBS blas)
list(APPEND LINK_LIBS gfortran)
list(APPEND LINK_LIBS quadmath)

# Include directories
include_directories(.)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Build FluDAG library
add_library(fludag STATIC ${SRC_FILES})
set_target_properties(fludag PROPERTIES PUBLIC_HEADER ${PUB_HEADERS})
target_link_libraries(fludag ${LINK_LIBS})
install(TARGETS fludag LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                       ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
                       PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}
                       LIBRARY       PERMISSIONS ${PERMS})

# Build FluDAG
add_subdirectory(build)

# Build the FluDAG tests
add_subdirectory(tests)
