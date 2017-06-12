message(STATUS "Building fludag library")

enable_language(Fortran)

# Source files
set(SRC_FILES fluka_funcs.cpp)

# Public headers
set(PUB_HEADERS fluka_funcs.h)

# Libraries
set(LINK_LIBS)
set(LINK_LIBS_STATIC)
list(APPEND LINK_LIBS        ${FLUKA_DIR}/libflukahp.a)
list(APPEND LINK_LIBS        dagmc)
list(APPEND LINK_LIBS        pyne_dagmc)
list(APPEND LINK_LIBS        gfortran)
list(APPEND LINK_LIBS_STATIC ${FLUKA_DIR}/libflukahp.a)
list(APPEND LINK_LIBS_STATIC dagmc-static)
list(APPEND LINK_LIBS_STATIC pyne_dagmc-static)
list(APPEND LINK_LIBS_STATIC gfortran)
list(APPEND LINK_LIBS_STATIC quadmath)

# Include directories
include_directories(.)
include_directories(${CMAKE_SOURCE_DIR}/dagmc)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Build FluDAG library
add_library(fludag        SHARED ${SRC_FILES})
add_library(fludag-static STATIC ${SRC_FILES})
set_target_properties(fludag        PROPERTIES PUBLIC_HEADER ${PUB_HEADERS})
set_target_properties(fludag-static PROPERTIES OUTPUT_NAME fludag)
target_link_libraries(fludag        ${LINK_LIBS})
target_link_libraries(fludag-static ${LINK_LIBS_STATIC})
install(TARGETS fludag        LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                              PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR})
install(TARGETS fludag-static ARCHIVE       DESTINATION ${INSTALL_LIB_DIR})

# Build FluDAG
add_subdirectory(build)

# Build the FluDAG tests
add_subdirectory(tests)
