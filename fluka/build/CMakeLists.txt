message(STATUS "Building FluDAG")

enable_language(Fortran)

# Source files
set(SRC_FILES mainFluDAG.cpp)

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS fludag)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Include directories
include_directories(.)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Build FluDAG
add_executable (mainfludag ${SRC_FILES})
target_link_libraries(mainfludag ${LINK_LIBS})
install(TARGETS mainfludag DESTINATION ${INSTALL_BIN_DIR}
                           PERMISSIONS ${PERMS})
