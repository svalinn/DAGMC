message(STATUS "Building the DagSolid Library for Geant4")

# Find Geant4
if ("${LIB_TYPE}" STREQUAL "STATIC")
  find_package(Geant4 REQUIRED COMPONENTS static)
else ()
  find_package(Geant4 REQUIRED)
endif ("${LIB_TYPE}" STREQUAL "STATIC")

# Setup Geant4 include directories and compile definitions
include(${Geant4_USE_FILE})

# Source files
set(SRC_FILES DagSolid.cc)

# Public headers
set(PUB_HEADERS DagSolid.hh)

# Include directories
include_directories(.)
include_directories(build/include)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})
include_directories(${GEANT4_DIR}/include/Geant4)

# Libraries
set(LINK_LIBS ${Geant4_LIBRARIES})

# Build DagSolid library
add_library(dagsolid ${LIB_TYPE} ${SRC_FILES})
set_target_properties(dagsolid PROPERTIES PUBLIC_HEADER ${PUB_HEADERS})
target_link_libraries(dagsolid ${LINK_LIBS})
install(TARGETS dagsolid LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                         ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
                         PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}
                         LIBRARY       PERMISSIONS ${PERMS}
                         ARCHIVE       PERMISSIONS ${PERMS})

# Build DAG-Geant4
add_subdirectory(build)

# Build the DAG-Geant4 tests
add_subdirectory(tests)
