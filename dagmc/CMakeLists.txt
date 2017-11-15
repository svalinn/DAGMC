message(STATUS "Building DAGMC library")

# Configure version header
configure_file(DagMCVersion.hpp.in DagMCVersion.hpp)
list(APPEND DAGMC_PUB_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/DagMCVersion.hpp)

# Source files
list(APPEND DAGMC_SRC_FILES DagMC.cpp)
list(APPEND DAGMC_SRC_FILES dagmcmetadata.cpp)
list(APPEND DAGMC_SRC_FILES sdf.cpp)

# Public headers
list(APPEND DAGMC_PUB_HEADERS DagMC.hpp)
list(APPEND DAGMC_PUB_HEADERS dagmcmetadata.hpp)

# Libraries
set(LINK_LIBS)
set(LINK_LIBS_STATIC)
list(APPEND LINK_LIBS        ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS        pyne_dagmc)
list(APPEND LINK_LIBS_STATIC ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS_STATIC pyne_dagmc-static)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# CXX compiler flags
if (BUILD_SHARED AND ${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif (BUILD_SHARED AND ${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)

# Build DAGMC library
add_library(dagmc        SHARED ${DAGMC_SRC_FILES})
add_library(dagmc-static STATIC ${DAGMC_SRC_FILES})
set_target_properties(dagmc        PROPERTIES PUBLIC_HEADER "${DAGMC_PUB_HEADERS}")
set_target_properties(dagmc-static PROPERTIES OUTPUT_NAME dagmc)
install(TARGETS dagmc        LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                             PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR})
install(TARGETS dagmc-static ARCHIVE       DESTINATION ${INSTALL_LIB_DIR})
target_link_libraries(dagmc        ${LINK_LIBS})
target_link_libraries(dagmc-static ${LINK_LIBS_STATIC})

add_subdirectory(tools)
add_subdirectory(test)
