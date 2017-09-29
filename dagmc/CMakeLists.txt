message(STATUS "Building DAGMC library")

# Source files
set(SRC_FILES)
list(APPEND SRC_FILES DagMC.cpp)
list(APPEND SRC_FILES dagmcmetadata.cpp)
list(APPEND SRC_FILES ${CMAKE_SOURCE_DIR}/uwuw/uwuw.cpp)
list(APPEND SRC_FILES ${CMAKE_SOURCE_DIR}/uwuw/uwuw_preprocessor.cpp)

# Public headers
set(PUB_HEADERS)
list(APPEND PUB_HEADERS DagMC.hpp)
list(APPEND PUB_HEADERS dagmcmetadata.hpp)
list(APPEND PUB_HEADERS ${CMAKE_SOURCE_DIR}/uwuw/uwuw.hpp)
list(APPEND PUB_HEADERS ${CMAKE_SOURCE_DIR}/uwuw/uwuw_preprocessor.hpp)

# Libraries
set(LINK_LIBS)
set(LINK_LIBS_STATIC)
list(APPEND LINK_LIBS        ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS        pyne_dagmc)
list(APPEND LINK_LIBS_STATIC ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS_STATIC pyne_dagmc-static)

# Include directories
include_directories(.)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Build DAGMC library
add_library(dagmc        SHARED ${SRC_FILES})
add_library(dagmc-static STATIC ${SRC_FILES})
set_target_properties(dagmc        PROPERTIES PUBLIC_HEADER "${PUB_HEADERS}")
set_target_properties(dagmc-static PROPERTIES OUTPUT_NAME dagmc)
install(TARGETS dagmc        LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                             PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR})
install(TARGETS dagmc-static ARCHIVE       DESTINATION ${INSTALL_LIB_DIR})
target_link_libraries(dagmc        ${LINK_LIBS})
target_link_libraries(dagmc-static ${LINK_LIBS_STATIC})

add_subdirectory(tools)
add_subdirectory(test)
