message(STATUS "Building dagmciface library")

# Source files
list(APPEND DI_SRC_FILES dagmcmetadata.cpp)

# Public headers
list(APPEND DI_PUB_HEADERS dagmcmetadata.hpp)

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

# Build dagmciface library
add_library(dagmciface        SHARED ${DI_SRC_FILES})
add_library(dagmciface-static STATIC ${DI_SRC_FILES})
set_target_properties(dagmciface        PROPERTIES PUBLIC_HEADER "${DI_PUB_HEADERS}")
set_target_properties(dagmciface-static PROPERTIES OUTPUT_NAME dagmciface)
install(TARGETS dagmciface        LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                                  PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR})
install(TARGETS dagmciface-static ARCHIVE       DESTINATION ${INSTALL_LIB_DIR})
target_link_libraries(dagmciface        ${LINK_LIBS})
target_link_libraries(dagmciface-static ${LINK_LIBS_STATIC})

add_subdirectory(test)
