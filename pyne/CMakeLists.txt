message(STATUS "Building pyne_dagmc library")

# Source files
set(SRC_FILES pyne.cpp)

# Public headers
set(PUB_HEADERS pyne.h)

# Linked libraries
set(LINK_LIBS        ${DAG_LINK_LIBS})
set(LINK_LIBS_STATIC ${DAG_LINK_LIBS})

# Include directories
include_directories(.)
include_directories(${HDF5_INCLUDE_DIRS})

# Build pyne_dagmc library
add_library(pyne_dagmc        SHARED ${SRC_FILES})
add_library(pyne_dagmc-static STATIC ${SRC_FILES})
set_target_properties(pyne_dagmc        PROPERTIES PUBLIC_HEADER "${PUB_HEADERS}")
set_target_properties(pyne_dagmc-static PROPERTIES OUTPUT_NAME pyne_dagmc)
target_link_libraries(pyne_dagmc        ${LINK_LIBS})
target_link_libraries(pyne_dagmc-static ${LINK_LIBS_STATIC})
install(TARGETS pyne_dagmc        LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                                  PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR})
install(TARGETS pyne_dagmc-static ARCHIVE       DESTINATION ${INSTALL_LIB_DIR})
