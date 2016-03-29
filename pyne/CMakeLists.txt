message(STATUS "Building pyne_dagmc library")

# Source files
set(SRC_FILES pyne.cpp)

# Public headers
set(PUB_HEADERS pyne.h)

# Linked libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${HDF5_hdf5_hl_LIBRARY})
list(APPEND LINK_LIBS ${HDF5_hdf5_LIBRARY})

# Include directories
include_directories(.)
include_directories(${HDF5_INCLUDE_DIRS})

# Build pyne_dagmc library
add_library(pyne_dagmc ${LIB_TYPE} ${SRC_FILES})
set_target_properties(pyne_dagmc PROPERTIES PUBLIC_HEADER ${PUB_HEADERS})
target_link_libraries(pyne_dagmc ${LINK_LIBS})
install(TARGETS pyne_dagmc LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                           ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
                           PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}
                           LIBRARY       PERMISSIONS ${PERMS}
                           ARCHIVE       PERMISSIONS ${PERMS})
