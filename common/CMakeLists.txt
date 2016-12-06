message(STATUS "Building dagmciface library")

# Source files
list(APPEND DI_SRC_FILES dagmcmetadata.cpp)
                                
# Public headers
list(APPEND DI_PUB_HEADERS dagmcmetadata.hpp)

# Include directories
include_directories(.)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Build dagmciface library
add_library(dagmciface ${LIB_TYPE} ${DI_SRC_FILES})
set_target_properties(dagmciface PROPERTIES PUBLIC_HEADER "${DI_PUB_HEADERS}")
install(TARGETS dagmciface LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                           ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
                           PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}
                           LIBRARY       PERMISSIONS ${PERMS}
                           ARCHIVE       PERMISSIONS ${PERMS})
list(APPEND LINK_LIBS pyne_dagmc)
list(APPEND LINK_LIBS ${HDF5_hdf5_hl_LIBRARY})
list(APPEND LINK_LIBS ${HDF5_hdf5_LIBRARY})
list(APPEND LINK_LIBS ${DAGMC_LIBRARIES})
target_link_libraries(dagmciface ${LINK_LIBS})
add_subdirectory(test)
