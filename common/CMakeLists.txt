message(STATUS "Building dagmciface library")

# Source files
set(DI_SRC_FILES ${DI_SRC_FILES} dagmciface.cpp PARENT_SCOPE)

# Public headers
set(DI_PUB_HEADERS ${DI_PUB_HEADERS} dagmciface.h PARENT_SCOPE)

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
