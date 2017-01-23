message(STATUS "Building dagtally library")

enable_language(Fortran)

# Source files
set(SRC_FILES)
file(GLOB SRC_FILES ${SRC_FILES} "*.cpp")
file(GLOB SRC_FILES ${SRC_FILES} "lapack/*.f")

# Public headers
set(PUB_HEADERS)
file(GLOB PUB_HEADERS ${PUB_HEADERS} "*.hpp")

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Include directories
include_directories(.)
include_directories(${MOAB_INCLUDE_DIRS})

# Build dagtally library
add_library(dagtally ${LIB_TYPE} ${SRC_FILES})
set_target_properties(dagtally PROPERTIES PUBLIC_HEADER "${PUB_HEADERS}")
target_link_libraries(dagtally ${LINK_LIBS})
install(TARGETS dagtally LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                         ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
                         PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}
                         LIBRARY       PERMISSIONS ${PERMS}
                         ARCHIVE       PERMISSIONS ${PERMS})
# Build dagtally tests
add_subdirectory(tests)
