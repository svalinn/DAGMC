message(STATUS "Building build_obb")

# Source files
set(SRC_FILES)
list(APPEND SRC_FILES build_obb.cpp)

# Include directories
include_directories(.)
include_directories(${MOAB_INCLUDE_DIRS})

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Build build_obb
add_executable(build_obb ${SRC_FILES})
target_link_libraries(build_obb ${LINK_LIBS})
install(TARGETS build_obb DESTINATION ${INSTALL_BIN_DIR} PERMISSIONS ${PERMS})
