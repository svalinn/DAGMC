message(STATUS "Building build_obb")

# Source files
set(SRC_FILES)
list(APPEND SRC_FILES build_obb.cpp)

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Define executables
add_executable(build_obb ${SRC_FILES})

# Include directories
include_directories(.)
include_directories(${MOAB_INCLUDE_DIRS})

# Link executables to moab/dagmc libraries
target_link_libraries(build_obb ${LINK_LIBS})

# Install
install(TARGETS build_obb DESTINATION ${INSTALL_BIN_DIR} PERMISSIONS ${PERMS})
