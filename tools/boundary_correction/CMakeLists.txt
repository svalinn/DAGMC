message(STATUS "Building KDE boundary correction tool")

# Source files
set(SRC_FILES boundary.cpp)

# Libraries
set(LINK_LIBS ${DAG_LINK_LIBS})

# Include directories
include_directories(${MOAB_INCLUDE_DIRS})

# Build boundary correction tool
add_executable(boundary.exe ${SRC_FILES})
target_link_libraries(boundary.exe ${LINK_LIBS})
install(TARGETS boundary.exe DESTINATION tools
                             PERMISSIONS ${PERMS})
