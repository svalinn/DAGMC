message(STATUS "Building KDE boundary correction tool")

# Source files
set(SRC_FILES boundary.cpp)

# Libraries
set(LINK_LIBS ${MOAB_DIR}/lib/libMOAB.so)
set(LINK_LIBS ${MOAB_DIR}/lib/libdagmc.so)

# Include directories
include_directories(${MOAB_DIR}/include)

# Build boundary correction tool
add_executable(boundary.exe ${SRC_FILES})
target_link_libraries(boundary.exe ${LINK_LIBS})
install(TARGETS boundary.exe DESTINATION ${INSTALL_SHARE_DIR}/tools
                             PERMISSIONS ${PERMS})
