message(STATUS "Building astyle_google")

# Source files
set(SRC_FILES)
list(APPEND SRC_FILES ASBeautifier.cpp)
list(APPEND SRC_FILES ASEnhancer.cpp)
list(APPEND SRC_FILES ASFormatter.cpp)
list(APPEND SRC_FILES ASLocalizer.cpp)
list(APPEND SRC_FILES ASResource.cpp)
list(APPEND SRC_FILES astyle_main.cpp)

# Include directories
include_directories(.)

# Build astyle_google
add_executable(astyle_google ${SRC_FILES})
install(TARGETS astyle_google DESTINATION tools
                              PERMISSIONS ${PERMS})
