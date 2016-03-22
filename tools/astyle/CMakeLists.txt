message(STATUS "Building astyle_google")

# Source files
set(SRC_FILES)
list(APPEND SRC_FILES ASBeautifier.cpp)
list(APPEND SRC_FILES ASEnhancer.cpp)
list(APPEND SRC_FILES ASFormatter.cpp)
list(APPEND SRC_FILES ASLocalizer.cpp)
list(APPEND SRC_FILES ASResource.cpp)
list(APPEND SRC_FILES astyle_main.cpp)

# Public headers
set(PUB_HEADERS)
list(APPEND PUB_HEADERS ASLocalizer.h)
list(APPEND PUB_HEADERS astyle.h)
list(APPEND PUB_HEADERS astyle_main.h)

# Includes
include_directories(.)

# Build astyle_google
add_executable(astyle_google ${SRC_FILES})
set_target_properties(astyle_google PROPERTIES PUBLIC_HEADER ${PUB_HEADERS})
install(TARGETS astyle_google RUNTIME       DESTINATION ${INSTALL_SHARE_DIR}/tools
                              PUBLIC_HEADER DESTINATION ${INSTALL_SHARE_DIR}/tools/include
                              RUNTIME       PERMISSIONS ${PERMS})
