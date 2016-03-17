message(STATUS "Building astyle_google")

# Source files
set(CODE_SRC_FILES)
list(APPEND CODE_SRC_FILES "ASBeautifier.cpp")
list(APPEND CODE_SRC_FILES "ASEnhancer.cpp")
list(APPEND CODE_SRC_FILES "ASFormatter.cpp")
list(APPEND CODE_SRC_FILES "ASLocalizer.cpp")
list(APPEND CODE_SRC_FILES "ASResource.cpp")
list(APPEND CODE_SRC_FILES "astyle_main.cpp")

# Headers
set(CODE_PUB_HEADERS)
list(APPEND CODE_PUB_HEADERS "ASLocalizer.h")
list(APPEND CODE_PUB_HEADERS "astyle.h")
list(APPEND CODE_PUB_HEADERS "astyle_main.h")

# Includes
include_directories(.)

# Build astyle_google
add_executable(astyle_google ${CODE_SRC_FILES})
set_target_properties(astyle_google PROPERTIES PUBLIC_HEADER "${CODE_PUB_HEADERS}")
install(TARGETS astyle_google RUNTIME       DESTINATION ${INSTALL_SHARE_DIR}/tools
                              PUBLIC_HEADER DESTINATION ${INSTALL_SHARE_DIR}/tools/include
                              RUNTIME       PERMISSIONS ${PERMS})
