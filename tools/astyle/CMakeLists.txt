message(STATUS "Building astyle_google")

# Source files
set(SRC_FILES ASBeautifier.cpp
              ASEnhancer.cpp
              ASFormatter.cpp
              ASLocalizer.cpp
              ASResource.cpp
              astyle_main.cpp
)

# Public headers
set(PUB_HEADERS ASLocalizer.h
                astyle.h
                astyle_main.h
)

# Includes
include_directories(.)

# Build astyle_google
add_executable(astyle_google ${SRC_FILES})
set_target_properties(astyle_google PROPERTIES PUBLIC_HEADER ${PUB_HEADERS})
install(TARGETS astyle_google RUNTIME       DESTINATION ${INSTALL_SHARE_DIR}/tools
                              PUBLIC_HEADER DESTINATION ${INSTALL_SHARE_DIR}/tools/include
                              RUNTIME       PERMISSIONS ${PERMS})
