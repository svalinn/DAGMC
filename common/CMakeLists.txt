message(STATUS "Building dagmciface library")

set(CODE_SRC_DIR     "common")
set(CODE_SRC_FILES   "dagmciface.cpp")
set(CODE_PUB_HEADERS "dagmciface.h")

set(DI_PUB_HEADERS ${DI_PUB_HEADERS} "${CMAKE_CURRENT_LIST_DIR}/${CODE_PUB_HEADERS}" PARENT_SCOPE)
set(DI_SRC_FILES   ${DI_SRC_FILES}   "${CMAKE_CURRENT_LIST_DIR}/${CODE_SRC_FILES}"   PARENT_SCOPE)

# Build static library
#if (STATIC_LIB)
# set(BUILD_SHARED_LIBRARIES OFF)
# set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# set(CMAKE_EXE_LINKER_FLAGS "-static")
#endif (STATIC_LIB)

add_library(dagmciface ${LIB_TYPE} "${DI_SRC_FILES}")
set_target_properties(dagmciface PROPERTIES PUBLIC_HEADER "${DI_PUB_HEADERS}")
install(TARGETS dagmciface LIBRARY DESTINATION       ${INSTALL_LIB_DIR}
                           ARCHIVE DESTINATION       ${INSTALL_LIB_DIR}
                           PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR})
