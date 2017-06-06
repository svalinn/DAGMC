message(STATUS "Building dagtally library")

enable_language(Fortran)

# Source files
set(SRC_FILES)
file(GLOB SRC_FILES ${SRC_FILES} "*.cpp")

# Public headers
set(PUB_HEADERS)
file(GLOB PUB_HEADERS ${PUB_HEADERS} "*.hpp")

# Libraries
set(LINK_LIBS)
set(LINK_LIBS_STATIC)
list(APPEND LINK_LIBS        ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS_STATIC ${DAG_LINK_LIBS})

# Include directories
include_directories(.)
include_directories(${MOAB_INCLUDE_DIRS})

# CXX compiler flags
if (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)

# Build dagtally library
add_library(dagtally        SHARED ${SRC_FILES})
add_library(dagtally-static STATIC ${SRC_FILES})
set_target_properties(dagtally        PROPERTIES PUBLIC_HEADER "${PUB_HEADERS}")
set_target_properties(dagtally-static PROPERTIES OUTPUT_NAME dagtally)
target_link_libraries(dagtally ${LINK_LIBS})
target_link_libraries(dagtally-static ${LINK_LIBS_STATIC})
install(TARGETS dagtally        LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                                PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR})
install(TARGETS dagtally-static ARCHIVE       DESTINATION ${INSTALL_LIB_DIR})

# Build dagtally tests
add_subdirectory(tests)
