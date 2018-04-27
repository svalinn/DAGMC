message(STATUS "Building make_watertight")

# Source files
set(SRC_FILES)
list(APPEND SRC_FILES Gen.cpp)
list(APPEND SRC_FILES Arc.cpp)
list(APPEND SRC_FILES Zip.cpp)
list(APPEND SRC_FILES Cleanup.cpp)
list(APPEND SRC_FILES MakeWatertight.cpp)
list(APPEND SRC_FILES CheckWatertight.cpp)

# Public headers
set(PUB_HEADERS)
list(APPEND PUB_HEADERS Gen.hpp)
list(APPEND PUB_HEADERS Arc.hpp)
list(APPEND PUB_HEADERS Zip.hpp)
list(APPEND PUB_HEADERS Cleanup.hpp)
list(APPEND PUB_HEADERS CheckWatertight.hpp)
list(APPEND PUB_HEADERS MakeWatertight.hpp)

# Libraries
set(LINK_LIBS ${DAG_LINK_LIBS})

# Define executables
add_executable(make_watertight  make_watertight.cpp  ${SRC_FILES})
add_executable(check_watertight check_watertight.cpp ${SRC_FILES})

# Include directories
include_directories(.)
include_directories(${MOAB_INCLUDE_DIRS})

# Build library with make_watertight functionality
add_library(makeWatertight        SHARED ${SRC_FILES})
add_library(makeWatertight-static STATIC ${SRC_FILES})
set_target_properties(makeWatertight        PROPERTIES PUBLIC_HEADER "${PUB_HEADERS}")
set_target_properties(makeWatertight-static PROPERTIES OUTPUT_NAME makeWatertight)
target_link_libraries(makeWatertight        ${LINK_LIBS})
target_link_libraries(makeWatertight-static ${LINK_LIBS})
install(TARGETS makeWatertight        LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                                      PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}/make_watertight)
install(TARGETS makeWatertight-static ARCHIVE       DESTINATION ${INSTALL_LIB_DIR})

# Build executables
target_link_libraries(make_watertight ${LINK_LIBS})
target_link_libraries(check_watertight ${LINK_LIBS})
install(TARGETS make_watertight  DESTINATION ${INSTALL_BIN_DIR})
install(TARGETS check_watertight DESTINATION ${INSTALL_BIN_DIR})

# build tests
add_subdirectory(tests)
