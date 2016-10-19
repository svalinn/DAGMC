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
set(MW_HEADERS)
list(APPEND MW_HEADERS Gen.hpp)
list(APPEND MW_HEADERS Arc.hpp)
list(APPEND MW_HEADERS Zip.hpp)
list(APPEND MW_HEADERS Cleanup.hpp)
list(APPEND MW_HEADERS CheckWatertight.hpp)
list(APPEND MW_HEADERS MakeWatertight.hpp)

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Define executables
add_executable(make_watertight make_watertight.cpp ${SRC_FILES})
add_executable(check_watertight check_watertight.cpp ${SRC_FILES})

# Include directories
include_directories(.)
include_directories(${MOAB_INCLUDE_DIRS})

# Link executables to moab/dagmc libraries
target_link_libraries(make_watertight ${LINK_LIBS})
target_link_libraries(check_watertight ${LINK_LIBS})

#Create library with make_watertight functionality
add_library(makeWatertight ${LIB_TYPE} Zip.cpp Arc.cpp Cleanup.cpp Gen.cpp CheckWatertight.cpp MakeWatertight.cpp)
set_target_properties(makeWatertight PROPERTIES PUBLIC_HEADER "${MW_HEADERS}")

# Install
install(TARGETS check_watertight make_watertight DESTINATION ${INSTALL_BIN_DIR} PERMISSIONS ${PERMS})
install(TARGETS makeWatertight LIBRARY 	     DESTINATION ${INSTALL_LIB_DIR}
			       ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
			       PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}/make_watertight
			       LIBRARY       PERMISSIONS ${PERMS}
			       ARCHIVE 	     PERMISSIONS ${PERMS})


# build tests
add_subdirectory(tests)
