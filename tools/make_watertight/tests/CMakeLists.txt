message(STATUS "Building make_watertight tests")

# Source files
set(SRC_FILES)
list(APPEND SRC_FILES ../gen.cpp)
list(APPEND SRC_FILES ../arc.cpp)
list(APPEND SRC_FILES ../zip.cpp)
list(APPEND SRC_FILES ../cleanup.cpp)
list(APPEND SRC_FILES ../mw_func.cpp)
list(APPEND SRC_FILES ../cw_func.cpp)
list(APPEND SRC_FILES test_funcs.cpp)

# Public headers
set(PUB_HEADERS)
list(APPEND PUB_HEADERS ../gen.hpp)
list(APPEND PUB_HEADERS ../arc.hpp)
list(APPEND PUB_HEADERS ../zip.hpp)
list(APPEND PUB_HEADERS ../cleanup.hpp)
list(APPEND PUB_HEADERS ../cw_func.hpp)
list(APPEND PUB_HEADERS ../mw_func.hpp)
list(APPEND PUB_HEADERS test_funcs.hpp)

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Define executables
add_executable(test_cyl test_cyl.cpp ${SRC_FILES})
add_executable(test_cones test_cones.cpp ${SRC_FILES})

# Includes
include_directories(..)
include_directories(${MOAB_INCLUDE_DIRS})

# Link to libraries
target_link_libraries(test_cyl ${LINK_LIBS})
target_link_libraries(test_cones ${LINK_LIBS})

# Install tests and test models
install(TARGETS test_cyl test_cones DESTINATION tests PERMISSIONS ${PERMS})
install(FILES cones.h5m cyl.h5m DESTINATION tests)