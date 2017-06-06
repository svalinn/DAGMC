message(STATUS "Building DagSolid unit tests")

# Linker flags
if (BUILD_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
endif (BUILD_STATIC)

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS pthread)
list(APPEND LINK_LIBS gtest)
if (BUILD_STATIC)
  list(APPEND LINK_LIBS ${DAG_LINK_LIBS})
  list(APPEND LINK_LIBS dagsolid-static)
else ()
  list(APPEND LINK_LIBS ${DAG_LINK_LIBS})
  list(APPEND LINK_LIBS dagsolid)
endif (BUILD_STATIC)

# Include directories
include_directories(..)
include_directories(${GTEST_INCLUDE_DIR})

# Setup unit tests
set(DRIVER dagsolid_unit_test_driver.cpp)
setup_test(dagsolid_unit_tests cpp ${DRIVER} "${LINK_LIBS}")

# Install h5m files
install(FILES test_geom.h5m DESTINATION tests)

# Enable unit tests
enable_testing()
