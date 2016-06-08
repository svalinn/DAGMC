message(STATUS "Building DagSolid unit tests")

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS gtest)
list(APPEND LINK_LIBS pthread)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS dagsolid)
list(APPEND LINK_LIBS ${Geant4_LIBRARIES})

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
