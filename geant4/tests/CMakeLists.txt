message(STATUS "Building DagSolid unit tests")

# CXX compiler flags
if (BUILD_STATIC)
  # Needed for geant4 v10.02 and earlier
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--allow-multiple-definition")
endif (BUILD_STATIC)

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
include_directories(${CMAKE_SOURCE_DIR}/geant4)
include_directories(${GTEST_INCLUDE_DIR})

# Setup unit tests
set(DRIVER dagsolid_unit_test_driver.cpp)
setup_test(dagsolid_unit_tests cpp ${DRIVER} "${LINK_LIBS}")

# Install h5m files
install(FILES test_geom.h5m DESTINATION tests)
