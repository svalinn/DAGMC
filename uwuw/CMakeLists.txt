message(STATUS "Building uwuw_preproc")

# Source files
set(SRC_FILES uwuw_preproc.cpp)

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})
if (BUILD_STATIC)
  list(APPEND LINK_LIBS dagmc-static)
  list(APPEND LINK_LIBS pyne_dagmc-static)
else()
  list(APPEND LINK_LIBS dagmc)
  list(APPEND LINK_LIBS pyne_dagmc)
endif (BUILD_STATIC)

# Include directories
include_directories(.)
include_directories(${CMAKE_SOURCE_DIR}/dagmc)
include_directories(${GTEST_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Build uwuw_preproc executable
add_executable(uwuw_preproc ${SRC_FILES})
target_link_libraries(uwuw_preproc ${LINK_LIBS})
install(TARGETS uwuw_preproc DESTINATION ${INSTALL_BIN_DIR})

# Build UWUW tests
add_subdirectory(tests)
