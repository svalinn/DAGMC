message(STATUS "Including UWUW in dagmciface library")

# libdagmc source files
set(DAGMC_SRC_FILES ${DAGMC_SRC_FILES} ${CMAKE_CURRENT_LIST_DIR}/uwuw.cpp PARENT_SCOPE)

# libdagmc public headers
set(DAGMC_PUB_HEADERS ${DAGMC_PUB_HEADERS} ${CMAKE_CURRENT_LIST_DIR}/uwuw.hpp PARENT_SCOPE)

# Libraries
set(LINK_LIBS_STATIC)
list(APPEND LINK_LIBS_STATIC ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS_STATIC dagmc-static)
list(APPEND LINK_LIBS_STATIC pyne_dagmc-static)

# Include directories
include_directories(.)
include_directories(${CMAKE_SOURCE_DIR}/dagmc)
include_directories(${GTEST_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Build uwuw_preprocessor internal library
set(UWUW_PREPROCESSOR_SRC "uwuw_preprocessor.cpp")
add_library(uwuw_preprocessor-internal STATIC ${UWUW_PREPROCESSOR_SRC})
set_target_properties(uwuw_preprocessor-internal PROPERTIES OUTPUT_NAME uwuw_preprocessor)
target_link_libraries(uwuw_preprocessor-internal ${LINK_LIBS_STATIC})

# Build uwuw_preproc executable
set(UWUW_PREPROC_SRC "uwuw_preproc.cpp")
list(APPEND LINK_LIBS        uwuw_preprocessor-internal)
list(APPEND LINK_LIBS_STATIC uwuw_preprocessor-internal)
add_executable(uwuw_preproc ${UWUW_PREPROC_SRC})
if (BUILD_STATIC)
  target_link_libraries(uwuw_preproc ${LINK_LIBS_STATIC})
else ()
  target_link_libraries(uwuw_preproc ${LINK_LIBS})
endif (BUILD_STATIC)
install(TARGETS uwuw_preproc DESTINATION ${INSTALL_BIN_DIR})

# Build UWUW tests
add_subdirectory(tests)
