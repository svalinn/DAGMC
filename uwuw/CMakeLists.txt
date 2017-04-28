message(STATUS "Including UWUW in dagmciface library")

# dagmciface source files
set(DI_SRC_FILES ${DI_SRC_FILES} ${CMAKE_CURRENT_LIST_DIR}/uwuw.cpp PARENT_SCOPE)

# dagmciface public headers
set(DI_PUB_HEADERS ${DI_PUB_HEADERS} ${CMAKE_CURRENT_LIST_DIR}/uwuw.hpp PARENT_SCOPE)

# Libraries
set(LINK_LIBS)
set(LINK_LIBS_STATIC)
list(APPEND LINK_LIBS        ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS        dagmciface)
list(APPEND LINK_LIBS        pyne_dagmc)
list(APPEND LINK_LIBS_STATIC ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS_STATIC dagmciface-static)
list(APPEND LINK_LIBS_STATIC pyne_dagmc-static)

# Include directories
include_directories(..)
include_directories(${CMAKE_SOURCE_DIR}/common)
include_directories(${GTEST_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Build uwuw_preprocessor library
set(UWUW_PREPROCESSOR_SRC "uwuw_preprocessor.cpp")
add_library(uwuw_preprocessor        ${UWUW_PREPROCESSOR_SRC})
add_library(uwuw_preprocessor-static ${UWUW_PREPROCESSOR_SRC})
set_target_properties(uwuw_preprocessor-static PROPERTIES OUTPUT_NAME uwuw_preprocessor)
target_link_libraries(uwuw_preprocessor        ${LINK_LIBS})
target_link_libraries(uwuw_preprocessor-static ${LINK_LIBS_STATIC})

# Build uwuw_preproc executable
set(UWUW_PREPROC_SRC "uwuw_preproc.cpp")
list(APPEND LINK_LIBS        uwuw_preprocessor)
list(APPEND LINK_LIBS_STATIC uwuw_preprocessor-static)
add_executable(uwuw_preproc ${UWUW_PREPROC_SRC})
if (BUILD_STATIC)
  target_link_libraries(uwuw_preproc ${LINK_LIBS_STATIC})
else ()
  target_link_libraries(uwuw_preproc ${LINK_LIBS})
endif (BUILD_STATIC)
install(TARGETS uwuw_preproc DESTINATION ${INSTALL_BIN_DIR})

# Build UWUW tests
add_subdirectory(tests)
