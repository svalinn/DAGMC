message(STATUS "Including UWUW in dagmciface library")

# dagmciface source files
set(DI_SRC_FILES ${DI_SRC_FILES} ${CMAKE_CURRENT_LIST_DIR}/uwuw.cpp PARENT_SCOPE)

# dagmciface public headers
set(DI_PUB_HEADERS ${DI_PUB_HEADERS} ${CMAKE_CURRENT_LIST_DIR}/uwuw.hpp PARENT_SCOPE)

# Include directories
include_directories(..)
include_directories(${CMAKE_SOURCE_DIR}/common)
include_directories(${GTEST_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${HDF5_hdf5_hl_LIBRARY})
list(APPEND LINK_LIBS ${HDF5_hdf5_LIBRARY})
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS dagmciface)
list(APPEND LINK_LIBS pyne_dagmc)

set(UWUW_PREPROCESSOR_SRC "uwuw_preprocessor.cpp")
add_library(uwuw_preprocessor ${UWUW_PREPROCESSOR_SRC})
target_link_libraries(uwuw_preprocessor ${LINK_LIBS})

set(UWUW_PREPROC_SRC "uwuw_preproc.cpp")
list(APPEND LINK_LIBS uwuw_preprocessor)
add_executable(uwuw_preproc ${UWUW_PREPROC_SRC})
target_link_libraries(uwuw_preproc ${LINK_LIBS})

# Build UWUW tests
add_subdirectory(tests)

install(TARGETS uwuw_preproc DESTINATION ${INSTALL_BIN_DIR} PERMISSIONS ${PERMS})
