message(STATUS "Building pyne_dagmc library")

set(CODE_SRC_DIR     "pyne")
set(CODE_SRC_FILES   "pyne.cpp")
set(CODE_PUB_HEADERS "pyne.h")

# Includes
include_directories(.)

# Build static library
#if (STATIC_LIB)
# set(BUILD_SHARED_LIBRARIES OFF)
# set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# set(CMAKE_EXE_LINKER_FLAGS "-static")
#endif (STATIC_LIB)

# Build pyne_dagmc library
add_library(pyne_dagmc ${LIB_TYPE} "${CODE_SRC_FILES}")
set_target_properties(pyne_dagmc PROPERTIES PUBLIC_HEADER "${CODE_PUB_HEADERS}")
target_link_libraries(pyne_dagmc ${HDF5_hdf5_hl_LIBRARY})
target_link_libraries(pyne_dagmc ${HDF5_hdf5_LIBRARY})
install(TARGETS pyne_dagmc LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                           ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
                           PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR})
