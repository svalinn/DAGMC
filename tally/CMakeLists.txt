message(STATUS "Building dagtally library")

enable_language(Fortran)

# Source files
set(SRC_FILES)
list(APPEND SRC_FILES CellTally.cpp)
list(APPEND SRC_FILES KDEKernel.cpp)
list(APPEND SRC_FILES KDEMeshTally.cpp)
list(APPEND SRC_FILES KDENeighborhood.cpp)
list(APPEND SRC_FILES MeshTally.cpp)
list(APPEND SRC_FILES PolynomialKernel.cpp)
list(APPEND SRC_FILES Quadrature.cpp)
list(APPEND SRC_FILES Tally.cpp)
list(APPEND SRC_FILES TallyData.cpp)
list(APPEND SRC_FILES TallyManager.cpp)
list(APPEND SRC_FILES TrackLengthMeshTally.cpp)

# Public headers
set(PUB_HEADERS)
list(APPEND PUB_HEADERS CellTally.hpp)
list(APPEND PUB_HEADERS KDEKernel.hpp)
list(APPEND PUB_HEADERS KDEMeshTally.hpp)
list(APPEND PUB_HEADERS KDENeighborhood.hpp)
list(APPEND PUB_HEADERS Matrix3.hpp)
list(APPEND PUB_HEADERS MeshTally.hpp)
list(APPEND PUB_HEADERS PolynomialKernel.hpp)
list(APPEND PUB_HEADERS Quadrature.hpp)
list(APPEND PUB_HEADERS TallyData.hpp)
list(APPEND PUB_HEADERS TallyEvent.hpp)
list(APPEND PUB_HEADERS Tally.hpp)
list(APPEND PUB_HEADERS TallyManager.hpp)
list(APPEND PUB_HEADERS TrackLengthMeshTally.hpp)

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS lapack)
list(APPEND LINK_LIBS blas)
list(APPEND LINK_LIBS gfortran)
list(APPEND LINK_LIBS quadmath)

# Include directories
include_directories(.)
include_directories(${MOAB_INCLUDE_DIRS})

# Build dagtally library
add_library(dagtally ${LIB_TYPE} ${SRC_FILES})
set_target_properties(dagtally PROPERTIES PUBLIC_HEADER "${PUB_HEADERS}")
target_link_libraries(dagtally ${LINK_LIBS})
install(TARGETS dagtally LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                         ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
                         PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}
                         LIBRARY       PERMISSIONS ${PERMS}
                         ARCHIVE       PERMISSIONS ${PERMS})
# Build dagtally tests
add_subdirectory(tests)
