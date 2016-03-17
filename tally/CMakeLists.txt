message(STATUS "Building dagtally library")

enable_language(Fortran)

set(CODE_SRC_DIR "tally")

# Source files
list(APPEND CODE_SRC_FILES "CellTally.cpp")
list(APPEND CODE_SRC_FILES "KDEKernel.cpp")
list(APPEND CODE_SRC_FILES "KDEMeshTally.cpp")
list(APPEND CODE_SRC_FILES "KDENeighborhood.cpp")
list(APPEND CODE_SRC_FILES "MeshTally.cpp")
list(APPEND CODE_SRC_FILES "PolynomialKernel.cpp")
list(APPEND CODE_SRC_FILES "Quadrature.cpp")
list(APPEND CODE_SRC_FILES "Tally.cpp")
list(APPEND CODE_SRC_FILES "TallyData.cpp")
list(APPEND CODE_SRC_FILES "TallyManager.cpp")
list(APPEND CODE_SRC_FILES "TrackLengthMeshTally.cpp")

# Headers
list(APPEND CODE_PUB_HEADERS "CellTally.hpp")
list(APPEND CODE_PUB_HEADERS "KDEKernel.hpp")
list(APPEND CODE_PUB_HEADERS "KDEMeshTally.hpp")
list(APPEND CODE_PUB_HEADERS "KDENeighborhood.hpp")
list(APPEND CODE_PUB_HEADERS "Matrix3.hpp")
list(APPEND CODE_PUB_HEADERS "MeshTally.hpp")
list(APPEND CODE_PUB_HEADERS "PolynomialKernel.hpp")
list(APPEND CODE_PUB_HEADERS "Quadrature.hpp")
list(APPEND CODE_PUB_HEADERS "TallyData.hpp")
list(APPEND CODE_PUB_HEADERS "TallyEvent.hpp")
list(APPEND CODE_PUB_HEADERS "Tally.hpp")
list(APPEND CODE_PUB_HEADERS "TallyManager.hpp")
list(APPEND CODE_PUB_HEADERS "TrackLengthMeshTally.hpp")

# Includes
include_directories(.)

# Build static library
#if (STATIC_LIB)
# set(BUILD_SHARED_LIBRARIES OFF)
# set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# set(CMAKE_EXE_LINKER_FLAGS "-static")
#endif (STATIC_LIB)

# Build dagtally library
add_library(dagtally ${LIB_TYPE} "${CODE_SRC_FILES}")
set_target_properties(dagtally PROPERTIES PUBLIC_HEADER "${CODE_PUB_HEADERS}")
target_link_libraries(dagtally gfortran)
target_link_libraries(dagtally lapack)
target_link_libraries(dagtally blas)
target_link_libraries(dagtally ${DAG_LINK_LIBS})
install(TARGETS dagtally LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                         ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
                         PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR})

# Build dagtally tests
add_subdirectory(tests)
