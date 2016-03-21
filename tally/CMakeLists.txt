message(STATUS "Building dagtally library")

enable_language(Fortran)

# Source files
set(SRC_FILES CellTally.cpp
              KDEKernel.cpp
              KDEMeshTally.cpp
              KDENeighborhood.cpp
              MeshTally.cpp
              PolynomialKernel.cpp
              Quadrature.cpp
              Tally.cpp
              TallyData.cpp
              TallyManager.cpp
              TrackLengthMeshTally.cpp
)

# Public headers
set(PUB_HEADERS CellTally.hpp
                KDEKernel.hpp
                KDEMeshTally.hpp
                KDENeighborhood.hpp
                Matrix3.hpp
                MeshTally.hpp
                PolynomialKernel.hpp
                Quadrature.hpp
                TallyData.hpp
                TallyEvent.hpp
                Tally.hpp
                TallyManager.hpp
                TrackLengthMeshTally.hpp
)

# Libraries
set(LINK_LIBS gfortran
              blas
              lapack
              ${DAG_LINK_LIBS}
)

# Includes
include_directories(.)
include_directories(${MOAB_INCLUDE_DIRS})

# Build dagtally library
add_library(dagtally ${LIB_TYPE} ${SRC_FILES})
set_target_properties(dagtally PROPERTIES PUBLIC_HEADER ${PUB_HEADERS})
target_link_libraries(dagtally ${LINK_LIBS})
install(TARGETS dagtally LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                         ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
                         PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}
                         LIBRARY       PERMISSIONS ${PERMS}
)

# Build dagtally tests
add_subdirectory(tests)
