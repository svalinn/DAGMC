

# CMAKE build script for DAGMC Tally
SET ( DAG_TALLY_FILES 
    Tally.cpp
    TallyData.cpp
    TallyManager.cpp
    MeshTally.cpp
    TrackLengthMeshTally.cpp
    KDEMeshTally.cpp
    CellTally.cpp
    KDEKernel.cpp
    KDENeighborhood.cpp
    PolynomialKernel.cpp
    Quadrature.cpp
)

IF ( BUILD_MCNP5 OR BUILD_ALL)
  SET ( DAG_TALLY_FILES ${DAG_TALLY_FILES}
        ${DAG_TALLY_SRCDIR}/mcnp5/meshtal_funcs.cpp )
  INCLUDE_DIRECTORIES(${DAG_TALLY_SRCDIR}/mcnp5)
ENDIF ()


# Create DAGMC Tally library
ADD_LIBRARY(dagtally SHARED ${DAG_TALLY_FILES} )

TARGET_LINK_LIBRARIES(dagtally lapack)
