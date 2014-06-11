# CMAKE build script for DAGMC Tally
SET ( DAG_TALLY_SRCDIR "tally" )

# Create DAGMC Tally library
ADD_LIBRARY(dagtally SHARED
    ${DAG_TALLY_SRCDIR}/Tally.cpp
    ${DAG_TALLY_SRCDIR}/TallyData.cpp
    ${DAG_TALLY_SRCDIR}/TallyManager.cpp
    ${DAG_TALLY_SRCDIR}/MeshTally.cpp
    ${DAG_TALLY_SRCDIR}/TrackLengthMeshTally.cpp
    ${DAG_TALLY_SRCDIR}/KDEMeshTally.cpp
    ${DAG_TALLY_SRCDIR}/CellTally.cpp
    ${DAG_TALLY_SRCDIR}/KDEKernel.cpp
    ${DAG_TALLY_SRCDIR}/KDENeighborhood.cpp
    ${DAG_TALLY_SRCDIR}/PolynomialKernel.cpp
    ${DAG_TALLY_SRCDIR}/Quadrature.cpp
)

TARGET_LINK_LIBRARIES(dagtally lapack)
