# CMAKE build script for DAGMC Tally
SET ( DAG_TALLY_SRCDIR "tally" )

SET ( DAG_TALLY_FILES 
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

INCLUDE_DIRECTORIES( ${DAG_TALLY_SRCDIR} )

IF ( BUILD_MCNP5 OR BUILD_ALL)
  SET ( DAG_TALLY_FILES ${DAG_TALLY_FILES}
        ${DAG_TALLY_SRCDIR}/mcnp5/meshtal_funcs.cpp )
  INCLUDE_DIRECTORIES(${DAG_TALLY_SRCDIR}/mcnp5)
ENDIF ()


# Create DAGMC Tally library
ADD_LIBRARY(dagtally SHARED ${DAG_TALLY_FILES} )

TARGET_LINK_LIBRARIES(dagtally lapack)
