# CMAKE build script for DAGMC Tally
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
MESSAGE(STATUS "Creating DAGTally Library")
SET ( DAG_TALLY_SRC_FILES
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

LIST( APPEND CODE_PUB_HEADERS "CellTally.hpp" )
LIST( APPEND CODE_PUB_HEADERS "KDEKernel.hpp" )
LIST( APPEND CODE_PUB_HEADERS "KDEMeshTally.hpp" )
LIST( APPEND CODE_PUB_HEADERS "KDENeighborhood.hpp" )
LIST( APPEND CODE_PUB_HEADERS "Matrix3.hpp" )
LIST( APPEND CODE_PUB_HEADERS "MeshTally.hpp" )
LIST( APPEND CODE_PUB_HEADERS "PolynomialKernel.hpp" )
LIST( APPEND CODE_PUB_HEADERS "Quadrature.hpp" )
LIST( APPEND CODE_PUB_HEADERS "TallyData.hpp" )
LIST( APPEND CODE_PUB_HEADERS "TallyEvent.hpp" )
LIST( APPEND CODE_PUB_HEADERS "Tally.hpp" )
LIST( APPEND CODE_PUB_HEADERS "TallyManager.hpp" )
LIST( APPEND CODE_PUB_HEADERS "TrackLengthMeshTally.hpp" )


INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/tally)
#${PROJECT_SOURCE_DIR}/tally/mcnp5)


# Create DAGMC Tally library
IF ( STATIC_LIB )
  ADD_LIBRARY( dagtally STATIC ${DAG_TALLY_SRC_FILES} ${CODE_SRC_FILES} )
ELSE( STATIC_LIB )
  ADD_LIBRARY( dagtally SHARED ${DAG_TALLY_SRC_FILES} ${CODE_SRC_FILES} )
ENDIF( STATIC_LIB )


SET_TARGET_PROPERTIES( dagtally PROPERTIES PUBLIC_HEADER "${CODE_PUB_HEADERS}" )

INSTALL( TARGETS       dagtally 
  LIBRARY        DESTINATION   "${INSTALL_LIB_DIR}"
  ARCHIVE        DESTINATION   "${INSTALL_LIB_DIR}"
  PUBLIC_HEADER  DESTINATION   "${INSTALL_INCLUDE_DIR}"
  )


TARGET_LINK_LIBRARIES(dagtally lapack ${DAGMC_LIBRARIES})

ADD_SUBDIRECTORY(tests)
