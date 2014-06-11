# CMAKE build script for DAGMC Tally

# Create DAGMC Tally library
ADD_LIBRARY(dagtally SHARED
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

TARGET_LINK_LIBRARIES(dagtally lapack)

# Add unit testing subdirectory
# TODO: omitting test files for now
#ADD_SUBDIRECTORY(test)
