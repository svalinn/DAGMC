message(STATUS "Building DAGMC library")

set(SRC_FILES
  DagMC.cpp
  dagmcmetadata.cpp
  ${CMAKE_SOURCE_DIR}/uwuw/uwuw.cpp
  ${CMAKE_SOURCE_DIR}/uwuw/uwuw_preprocessor.cpp
)

set(PUB_HEADERS
  DagMC.hpp
  dagmcmetadata.hpp
  ${CMAKE_SOURCE_DIR}/uwuw/uwuw.hpp
  ${CMAKE_SOURCE_DIR}/uwuw/uwuw_preprocessor.hpp
)

set(LINK_LIBS_SHARED
  ${MOAB_LIBRARIES_SHARED}
  pyne_dagmc
)
set(LINK_LIBS_STATIC
  ${MOAB_LIBRARIES_STATIC}
  lapack blas gfortran quadmath  # Needed for static exes for some reason
  pyne_dagmc-static
)

include_directories(${CMAKE_CURRENT_LIST_DIR})

dagmc_install_library(dagmc)

add_subdirectory(tools)
add_subdirectory(tests)
