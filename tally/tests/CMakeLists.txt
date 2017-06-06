message(STATUS "Building dagtally unit tests")

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS pthread)
list(APPEND LINK_LIBS gtest)
if (BUILD_STATIC)
  list(APPEND LINK_LIBS dagtally-static)
else ()
  list(APPEND LINK_LIBS dagtally)
endif (BUILD_STATIC)

# Include directories
include_directories(..)
include_directories(${GTEST_INCLUDE_DIR})

# CXX compiler flags
if (BUILD_STATIC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--allow-multiple-definition")
endif (BUILD_STATIC)
if (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)

# Setup unit tests
set(DRIVER tally_unit_test_driver.cc)
setup_test(test_KDEKernel            cpp ${DRIVER} "${LINK_LIBS}")
setup_test(test_KDEMeshTally         cpp ${DRIVER} "${LINK_LIBS}")
setup_test(test_KDENeighborhood      cpp ${DRIVER} "${LINK_LIBS}")
setup_test(test_PolynomialKernel     cpp ${DRIVER} "${LINK_LIBS}")
setup_test(test_Quadrature           cpp ${DRIVER} "${LINK_LIBS}")
setup_test(test_CellTally            cpp ${DRIVER} "${LINK_LIBS}")
setup_test(test_TallyEvent           cpp ${DRIVER} "${LINK_LIBS}")
setup_test(test_TallyData            cpp ${DRIVER} "${LINK_LIBS}")
setup_test(test_Tally                cpp ${DRIVER} "${LINK_LIBS}")
setup_test(test_TrackLengthMeshTally cpp ${DRIVER} "${LINK_LIBS}")

# Install h5m files
install(FILES hashtag_mesh.h5m      DESTINATION tests)
install(FILES rune_mesh.h5m         DESTINATION tests)
install(FILES structured_mesh.h5m   DESTINATION tests)
install(FILES unstr_mesh_split.h5m  DESTINATION tests)
install(FILES unstructured_mesh.h5m DESTINATION tests)

# Enable unit tests
enable_testing()
