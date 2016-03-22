message(STATUS "Building dagtally unit tests")

# Includes
include_directories(..)

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS gtest)
list(APPEND LINK_LIBS dagtally)
list(APPEND LINK_LIBS pthread)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

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
install(FILES hashtag_mesh.h5m      DESTINATION ${INSTALL_SHARE_DIR}/tests)
install(FILES rune_mesh.h5m         DESTINATION ${INSTALL_SHARE_DIR}/tests)
install(FILES structured_mesh.h5m   DESTINATION ${INSTALL_SHARE_DIR}/tests)
install(FILES unstr_mesh_split.h5m  DESTINATION ${INSTALL_SHARE_DIR}/tests)
install(FILES unstructured_mesh.h5m DESTINATION ${INSTALL_SHARE_DIR}/tests)

# Enable unit tests
enable_testing()
