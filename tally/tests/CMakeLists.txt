message(STATUS "Building dagtally unit tests")

# Includes
include_directories(..)

# Libraries
set(TALLY_TEST_LIBS)
list(APPEND TALLY_TEST_LIBS gtest)
list(APPEND TALLY_TEST_LIBS dagtally)
list(APPEND TALLY_TEST_LIBS pthread)
list(APPEND TALLY_TEST_LIBS ${DAG_LINK_LIBS})

# Setup unit tests
set(DRIVER tally_unit_test_driver.cc)
setup_test(test_KDEKernel            cpp ${DRIVER} "${TALLY_TEST_LIBS}")
setup_test(test_KDEMeshTally         cpp ${DRIVER} "${TALLY_TEST_LIBS}")
setup_test(test_KDENeighborhood      cpp ${DRIVER} "${TALLY_TEST_LIBS}")
setup_test(test_PolynomialKernel     cpp ${DRIVER} "${TALLY_TEST_LIBS}")
setup_test(test_Quadrature           cpp ${DRIVER} "${TALLY_TEST_LIBS}")
setup_test(test_CellTally            cpp ${DRIVER} "${TALLY_TEST_LIBS}")
setup_test(test_TallyEvent           cpp ${DRIVER} "${TALLY_TEST_LIBS}")
setup_test(test_TallyData            cpp ${DRIVER} "${TALLY_TEST_LIBS}")
setup_test(test_Tally                cpp ${DRIVER} "${TALLY_TEST_LIBS}")
setup_test(test_TrackLengthMeshTally cpp ${DRIVER} "${TALLY_TEST_LIBS}")

# Install h5m files
install(FILES hashtag_mesh.h5m      DESTINATION ${INSTALL_SHARE_DIR}/tests)
install(FILES rune_mesh.h5m         DESTINATION ${INSTALL_SHARE_DIR}/tests)
install(FILES structured_mesh.h5m   DESTINATION ${INSTALL_SHARE_DIR}/tests)
install(FILES unstr_mesh_split.h5m  DESTINATION ${INSTALL_SHARE_DIR}/tests)
install(FILES unstructured_mesh.h5m DESTINATION ${INSTALL_SHARE_DIR}/tests)

# enable unit tests
enable_testing()
