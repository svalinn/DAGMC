set(DRIVERS tally_unit_test_driver.cc)

set(LINK_LIBS dagtally)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${GTEST_INCLUDE_DIR})

dagmc_install_test(test_KDEKernel            cpp)
dagmc_install_test(test_KDEMeshTally         cpp)
dagmc_install_test(test_KDENeighborhood      cpp)
dagmc_install_test(test_PolynomialKernel     cpp)
dagmc_install_test(test_Quadrature           cpp)
dagmc_install_test(test_CellTally            cpp)
dagmc_install_test(test_TallyEvent           cpp)
dagmc_install_test(test_TallyData            cpp)
dagmc_install_test(test_Tally                cpp)
dagmc_install_test(test_TrackLengthMeshTally cpp)

install(FILES hashtag_mesh.h5m      DESTINATION tests)
install(FILES rune_mesh.h5m         DESTINATION tests)
install(FILES structured_mesh.h5m   DESTINATION tests)
install(FILES unstr_mesh_split.h5m  DESTINATION tests)
install(FILES unstructured_mesh.h5m DESTINATION tests)
