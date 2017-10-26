message(STATUS "Building DAGMC unit tests")

set(LINK_LIBS_SHARED gtest dagmc)
set(LINK_LIBS_STATIC gtest dagmc-static)

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_EXE ${LINK_LIBS_STATIC})
else ()
  set(LINK_LIBS_EXE ${LINK_LIBS_SHARED})
endif ()

include_directories(${GTEST_INCLUDE_DIR})

set(DRIVER dagmc_unit_test_driver.cc)
dagmc_setup_test(dagmc_unit_tests      cpp ${DRIVER} "${LINK_LIBS_EXE}")
dagmc_setup_test(dagmc_pointinvol_test cpp ${DRIVER} "${LINK_LIBS_EXE}")
dagmc_setup_test(dagmc_rayfire_test    cpp ${DRIVER} "${LINK_LIBS_EXE}")
dagmc_setup_test(dagmc_simple_test     cpp ${DRIVER} "${LINK_LIBS_EXE}")

install(FILES test_dagmc.h5m      DESTINATION tests)
install(FILES test_dagmc_impl.h5m DESTINATION tests)
install(FILES test_geom.h5m       DESTINATION tests)
