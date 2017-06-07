message(STATUS "Building DAGMC unit tests")

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS pthread)
list(APPEND LINK_LIBS gtest)
if (BUILD_STATIC)
  list(APPEND LINK_LIBS dagmc-static)
  list(APPEND LINK_LIBS pyne_dagmc-static)
else ()
  list(APPEND LINK_LIBS dagmc)
  list(APPEND LINK_LIBS pyne_dagmc)
endif (BUILD_STATIC)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Include directories
include_directories(.)
include_directories(${MOAB_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/dagmc)
include_directories(${GTEST_INCLUDE_DIR})

# Setup dagmc unit tests
set(DRIVER dagmc_unit_test_driver.cc)
setup_test(dagmc_unit_tests      cpp ${DRIVER} "${LINK_LIBS}")
setup_test(dagmc_pointinvol_test cpp ${DRIVER} "${LINK_LIBS}")
setup_test(dagmc_rayfire_test    cpp ${DRIVER} "${LINK_LIBS}")
setup_test(dagmc_simple_test     cpp ${DRIVER} "${LINK_LIBS}")
#setup_test(pt_vol_test           cc  ${DRIVER} "${LINK_LIBS}")
#setup_test(ray_fire_test         cc  ${DRIVER} "${LINK_LIBS}")
#setup_test(test_geom             cc  ${DRIVER} "${LINK_LIBS}")

# Install h5m files
install(FILES test_dagmc.h5m      DESTINATION tests)
install(FILES test_dagmc_impl.h5m DESTINATION tests)
install(FILES test_geom.h5m       DESTINATION tests)

# Enable unit tests
enable_testing()
