set(LINK_LIBS dagmc)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${CMAKE_SOURCE_DIR}/dagmc)

set(SRC_FILES pt_vol_test.cpp)
dagmc_install_exe(pt_vol_test)
set(SRC_FILES ray_fire_test.cpp)
dagmc_install_exe(ray_fire_test)
set(SRC_FILES test_geom.cpp)
dagmc_install_exe(test_geom)
