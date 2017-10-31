set(DRIVERS uwuw_unit_test_driver.cpp)

set(LINK_LIBS uwuw)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/pyne)
include_directories(${CMAKE_SOURCE_DIR}/uwuw)

dagmc_install_test(uwuw_unit_tests              cpp)
dagmc_install_test(uwuw_unit_tests_tally        cpp)
dagmc_install_test(uwuw_unit_tests_preprocessor cpp)

install(FILES test_uwuw.h5m DESTINATION tests)
install(FILES mat_lib.h5    DESTINATION tests)
install(FILES dag_file.h5m  DESTINATION tests)
