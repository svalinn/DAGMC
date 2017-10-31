if (BUILD_STATIC_EXE)
  set(LINK_LIBS_TEST uwuw-static)
else ()
  set(LINK_LIBS_TEST uwuw)
endif ()

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/pyne)
include_directories(${CMAKE_SOURCE_DIR}/uwuw)

set(DRIVERS uwuw_unit_test_driver.cpp)
dagmc_setup_test(uwuw_unit_tests              cpp)
dagmc_setup_test(uwuw_unit_tests_tally        cpp)
dagmc_setup_test(uwuw_unit_tests_preprocessor cpp)

install(FILES test_uwuw.h5m DESTINATION tests)
install(FILES mat_lib.h5    DESTINATION tests)
install(FILES dag_file.h5m  DESTINATION tests)
