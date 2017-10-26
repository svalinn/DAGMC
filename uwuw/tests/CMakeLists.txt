message(STATUS "Building UWUW unit tests")

set(LINK_LIBS_SHARED gtest dagmc)
set(LINK_LIBS_STATIC gtest dagmc-static)

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_EXE ${LINK_LIBS_STATIC})
else ()
  set(LINK_LIBS_EXE ${LINK_LIBS_SHARED})
endif ()

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/pyne)
include_directories(${CMAKE_SOURCE_DIR}/uwuw)

set(DRIVER uwuw_unit_test_driver.cpp)
dagmc_setup_test(uwuw_unit_tests              cpp ${DRIVER} "${LINK_LIBS_EXE}")
dagmc_setup_test(uwuw_unit_tests_tally        cpp ${DRIVER} "${LINK_LIBS_EXE}")
dagmc_setup_test(uwuw_unit_tests_preprocessor cpp ${DRIVER} "${LINK_LIBS_EXE}")

install(FILES test_uwuw.h5m DESTINATION tests)
install(FILES mat_lib.h5    DESTINATION tests)
install(FILES dag_file.h5m  DESTINATION tests)
