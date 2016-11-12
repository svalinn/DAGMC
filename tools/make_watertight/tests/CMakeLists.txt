message(STATUS "Building make_watertight tests")

# Source files
set(SRC_FILES)
list(APPEND SRC_FILES ../Gen.cpp)
list(APPEND SRC_FILES ../Arc.cpp)
list(APPEND SRC_FILES ../Zip.cpp)
list(APPEND SRC_FILES ../Cleanup.cpp)
list(APPEND SRC_FILES ../MakeWatertight.cpp)
list(APPEND SRC_FILES ../CheckWatertight.cpp)
list(APPEND SRC_FILES test_classes.cpp)

# Public headers
set(PUB_HEADERS)
list(APPEND PUB_HEADERS ../Gen.hpp)
list(APPEND PUB_HEADERS ../Arc.hpp)
list(APPEND PUB_HEADERS ../Zip.hpp)
list(APPEND PUB_HEADERS ../Cleanup.hpp)
list(APPEND PUB_HEADERS ../CheckWatertight.hpp)
list(APPEND PUB_HEADERS ../MakeWatertight.hpp)
list(APPEND PUB_HEADERS test_classes.hpp)

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS gtest)

# Define executables
add_executable(make_watertight_cylinder_tests mw_unit_test_driver.cc test_cyl.cpp ${SRC_FILES})
add_executable(make_watertight_cone_tests mw_unit_test_driver.cc test_cones.cpp ${SRC_FILES})
add_executable(make_watertight_no_curve_sphere_tests mw_unit_test_driver.cc test_no_curve_sphere.cpp ${SRC_FILES})
add_executable(make_watertight_regression_tests mw_unit_test_driver.cc regression_tests.cpp ${SRC_FILES})

# Includes
include_directories(..)
include_directories(${GTEST_INCLUDE_DIR})
include_directories(${MOAB_INCLUDE_DIRS})

# Link to libraries
target_link_libraries(make_watertight_cylinder_tests ${LINK_LIBS})
target_link_libraries(make_watertight_cone_tests ${LINK_LIBS})
target_link_libraries(make_watertight_no_curve_sphere_tests ${LINK_LIBS})
target_link_libraries(make_watertight_regression_tests ${LINK_LIBS})

# Install tests and test models
install(TARGETS make_watertight_cylinder_tests make_watertight_cone_tests make_watertight_no_curve_sphere_tests DESTINATION tests PERMISSIONS ${PERMS})
install(TARGETS make_watertight_regression_tests DESTINATION tests PERMISSIONS ${PERMS})
install(FILES cones.h5m cyl.h5m no_curve_sphere.h5m DESTINATION tests)


# Add to unit test framework
add_test(make_watertight_cylinder_tests make_watertight_cylinder_tests make_watertight_no_curve_sphere_tests)
add_test(make_watertight_cone_tests make_watertight_cone_tests)

# enable tests
enable_testing()
