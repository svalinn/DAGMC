set(DRIVERS dagsolid_magfield_unit_test_driver.cpp 
 ${CMAKE_SOURCE_DIR}/src/geant4/app/src/DagSolidMagneticField.cc)

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/src/geant4)

set(LINK_LIBS_EXTERN_NAMES )

dagmc_install_test(dagsolid_magfield_unit_tests cpp)
dagmc_install_test_file("test_field.txt")
