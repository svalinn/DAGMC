# CMAKE build script for DAGMC Tally unit testing

# CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Where to look for includes
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/tally )

MESSAGE(STATUS "${GTEST_HOME}")

# if we want a static executable
IF(STATIC_LIB)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static")
ENDIF(STATIC_LIB)

# set libraries
SET(LIBRARIES
    pthread
    ${DAGMC_LIBRARIES}
    dagtally
    gtest
    pthread
)

#add_custom_target(check)
macro(run_test test_target)
  add_custom_target(${test_target}_runtest
      COMMAND ${test_target} #cmake 2.6 required
      DEPENDS ${test_target}
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  add_dependencies(${test_target}_runtest gtest)
endmacro()

# set up DAGMC Tally test cases
ADD_EXECUTABLE(test_KDEKernel test_KDEKernel.cpp tally_unit_test_driver.cc)
TARGET_LINK_LIBRARIES(test_KDEKernel ${LIBRARIES})

ADD_EXECUTABLE(test_KDEMeshTally test_KDEMeshTally.cpp tally_unit_test_driver.cc)
TARGET_LINK_LIBRARIES(test_KDEMeshTally ${LIBRARIES})

ADD_EXECUTABLE(test_KDENeighborhood test_KDENeighborhood.cpp tally_unit_test_driver.cc)
TARGET_LINK_LIBRARIES(test_KDENeighborhood ${LIBRARIES})

ADD_EXECUTABLE(test_PolynomialKernel test_PolynomialKernel.cpp tally_unit_test_driver.cc)
TARGET_LINK_LIBRARIES(test_PolynomialKernel ${LIBRARIES})

ADD_EXECUTABLE(test_Quadrature test_Quadrature.cpp tally_unit_test_driver.cc)
TARGET_LINK_LIBRARIES(test_Quadrature ${LIBRARIES})

ADD_EXECUTABLE(test_CellTally test_CellTally.cpp tally_unit_test_driver.cc)
TARGET_LINK_LIBRARIES(test_CellTally ${LIBRARIES})

ADD_EXECUTABLE(test_TallyEvent test_TallyEvent.cpp tally_unit_test_driver.cc)
TARGET_LINK_LIBRARIES(test_TallyEvent ${LIBRARIES})

ADD_EXECUTABLE(test_TallyData test_TallyData.cpp tally_unit_test_driver.cc)
TARGET_LINK_LIBRARIES(test_TallyData ${LIBRARIES})

ADD_EXECUTABLE(test_Tally test_Tally.cpp tally_unit_test_driver.cc)
TARGET_LINK_LIBRARIES(test_Tally ${LIBRARIES})

ADD_EXECUTABLE(test_TrackLengthMeshTally test_TrackLengthMeshTally.cpp tally_unit_test_driver.cc)
TARGET_LINK_LIBRARIES(test_TrackLengthMeshTally ${LIBRARIES})


INSTALL(TARGETS test_KDEMeshTally DESTINATION tests)
INSTALL(TARGETS test_KDEKernel DESTINATION tests)
INSTALL(TARGETS test_KDENeighborhood DESTINATION tests)
INSTALL(TARGETS test_PolynomialKernel DESTINATION tests)
INSTALL(TARGETS test_Quadrature DESTINATION tests)
INSTALL(TARGETS test_CellTally DESTINATION tests)
INSTALL(TARGETS test_TallyEvent DESTINATION tests)
INSTALL(TARGETS test_TallyData DESTINATION tests)
INSTALL(TARGETS test_Tally DESTINATION tests)
INSTALL(TARGETS test_TrackLengthMeshTally DESTINATION tests)

INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/structured_mesh.h5m" DESTINATION tests)
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/unstructured_mesh.h5m" DESTINATION tests)
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/unstr_mesh_split.h5m" DESTINATION tests)
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/rune_mesh.h5m" DESTINATION tests)



# enable DAGMC Tally test cases
ENABLE_TESTING()

# list DAGMC Tally test cases
ADD_TEST(test_KDEKernel test_KDEKernel)
ADD_TEST(test_KDEMeshTally test_KDEMeshTally)
ADD_TEST(test_KDENeighborhood test_KDENeighborhood)
ADD_TEST(test_PolynomialKernel test_PolynomialKernel)
ADD_TEST(test_Quadrature test_Quadrature)
ADD_TEST(test_CellTally test_CellTally)
ADD_TEST(test_TallyEvent test_TallyEvent)
ADD_TEST(test_TallyData test_TallyData)
ADD_TEST(test_Tally test_Tally)
ADD_TEST(test_TrackLengthMeshTally test_TrackLengthMeshTally)
