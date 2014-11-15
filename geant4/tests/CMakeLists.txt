# PROJECT(DagSolidUnitTests)

# dagsolid tests
set(dagsolid_test_source
  ${CMAKE_CURRENT_SOURCE_DIR}/DagSolid_test.cpp
)

# geant libs needed for test cases


# two pthreads are needed for final linking, dont know why
# set libraries needed to link to test cases
SET(LIBS
  pthread
  ${Geant4_LIBRARIES}
  ${DAGMC_LIBRARIES}
  gtest
  pthread
)

macro(run_test test_target)
  add_custom_target(${test_target}_runtest
      COMMAND ${test_target} #cmake 2.6 required
      DEPENDS ${test_target}
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  add_dependencies(${test_target}_runtest gtest)
endmacro()

# for the DagSolid includes
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/../")

ADD_EXECUTABLE (dagsolid_unit_tests  
    	dagsolid_unit_test_driver.cpp
	${dagsolid_test_source})

INSTALL(TARGETS dagsolid_unit_tests DESTINATION tests)
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/test_geom.h5m" DESTINATION tests)

TARGET_LINK_LIBRARIES (dagsolid_unit_tests ${LIBS} dagsolid )

#
RUN_TEST(dagsolid_unit_tests)
