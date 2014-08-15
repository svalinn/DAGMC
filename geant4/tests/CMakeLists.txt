PROJECT(DagSolidUnitTests)

SET (tgt "gtest")
ADD_CUSTOM_COMMAND (TARGET ${tgt}
	COMMAND ln "-sf"
		"${CMAKE_CURRENT_SOURCE_DIR}/test_geom.h5m"
		"${CMAKE_CURRENT_BINARY_DIR}/test_geom.h5m"
		COMMENT "...linking test file in ${CMAKE_CURRENT_SOURCE_DIR} to ${CMAKE_CURRENT_BINARY_DIR}/slabs.h5m"
	DEPENDS
	VERBATIM
 	)

# dagsolid tests
set(dagsolid_test_source
  ${CMAKE_CURRENT_SOURCE_DIR}/DagSolid_test.cpp
)

# geant libs needed for test cases
set(GEANT_LIBS
  ${GEANT4_DIR}/lib/libG4geometry.so
  ${GEANT4_DIR}/lib/libG4clhep.so
  ${GEANT4_DIR}/lib/libG4global.so
)


# two pthreads are needed for final linking, dont know why
# set libraries needed to link to test cases
SET(LIBS
  pthread
  ${GEANT_LIBS}
  ${MOAB_LIBRARIES}
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

# custom driver which calls tests
#ADD_CUSTOM_TARGET (dagsolid_tests  DEPENDS gtest
#    	dagsolid_unit_test_driver.cpp
#	${dagsolid_test_source})

ADD_EXECUTABLE (dagsolid_tests  
  EXCLUDE_FROM_ALL
    	dagsolid_unit_test_driver.cpp
	${dagsolid_test_source})
# INSTALL(TARGETS dagsolid_tests DESTINATION tests)

TARGET_LINK_LIBRARIES (dagsolid_tests ${LIBS} dagsolid )

# enable test cases
ENABLE_TESTING()

# list test cases for compatibility with BaTLab
ADD_TEST(DagSolid_test DagSolid_test)

#
RUN_TEST(dagsolid_tests)