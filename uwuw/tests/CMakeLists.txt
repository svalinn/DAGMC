# PROJECT(UWUWUnitTests)

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
set(uwuw_test_source
  ${CMAKE_CURRENT_SOURCE_DIR}/uwuw_test.cpp
)

# geant libs needed for test cases


# two pthreads are needed for final linking, dont know why
# set libraries needed to link to test cases
SET(LIBS
  pthread
  gtest
  ${PYNE_LIBS_DIR}/libpyne.so
  ${HDF5_LIBRARIES}
  dagmciface
  ${MOAB_LIBRARIES}
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

ADD_EXECUTABLE (uwuw_unit_tests  
    	uwuw_unit_test_driver.cpp
	${uwuw_test_source})

INSTALL(TARGETS uwuw_unit_tests DESTINATION tests)
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/test_geom.h5m" DESTINATION tests)

TARGET_LINK_LIBRARIES (uwuw_unit_tests ${LIBS} )

#
RUN_TEST(uwuw_unit_tests)