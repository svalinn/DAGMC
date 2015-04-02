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

set(uwuw_test_tally_source
  ${CMAKE_CURRENT_SOURCE_DIR}/uwuw_test_tally.cpp
)
# geant libs needed for test cases


# two pthreads are needed for final linking, dont know why
# set libraries needed to link to test cases
SET(LIBS
  pthread
  gtest
  dagmciface
  pyne_dagmc
  ${DAGMC_LIBRARIES}
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

# if we want a static executable
IF(STATIC_LIB)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static")
ENDIF(STATIC_LIB)

ADD_EXECUTABLE (uwuw_unit_tests  
    	uwuw_unit_test_driver.cpp
	${uwuw_test_source})

ADD_EXECUTABLE (uwuw_unit_tests_tally  
    	uwuw_unit_test_driver.cpp
	${uwuw_test_tally_source})

INSTALL(TARGETS uwuw_unit_tests DESTINATION tests)
INSTALL(TARGETS uwuw_unit_tests_tally DESTINATION tests)
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/test_uwuw.h5m" DESTINATION tests)
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/mat_lib.h5" DESTINATION tests)

TARGET_LINK_LIBRARIES (uwuw_unit_tests ${LIBS} )
TARGET_LINK_LIBRARIES (uwuw_unit_tests_tally ${LIBS} )

#
RUN_TEST(uwuw_unit_tests)
