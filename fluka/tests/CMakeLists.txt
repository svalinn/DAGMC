#PROJECT(FluDAGUnitTests)
#CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)

# set libraries needed to link to test cases
SET(LIBS
    pthread
    dagmciface
    ${MOAB_DIR}/libMOAB.so
    ${MOAB_DIR}/libdagmc.so
    ${MOAB_DIR}/libMOAB.so
    gtest
    pthread
)

SET (fludag_test_source ${CMAKE_CURRENT_SOURCE_DIR}/test_FlukaFuncs.cpp)

#add_custom_target(check)
macro(run_test test_target)
  add_custom_target(${test_target}_runtest
      COMMAND ${test_target} #cmake 2.6 required
      DEPENDS ${test_target}
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  add_dependencies(${test_target}_runtest gtest)
endmacro()

ADD_EXECUTABLE (fludag_tests 
	EXCLUDE_FROM_ALL
        fludag_unit_test_driver.cc
	${fludag_test_source}
    )

TARGET_LINK_LIBRARIES (fludag_tests ${LIBS})

RUN_TEST(fludag_tests)

	
