#PROJECT(FluDAGUnitTests)
#CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)

SET(pynelibs
	${HDF5_HOME}/lib/libhdf5.so
	${PYNE_HOME}/lib/libpyne.so
	${PYNE_HOME}/lib/libpyne_material.so
	${PYNE_HOME}/lib/libpyne_nucname.so
	${PYNE_HOME}/lib/libpyne_data.so
)
# set libraries needed to link to test cases
SET(LIBS
    pthread
    dagmciface
    ${pynelibs}
    ${MOAB_DIR}/libMOAB.so
    ${MOAB_DIR}/libdagmc.so
    ${MOAB_DIR}/libMOAB.so
    gtest
    pthread
)

SET (fludag_test_source ${CMAKE_CURRENT_SOURCE_DIR}/test_FlukaFuncs.cpp)

add_custom_target(check)
macro(run_test test_target)
  add_custom_target(${test_target}_runtest
      COMMAND ${test_target} #cmake 2.6 required
      DEPENDS ${test_target}
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  add_dependencies(check ${test_target}_runtest gtest)
endmacro()

ADD_EXECUTABLE (fludag_unit_tests 
	EXCLUDE_FROM_ALL
    fludag_unit_test_driver.cc
	${fludag_test_source}
    )

TARGET_LINK_LIBRARIES (fludag_unit_tests ${LIBS})

run_test(fludag_unit_tests)

	
