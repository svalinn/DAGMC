#PROJECT(FluDAGUnitTests)
# 
SET (tgt "gtest")
ADD_CUSTOM_COMMAND (TARGET ${tgt}
	COMMAND ln "-sf"
		"${CMAKE_CURRENT_SOURCE_DIR}/test_geom.h5m"
		"${CMAKE_CURRENT_BINARY_DIR}/test_geom.h5m"
		COMMENT "...linking test file in ${CMAKE_CURRENT_SOURCE_DIR} to ${CMAKE_CURRENT_BINARY_DIR}/slabs.h5m"
	DEPENDS
	VERBATIM
 	)

# fluka test source
SET (fludag_test_source ${CMAKE_CURRENT_SOURCE_DIR}/test_FlukaFuncs.cpp)

# set libraries needed to link to test cases
SET(LIBS
    pthread
    dagmciface
    ${DAG_LINK_LIBS} 
    gtest
    pthread
)

#includes

INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/fluka )
#MESSAGE ( STATUS "Not building dagmciface physics interface library" ${DI_PUB_HEADERS} )

#add_custom_target(check)
macro(run_test test_target)
  add_custom_target(${test_target}_runtest
      COMMAND ${test_target} #cmake 2.6 required
      DEPENDS ${test_target}
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  add_dependencies(${test_target}_runtest gtest)
endmacro()

ADD_EXECUTABLE (fludag_unit_tests 
        fludag_unit_test_driver.cc
	${fludag_test_source}
    )
# install into tests dir
INSTALL(TARGETS fludag_unit_tests DESTINATION tests)
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/slabs.h5m" DESTINATION tests)

TARGET_LINK_LIBRARIES (fludag_unit_tests 
	${FLUDAG_LIBS} 
        gtest
	${FLUKA_DIR}/libflukahp.a
	gfortran
        pthread )

RUN_TEST(fludag_unit_tests)

	
