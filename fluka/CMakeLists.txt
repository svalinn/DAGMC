CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
MESSAGE ( STATUS "Including FLUKA in dagmciface library.")

SET ( CODE_SRCDIR      "fluka" )
SET ( CODE_SRC_FILES   "${CODE_SRCDIR}/fluka_funcs.cpp" 
	               "${CODE_SRCDIR}/pyne.cpp")
SET ( CODE_PUB_HEADERS "${CODE_SRCDIR}/fluka_funcs.h" 
	               "${CODE_SRCDIR}/pyne.h")
SET ( CODE_CONFIG       )

SET ( DI_PUB_HEADERS ${DI_PUB_HEADERS}
                     ${CODE_PUB_HEADERS}
                     ${CODE_CONFIG}
      PARENT_SCOPE)

SET ( DI_SRC_FILES ${DI_SRC_FILES} 
                   ${CODE_SRC_FILES} 
      PARENT_SCOPE)

INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/${CODE_SRCDIR} )

ENABLE_LANGUAGE( Fortran )

SET(FLUDAG_LIBS
    dagmciface
    ${HDF5_HOME}/lib/libhdf5.so
    # ${PYNE_LIBS_DIR}/libpyne.so
    ${MOAB_DIR}/libMOAB.so
    ${MOAB_DIR}/libdagmc.so
    ${MOAB_DIR}/libMOAB.so
    gfortran
)
# If making mainfludag, define FLUPRO)
# add_custom_command (TARGET mainfludag
# PRE_BUILD COMMAND 
IF (NOT DEFINED FLUPRO) 
   set (FLUPRO "$ENV{FLUPRO}")
ENDIF (NOT DEFINED FLUPRO)
	   # 	COMMENT "Setting $FLUPRO"
	   # )
ADD_EXECUTABLE ( mainfludag EXCLUDE_FROM_ALL mainFluDAG.cpp )
TARGET_LINK_LIBRARIES ( mainfludag 
	${FLUDAG_LIBS}  
	${FLUPRO}/libflukahp.a )

ADD_SUBDIRECTORY(tests)

