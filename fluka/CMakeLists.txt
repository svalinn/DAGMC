CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
MESSAGE ( STATUS "Including FLUKA in dagmciface library.")

SET ( CODE_SRCDIR      "fluka" )
SET ( CODE_SRC_FILES   "${CODE_SRCDIR}/fluka_funcs.cpp"  )
SET ( CODE_PUB_HEADERS "${CODE_SRCDIR}/fluka_funcs.h" )
SET ( CODE_CONFIG       )

SET ( DI_PUB_HEADERS ${DI_PUB_HEADERS}
                     ${CODE_PUB_HEADERS}
                     ${CODE_CONFIG}
      PARENT_SCOPE)

SET ( DI_SRC_FILES ${DI_SRC_FILES} 
                   ${CODE_SRC_FILES} 
      PARENT_SCOPE)

INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/${CODE_SRCDIR} )

ADD_SUBDIRECTORY(tests)


