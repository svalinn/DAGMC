CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
MESSAGE ( STATUS "Including Common functions in dagmciface library.")

SET ( CODE_SRCDIR      "common" )
SET ( CODE_SRC_FILES   "${CODE_SRCDIR}/dagmciface.cpp"  )
SET ( CODE_PUB_HEADERS "${CODE_SRCDIR}/dagmciface.h" )
SET ( CODE_CONFIG       )

SET ( DI_PUB_HEADERS ${DI_PUB_HEADERS}
                     ${CODE_PUB_HEADERS}
                     ${CODE_CONFIG}
      PARENT_SCOPE)

SET ( DI_SRC_FILES ${DI_SRC_FILES} 
                   ${CODE_SRC_FILES} 
      PARENT_SCOPE)

#ADD_SUBDIRECTORY(tests)

