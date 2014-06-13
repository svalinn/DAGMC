CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
MESSAGE ( STATUS "Including MCNP5 in dagmciface library.")

SET ( CODE_SRCDIR      "mcnp5" )
SET ( CODE_SRC_FILES   "${CODE_SRCDIR}/mcnp_funcs.cpp"  )
SET ( CODE_PUB_HEADERS "${CODE_SRCDIR}/mcnp_funcs.h" )
SET ( CODE_CONFIG      "${CODE_SRCDIR}/dagmcnp5.make" )

SET ( DI_PUB_HEADERS ${DI_PUB_HEADERS}
                     ${CODE_PUB_HEADERS}
                     ${CODE_CONFIG}
      PARENT_SCOPE)

SET ( DI_SRC_FILES ${DI_SRC_FILES} 
                   ${CODE_SRC_FILES} 
      PARENT_SCOPE)

#ADD_SUBDIRECTORY(tests)
