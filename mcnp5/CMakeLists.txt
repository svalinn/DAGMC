CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
MESSAGE ( STATUS "Including MCNP5 in dagmciface library.")
SET ( DAG_MCNP_SRCDIR "mcnp5" )
SET ( DAGMCIFACE_MCNP_SRC "${DAG_MCNP_SRCDIR}/mcnp_funcs.cpp"  )
SET ( DAGMCIFACE_MCNP_HEADERS "${DAG_MCNP_SRCDIR}/mcnp_funcs.h" )
SET ( DAGMCIFACE_MCNP_CONFIG "${DAG_MCNP_SRCDIR}/dagmcnp5.make" )

SET( dagmciface_headers ${dagmciface_headers} 
                        ${DAGMCIFACE_MCNP_HEADERS} 
                        ${DAGMCIFACE_MCNP_CONFIG}
   )

SET ( DAGMCIFACE_SRC ${DAGMCIFACE_SRC} ${DAGMCIFACE_MCNP_SRC} )
INCLUDE_DIRECTORIES( "${DAG_MCNP_SRCDIR}")
