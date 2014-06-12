######################################################################
#
# CMakeLists.txt for dagmc fluka
#
######################################################################

MESSAGE ( STATUS "Including FLUKA in dagmciface library.")
SET ( DAG_FLUKA_SRCDIR "fluka" )

SET ( DAGMCIFACE_FLUKA_SRC "${DAG_FLUKA_SRCDIR}/fluka_funcs.cpp" )

SET ( DAGMCIFACE_FLUKA_HEADERS "${DAG_FLUKA_SRCDIR}/fluka_funcs.h")

SET( dagmciface_headers ${dagmciface_headers} ${DAGMCIFACE_FLUKA_HEADERS} PARENT_SCOPE)

SET ( DAGMCIFACE_SRC ${DAGMCIFACE_SRC} ${DAGMCIFACE_FLUKA_SRC} PARENT_SCOPE)
INCLUDE_DIRECTORIES( DAG_FLUKA_SRCDIR )

ADD_SUBDIRECTORY(tests EXCLUDE_FROM_ALL)
