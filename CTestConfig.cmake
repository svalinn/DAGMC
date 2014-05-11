######################################################################
#
# CTestConfig: Set variables for testing.
#
# $Id: CTestConfig.cmake 22535 2014-02-11 23:21:43Z veitzer $
#
# Copyright 2010-2013 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
######################################################################

#
# enable_testing() and include(CTest) must be called after this is included.
#
# This is actually set by Bilder
if (CTEST_DROP_SITE)
  include(${PROJECT_SOURCE_DIR}/scimake/SciTestConfig.cmake)
endif ()

