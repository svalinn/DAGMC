######################################################################
#
# CMakeLists.txt for Geant4 of dagmc
#
# $Id: CMakeLists.txt 6426 2014-02-21 17:56:17Z cary $
#
# Copyright &copy; 2002-2012, University of Colorado and Tech-X Corporation
# See LICENSE file for conditions of use.
#
######################################################################

# Project information
project(dagmcg4)
set(VERSION_MAJOR "0")
set(VERSION_MINOR "1")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

# Required version
cmake_minimum_required(VERSION 2.8)

######################################################################
#
# Grab the scimake repo
#
######################################################################

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/scimake)
  message(STATUS "${CMAKE_CURRENT_SOURCE_DIR}/scimake found.  Update needed?")
else ()
  execute_process(COMMAND svn co https://svn.code.sf.net/p/scimake/code/trunk scimake
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endif ()

######################################################################
#
# Usual initialization stuff, sciInit includes txfindpackage and
# takes care of MPI
#
######################################################################

if (FALSE)
set(NOFORTRAN TRUE)
set(NO_CONFIG_H TRUE)
include(${PROJECT_SOURCE_DIR}/scimake/SciInit.cmake)
include_directories(${SCIMAKE_DIR}/include)
set(CTK_VERSION "${VERSION}-${PROJECT_REV}")
install(CODE "message(STATUS \"Installing into \${CMAKE_INSTALL_PREFIX}.\")")
include(${PROJECT_SOURCE_DIR}/txcmake/TxInit.cmake)

if(CMAKE_COMPILER_IS_GNUCXX AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "4.4")
  # Only works on gcc > 4.4
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=return-type")
endif()

######################################################################
#
# Always use rpath to greatest extent.
# See: http://www.itk.org/Wiki/CMake_RPATH_handling
# Add -DCMAKE_SKIP_RPATH:BOOL=TRUE to prevent any rpath handling
#
######################################################################

# SET(CMAKE_SKIP_BUILD_RPATH FALSE)
# Use, i.e. don't skip the full RPATH for the build tree
# Not needed since we build static libs

# SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
# When building, don't link with the install RPATH at build,
# but add later on when installing.
# Not needed since we build static libs

# SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# The RPATH to be used when installing
# Not needed since we build static libs

if (NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_PATH)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()
# Add the automatically determined parts of the RPATH that
# point to directories outside the build tree to the install RPATH

######################################################################
#
# Find packages and set includes
# For libraries that are used for ALL builds
#
######################################################################

# To find the right installations
set(USE_CC4PY_LIBS TRUE)

# Find Qt
set(QT_REQUIRED_LIBRARIES
  QtCore QtGui QtOpenGL QtXml QtXmlPatterns QtNetwork QtWebKit
)
find_package(SciQt4 COMPONENTS ${QT_REQUIRED_LIBRARIES} REQUIRED)
include(SciFindQtPkg)
# SciFindQtPkg(PACKAGE Qt3D HEADERS qglview.h LIBRARIES Qt3D)
find_package(SciPythonLibs REQUIRED)

######################################################################
#
# Find visit related packages
#
######################################################################

SciGetInstSubdirs(visit4ctk instDirs)
set(VisIt_INSTALL_DIRS ${instDirs})
find_package(SciVisIt REQUIRED)

set(Vtk_INSTALL_DIRS vtk4ctk-cc4py vtk4ctk-sersh vtk4ctk)
find_package(SciVtk REQUIRED)

# Add VTK headers to the build
include_directories(${Vtk_INCLUDE_DIRS})

###################################################
#
# For remote operations
#
###################################################

#if (ENABLE_REMOTE)
# Make sure we always find the shared installations.
find_package(TxSsh  REQUIRED)
find_package(SciLibssh REQUIRED)
#endif ()

find_package(TxGml)
if (TXGML_FOUND)
  add_definitions(${TxGml_DEFINITIONS})
  find_package(SciOce COMPONENTS XdeIges XdeStep Xde Mesh Step Stl REQUIRED)
  include(SciFindQtPkg)
  SciFindQtPkg(PACKAGE Qt3D HEADERS qglview.h LIBRARIES Qt3D)
  find_package(TxUtils REQUIRED)
  include_directories(
    ${TxGml_INCLUDE_DIRS}
    ${Oce_INCLUDE_DIRS}
    ${QT_QT3D_INCLUDE_DIRS}
    ${TxUtils_INCLUDE_DIRS}
  )
endif ()

# Licensing
include(${PROJECT_SOURCE_DIR}/txcmake/TxEnableSecurity.cmake)

###################################################
#
# Boost and compression
#
###################################################

# Boost needed on windows, zlib is not
find_package(SciBoost REQUIRED thread chrono system)
find_package(SciSz)
find_package(SciZ REQUIRED)

###################################################
#
# For basic capabilities
#
###################################################

find_package(SciTxBase REQUIRED)
find_package(SciHdf5 REQUIRED)

###################################################
#
# Code subdirectories that depend on finding something.
#
###################################################

message(STATUS "Adding src subdirectory.")
add_subdirectory(src)
message(STATUS "Adding ctkqt subdirectory.")
add_subdirectory(ctkqt)
message(STATUS "Adding resources subdirectory.")
add_subdirectory(resources)
message(STATUS "Adding ctkscripts subdirectory.")
add_subdirectory(ctkscripts)

###################################################
#
# Documentation
#
###################################################

find_package(SciDoxygen)
if (${DOXYGEN_FOUND})
  message(STATUS "Adding ctkdevdocs subdirectory.")
  add_subdirectory(ctkdevdocs)
else ()
  message(WARNING "Doxygen not found. Not adding ctkdevdocs sub directory")
endif ()

###################################################
#
# Testing
#
###################################################

# Enable tests
include(CTestConfig.cmake)
include(CTest)
enable_testing()

######################################################################
#
# Configure files
#
######################################################################

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CtkConfig.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/CtkConfig.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CtkVersion.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/CtkVersion.h)

######################################################################
#
# Installations at this level
#
######################################################################

if (EXISTS ${PROJECT_SOURCE_DIR}/svninfo.txt)
  install(FILES svninfo.txt DESTINATION share)
  message(STATUS "--- svninfo.txt found. Will be installed")
else ()
  message(STATUS "--- svninfo.txt NOT found. Skipping install")
endif ()

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt
  DESTINATION share
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/CtkVersion.h
  ${CMAKE_CURRENT_BINARY_DIR}/CtkConfig.h
  ${CMAKE_CURRENT_SOURCE_DIR}/CtkConfig.h.in  # For configuration by composers
  DESTINATION include
)

# Install the entire resources directory
install(DIRECTORY resources DESTINATION . PATTERN ".svn" EXCLUDE)

######################################################################
#
# Packaging
#
######################################################################

# CPack version numbers for release tarball name.
set(CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH}-r${PROJECT_REV})

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Tech-X ComposerToolkit Library")
set(CPACK_PACKAGE_VENDOR "Tech-X Corporation")
set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README.txt)
set(CPACK_SOURCE_PACKAGE_FILE_NAME "ctk-${VERSION}-r${PROJECT_REV}"
    CACHE INTERNAL "tarball basename")
set(CPACK_SOURCE_GENERATOR TGZ)
set(CPACK_SOURCE_IGNORE_FILES
"/CVS/;/.svn/;.swp$;.#;/#;/build/;/serial/;/parallel/;~$;/autom4te.cache/;/.config;preconfig.out")
set(CPACK_GENERATOR TGZ)
include(CPack)

endif ()

