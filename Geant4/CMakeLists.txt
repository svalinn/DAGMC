######################################################################
#
# CMakeLists.txt for dagmc Geant4
#
# $Id: CMakeLists.txt 6346 2014-01-27 03:33:21Z cary $
#
######################################################################

# Add in directories in layering order, lowest to highest
set(SCI_INDENT "${SCI_INDENT}  ")
message(STATUS "${SCI_INDENT}Adding dagsolid subdir.")

find_package(SciGeant4 REQUIRED)

LIST(GET Geant4_LIBRARY_DIRS 0 Geant4_LIB_DIR)
FILE(GLOB Geant4_CONFIG "${Geant4_LIB_DIR}/Geant4*")

add_subdirectory(dagsolid)

include(ExternalProject)
ExternalProject_Add(radioprotection
  URL ${CMAKE_CURRENT_SOURCE_DIR}/radioprotection
  CMAKE_ARGS -DGeant4_DIR:PATH='${Geant4_CONFIG}'
     -DCMAKE_INSTALL_PREFIX:PATH='${CMAKE_INSTALL_PREFIX}'
  )
#add_subdirectory(radioprotection)
