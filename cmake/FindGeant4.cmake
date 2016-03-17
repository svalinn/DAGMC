# Try to find Geant4
#
# Once done this will define
#
#  GEANT4_FOUND - System has Geant4
#  GEANT4_INCLUDE_DIR - The Geant4 include dir
#  GEANT4_LIBRARIES - In case they are ever needed

find_path(GEANT4_CMAKE_CONFIG
          NAMES Geant4Config.cmake
          HINTS ${GEANT4_DIR}
                ${GEANT4_DIR}/lib/Geant4-10.0.0
                ${GEANT4_DIR}/lib/Geant4-10.0.1
                ${GEANT4_DIR}/lib/Geant4-10.0.2
                ${GEANT4_DIR}/lib/Geant4-10.0.3
                ${GEANT4_DIR}/lib/Geant4-10.0.4
                ${GEANT4_DIR}/lib/Geant4-10.1.0
                ${GEANT4_DIR}/lib/Geant4-10.1.1
                ${GEANT4_DIR}/lib/Geant4-10.1.2
                ${GEANT4_DIR}/lib/Geant4-10.1.3
                ${GEANT4_DIR}/lib/Geant4-10.2.0
          PATHS ENV LD_LIBRARY_PATH
          PATH_SUFFIXES lib Lib
          NO_DEFAULT_PATH)

message(STATUS "Found Geant4 in ${GEANT4_CMAKE_CONFIG}")

include(${GEANT4_CMAKE_CONFIG}/Geant4Config.cmake)
