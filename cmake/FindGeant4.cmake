# Try to find Geant4
#
# Once done this will define
#
#  Geant4_USE_FILE - CMake config for Geant4
#  Geant4_INCLUDE_DIRS - Geant4 include directories
#  Geant4_LIBRARIES - In case they are ever needed

# Try to find Geant4Config.cmake
file(GLOB SEARCH_DIRS "${GEANT4_DIR}/lib*/Geant4-*")
string(REPLACE "\n" ";" SEARCH_DIRS ${SEARCH_DIRS})
find_path(GEANT4_CMAKE_CONFIG
          NAMES Geant4Config.cmake
          PATHS ${SEARCH_DIRS}
          NO_DEFAULT_PATH)
if (${GEANT4_CMAKE_CONFIG} MATCHES GEANT4_CMAKE_CONFIG-NOTFOUND)
  message(FATAL_ERROR "Geant4Config.cmake not found. Make sure to set -DGEANT4_DIR when running CMake.")
endif (${GEANT4_CMAKE_CONFIG} MATCHES GEANT4_CMAKE_CONFIG-NOTFOUND)
message(STATUS "Found Geant4Config.cmake at ${GEANT4_CMAKE_CONFIG}")

# Geant4 CMake config
include(${GEANT4_CMAKE_CONFIG}/Geant4Config.cmake)
