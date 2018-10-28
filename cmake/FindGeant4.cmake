# Find Geant4 cmake config file
set(Geant4_SEARCH_DIRS)
file(GLOB Geant4_SEARCH_DIRS "${GEANT4_DIR}/lib*/Geant4-*")
string(REPLACE "\n" ";" Geant4_SEARCH_DIRS ${Geant4_SEARCH_DIRS})
find_path(Geant4_CMAKE_CONFIG
  NAMES Geant4Config.cmake
  PATHS ${Geant4_SEARCH_DIRS}
  NO_DEFAULT_PATH
)
if (Geant4_CMAKE_CONFIG)
  set(Geant4_CMAKE_CONFIG ${Geant4_CMAKE_CONFIG}/Geant4Config.cmake)
  message(STATUS "Geant4_CMAKE_CONFIG: ${Geant4_CMAKE_CONFIG}")
else ()
  message(FATAL_ERROR "Could not find Geant4")
endif ()

# Activate all available UI and Vis drivers by default
if (WITH_GEANT4_UIVIS)
  set(Geant4_FIND_REQUIRED_ui_all ON)
  set(Geant4_FIND_REQUIRED_vis_all ON)
endif ()

# Get Geant4 include directories and libraries
set(Geant4_FIND_REQUIRED_static ON)
include(${Geant4_CMAKE_CONFIG})
set(Geant4_LIBRARIES_STATIC ${Geant4_LIBRARIES})
set(Geant4_FIND_REQUIRED_static OFF)
include(${Geant4_CMAKE_CONFIG})
set(Geant4_LIBRARIES_SHARED ${Geant4_LIBRARIES})
set(Geant4_LIBRARIES)

# Get Geant4 compile definitions and flags
set(CMAKE_EXE_LINKER_FLAGS_SAVE ${CMAKE_EXE_LINKER_FLAGS})
include(${Geant4_USE_FILE})
set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS_SAVE})

message(STATUS "Geant4_INCLUDE_DIRS: ${Geant4_INCLUDE_DIRS}")
message(STATUS "Geant4_LIBRARIES_SHARED: ${Geant4_LIBRARIES_SHARED}")
message(STATUS "Geant4_LIBRARIES_STATIC: ${Geant4_LIBRARIES_STATIC}")
message(STATUS "Geant4_DEFINITIONS: ${Geant4_DEFINITIONS}")
message(STATUS "Geant4_CXX_FLAGS: ${Geant4_CXX_FLAGS}")
message(STATUS "Geant4_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}: ${Geant4_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
message(STATUS "Geant4_EXE_LINKER_FLAGS: ${Geant4_EXE_LINKER_FLAGS}")

if (Geant4_INCLUDE_DIRS AND Geant4_LIBRARIES_SHARED AND Geant4_LIBRARIES_STATIC AND Geant4_CXX_FLAGS)
  message(STATUS "Found Geant4")
else ()
  message(FATAL_ERROR "Could not find Geant4")
endif ()
