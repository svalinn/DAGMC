message(STATUS "Building DAG-Geant4")

# CXX compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Geant4_CXX_FLAGS}")

# Linker flags
if (BUILD_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
endif (BUILD_STATIC)

# Source files
set(SRC_FILES)
file(GLOB SRC_FILES ${SRC_FILES} "*.cc")
file(GLOB SRC_FILES ${SRC_FILES} "src/*.cc")

# Libraries
set(LINK_LIBS)
if (BUILD_STATIC)
  list(APPEND LINK_LIBS dagsolid-static)
  list(APPEND LINK_LIBS dagmc-static)
  list(APPEND LINK_LIBS pyne_dagmc-static)
  list(APPEND LINK_LIBS ${Geant4_LIBRARIES_STATIC})
else ()
  list(APPEND LINK_LIBS dagsolid)
  list(APPEND LINK_LIBS dagmc)
  list(APPEND LINK_LIBS pyne_dagmc)
  list(APPEND LINK_LIBS ${Geant4_LIBRARIES_SHARED})
endif (BUILD_STATIC)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Include directories
include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/uwuw)
include_directories(${CMAKE_SOURCE_DIR}/dagmc)

# Build DAG-Geant4
add_executable(DagGeant4 ${SRC_FILES})
target_link_libraries(DagGeant4 ${LINK_LIBS})
install(TARGETS DagGeant4 DESTINATION ${INSTALL_BIN_DIR})
