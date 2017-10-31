message(STATUS "Building the DagSolid Library for Geant4")

find_package(Geant4 REQUIRED)

set(SRC_FILES DagSolid.cc)

set(PUB_HEADERS DagSolid.hh)

include_directories(${CMAKE_CURRENT_LIST_DIR})
include_directories(${CMAKE_SOURCE_DIR}/dagmc)

set(LINK_LIBS_SHARED dagmc        ${Geant4_LIBRARIES_SHARED})
set(LINK_LIBS_STATIC dagmc-static ${Geant4_LIBRARIES_STATIC})

add_library(dagsolid        SHARED ${SRC_FILES})
add_library(dagsolid-static STATIC ${SRC_FILES})
set_target_properties(dagsolid        PROPERTIES PUBLIC_HEADER "${PUB_HEADERS}")
set_target_properties(dagsolid-static PROPERTIES OUTPUT_NAME dagsolid)
target_link_libraries(dagsolid        ${LINK_LIBS_SHARED})
target_link_libraries(dagsolid-static ${LINK_LIBS_STATIC})
install(TARGETS dagsolid        LIBRARY       DESTINATION lib
                                PUBLIC_HEADER DESTINATION include)
install(TARGETS dagsolid-static ARCHIVE       DESTINATION lib)

add_subdirectory(build)
add_subdirectory(tests)
