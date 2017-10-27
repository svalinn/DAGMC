message(STATUS "Building DAGMC library")

set(SRC_FILES
  DagMC.cpp
  dagmcmetadata.cpp
  ${CMAKE_SOURCE_DIR}/uwuw/uwuw.cpp
  ${CMAKE_SOURCE_DIR}/uwuw/uwuw_preprocessor.cpp
)

set(PUB_HEADERS
  DagMC.hpp
  dagmcmetadata.hpp
  ${CMAKE_SOURCE_DIR}/uwuw/uwuw.hpp
  ${CMAKE_SOURCE_DIR}/uwuw/uwuw_preprocessor.hpp
)

set(LINK_LIBS_SHARED
  ${HDF5_LIBRARIES_SHARED}
  ${MOAB_LIBRARIES_SHARED}
  pyne_dagmc
)
set(LINK_LIBS_STATIC
  ${HDF5_LIBRARIES_STATIC}
  ${MOAB_LIBRARIES_STATIC}
  pyne_dagmc-static
)

include_directories(${CMAKE_CURRENT_LIST_DIR})

add_library(dagmc        SHARED ${SRC_FILES})
add_library(dagmc-static STATIC ${SRC_FILES})
set_target_properties(dagmc        PROPERTIES PUBLIC_HEADER "${PUB_HEADERS}")
set_target_properties(dagmc-static PROPERTIES OUTPUT_NAME dagmc)
install(TARGETS dagmc        LIBRARY       DESTINATION lib
                             PUBLIC_HEADER DESTINATION include)
install(TARGETS dagmc-static ARCHIVE       DESTINATION lib)
target_link_libraries(dagmc        ${LINK_LIBS_SHARED})
target_link_libraries(dagmc-static ${LINK_LIBS_STATIC})

add_subdirectory(tools)
add_subdirectory(tests)
