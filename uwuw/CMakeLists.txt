message(STATUS "Including UWUW in dagmciface library")

# dagmciface source files
set(DI_SRC_FILES ${DI_SRC_FILES} ${CMAKE_CURRENT_LIST_DIR}/uwuw.cpp PARENT_SCOPE)

# dagmciface public headers
set(DI_PUB_HEADERS ${DI_PUB_HEADERS} ${CMAKE_CURRENT_LIST_DIR}/uwuw.hpp PARENT_SCOPE)

# Build UWUW tests
add_subdirectory(tests)

set (UWUW_PREPROC_SRC "uwuw_preproc.cpp")
add_executable(uwuw_preproc ${UWUW_PREPROC_SRC})
target_link_libraries(uwuw_preproc dagmciface  
  ${DAG_LINK_LIBS} ${HDF5_hdf5_hl_LIBRARY} ${HDF5_hdf5_LIBRARY} 
  pyne_dagmc )
