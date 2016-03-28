message(STATUS "Including UWUW in dagmciface library")

# dagmciface source files
set(DI_SRC_FILES ${DI_SRC_FILES} ${CMAKE_CURRENT_LIST_DIR}/uwuw.cpp PARENT_SCOPE)
# if we want a static executable
IF(STATIC_LIB)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static")
ENDIF(STATIC_LIB)

CMAKE_POLICY(SET CMP0022 OLD)

# dagmciface public headers
set(DI_PUB_HEADERS ${DI_PUB_HEADERS} ${CMAKE_CURRENT_LIST_DIR}/uwuw.hpp PARENT_SCOPE)

set (UWUW_PREPROC_SRC "uwuw_preproc.cpp")
add_executable(uwuw_preproc ${UWUW_PREPROC_SRC})
target_link_libraries(uwuw_preproc dagmciface  
  ${DAG_LINK_LIBS} ${HDF5_hdf5_hl_LIBRARY} ${HDF5_hdf5_LIBRARY} 
  pyne_dagmc)

set(UWUW_PREPROCESSOR_SRC "uwuw_preprocessor.cpp")
add_library(uwuw_preprocessor ${UWUW_PREPROCESSOR_SRC})
target_link_libraries(uwuw_preprocessor dagmciface  
  ${DAG_LINK_LIBS} ${HDF5_hdf5_hl_LIBRARY} ${HDF5_hdf5_LIBRARY} 
  pyne_dagmc)

# Build UWUW tests
add_subdirectory(tests)

install(TARGETS uwuw_preproc DESTINATION bin)

