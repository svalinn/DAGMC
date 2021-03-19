message("")

find_package(Git)
  if (GIT_FOUND)
    message(STATUS "Updating PyNE submodule")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --force
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if (NOT GIT_SUBMOD_RESULT EQUAL 0)
      message(FATAL_ERROR "git submodule update --init failed with \
          ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif ()
  endif ()


MESSAGE("${CMAKE_CURRENT_SOURCE_DIR}/amalgamate_pyne.sh")
execute_process(COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/amalgamate_pyne.sh ${CMAKE_CURRENT_SOURCE_DIR}/pyne ${CMAKE_CURRENT_BINARY_DIR} RESULT_VARIABLE STATUS  OUTPUT_VARIABLE OUTPUT2)
if(STATUS AND NOT STATUS EQUAL 0)
  message(STATUS "PyNE AMALGAMATION FAILED: ${STATUS}")
else()
  message(STATUS "PyNE AMALGAMATION SUCCESSFULL ${OUTPUT2}")
endif()

set(SRC_FILES ${CMAKE_CURRENT_BINARY_DIR}/pyne.cpp)
set(PUB_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/pyne.h)

set(LINK_LIBS)
set(LINK_LIBS_EXTERN_NAMES HDF5_LIBRARIES MOAB_LIBRARIES)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE        "-O0")
  set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0 -g")
endif ()

if(NOT (${HDF5_VERSION} VERSION_LESS 1.12.0))
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DH5Oget_info_by_idx_vers=1 -DH5O_info_t_vers=1 -DH5Oget_info_by_name_vers=1")
 endif()

dagmc_install_library(pyne_dagmc)
