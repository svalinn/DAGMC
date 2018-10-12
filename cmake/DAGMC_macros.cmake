# All DAGMC libraries
set(DAGMC_LIBRARY_LIST dagmc pyne_dagmc uwuw dagtally makeWatertight dagsolid fludag)

macro (dagmc_setup_build)
  message("")

  # Default to a release build
  if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "CMAKE_BUILD_TYPE not specified, defaulting to Release")
    set(CMAKE_BUILD_TYPE Release)
  endif ()
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Release" AND
      NOT CMAKE_BUILD_TYPE STREQUAL "Debug" AND
      NOT CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    message(FATAL_ERROR "Specified CMAKE_BUILD_TYPE is invalid; valid options are Release, Debug, RelWithDebInfo")
  endif ()
  string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
  message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

  # Installation directories
  set(INSTALL_BIN_DIR     bin)
  set(INSTALL_LIB_DIR     lib)
  set(INSTALL_INCLUDE_DIR include)
  set(INSTALL_TESTS_DIR   tests)
  set(INSTALL_TOOLS_DIR   tools)
  set(INSTALL_SHARE_DIR   share)

  # Get some environment variables
  set(ENV_USER "$ENV{USER}")
  execute_process(COMMAND hostname       OUTPUT_VARIABLE ENV_HOST OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND uname -s       OUTPUT_VARIABLE ENV_OS   OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND date +%m/%d/%y OUTPUT_VARIABLE ENV_DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND date +%H:%M:%S OUTPUT_VARIABLE ENV_TIME OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro ()

macro (dagmc_setup_options)
  message("")

  option(BUILD_MCNP5       "Build DAG-MCNP5"                         OFF)
  option(BUILD_MCNP6       "Build DAG-MCNP6"                         OFF)
  option(BUILD_MCNP_PLOT   "Build DAG-MCNP5/6 with plotting support" OFF)
  option(BUILD_MCNP_OPENMP "Build DAG-MCNP5/6 with OpenMP support"   OFF)
  option(BUILD_MCNP_MPI    "Build DAG-MCNP5/6 with MPI support"      OFF)

  option(BUILD_GEANT4      "Build DAG-Geant4" OFF)
  option(WITH_GEANT4_UIVIS "Build DAG-Geant4 with visualization support" ${BUILD_GEANT4})

  option(BUILD_FLUKA "Build FluDAG" OFF)

  option(BUILD_UWUW "Build UWUW library and uwuw_preproc" ON)
  option(BUILD_TALLY "Build dagtally library"              ON)

  option(BUILD_BUILD_OBB       "Build build_obb tool"       ON)
  option(BUILD_MAKE_WATERTIGHT "Build make_watertight tool" ON)

  option(BUILD_TESTS    "Build unit tests" ON)
  option(BUILD_CI_TESTS "Build everything needed to run the CI tests" OFF)

  option(BUILD_STATIC_EXE "Build static executables" OFF)
  option(BUILD_PIC        "Build with PIC"           OFF)

  option(BUILD_RPATH "Build libraries and executables with RPATH" ON)

  if (BUILD_ALL)
    set(BUILD_MCNP5  ON)
    set(BUILD_MCNP6  ON)
    set(BUILD_GEANT4 ON)
    set(BUILD_FLUKA  ON)
  endif ()
endmacro ()

macro (dagmc_setup_flags)
  message("")

  if (BUILD_PIC)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  endif ()

  set(CXX_LIBRARY)
  foreach (library IN LISTS CMAKE_CXX_IMPLICIT_LINK_LIBRARIES)
    if (library MATCHES "c\\+\\+")
      set(CXX_LIBRARY ${library})
      break()
    endif ()
  endforeach ()

  set(CMAKE_C_IMPLICIT_LINK_LIBRARIES         "")
  set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES       "")
  set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES       "${CXX_LIBRARY}")
  set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES     "")
  set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES   "")
  set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "")

  if (BUILD_STATIC_EXE)
    message(STATUS "Building static executables")
    set(BUILD_SHARED_EXE OFF)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
    set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
    set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
    set(CMAKE_EXE_LINK_DYNAMIC_Fortran_FLAGS)
  else ()
    message(STATUS "Building shared executables")
    set(BUILD_SHARED_EXE ON)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
  endif ()

  if (BUILD_RPATH)
    set(INSTALL_RPATH_DIRS "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}")
    message(STATUS "INSTALL_RPATH_DIRS: ${INSTALL_RPATH_DIRS}")
  endif ()

  message(STATUS "CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")
  message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
  message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
  message(STATUS "CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}: ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
  message(STATUS "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
  message(STATUS "CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}: ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
  message(STATUS "CMAKE_C_IMPLICIT_LINK_LIBRARIES: ${CMAKE_C_IMPLICIT_LINK_LIBRARIES}")
  message(STATUS "CMAKE_CXX_IMPLICIT_LINK_LIBRARIES: ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES}")
  message(STATUS "CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES: ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
  message(STATUS "CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")
endmacro ()

# Figure out what LINK_LIBS_SHARED and LINK_LIBS_STATIC should be based on the
# values of LINK_LIBS and LINK_LIBS_EXTERN_NAMES
macro (dagmc_get_link_libs)
  set(LINK_LIBS_SHARED)
  set(LINK_LIBS_STATIC)

  foreach (extern_name IN LISTS LINK_LIBS_EXTERN_NAMES)
    list(APPEND LINK_LIBS_SHARED ${${extern_name}_SHARED})
    list(APPEND LINK_LIBS_STATIC ${${extern_name}_STATIC})
  endforeach ()

  foreach (link_lib IN LISTS LINK_LIBS)
    list(FIND DAGMC_LIBRARY_LIST ${link_lib} index)
    if (index STREQUAL "-1")
      list(APPEND LINK_LIBS_SHARED ${link_lib})
      list(APPEND LINK_LIBS_STATIC ${link_lib})
    else ()
      list(APPEND LINK_LIBS_SHARED ${link_lib}-shared)
      list(APPEND LINK_LIBS_STATIC ${link_lib}-static)
    endif ()
  endforeach ()
endmacro ()

# Setup the configuration file and install
macro (dagmc_make_configure_file)
  message("")

  message(STATUS "DAGMC cmake config file: ${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}/cmake/DAGMCConfig.cmake")
  configure_file(cmake/DAGMCConfig.cmake.in DAGMCConfig.cmake @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/DAGMCConfig.cmake DESTINATION ${INSTALL_LIB_DIR}/cmake/)
endmacro ()

# To use the dagmc_install macros, the following lists must be defined:
#   SRC_FILES: source files
#   PUB_HEADERS: public header files
#   LINK_LIBS: e.g. dagmc, pyne_dagmc, uwuw, lapack, gfortran
#   LINK_LIBS_EXTERN_NAMES: e.g. MOAB_LIBRARIES

# Install a library in both shared and static mode
macro (dagmc_install_library lib_name)
  message(STATUS "Building library: ${lib_name}")

  dagmc_get_link_libs()

  add_library(${lib_name}-shared SHARED ${SRC_FILES})
  add_library(${lib_name}-static STATIC ${SRC_FILES})
  if (BUILD_RPATH)
    set_target_properties(${lib_name}-shared
      PROPERTIES OUTPUT_NAME ${lib_name}
                 PUBLIC_HEADER "${PUB_HEADERS}"
                 INSTALL_RPATH "${INSTALL_RPATH_DIRS}"
                 INSTALL_RPATH_USE_LINK_PATH TRUE)
    set_target_properties(${lib_name}-static
      PROPERTIES OUTPUT_NAME ${lib_name}
                 INSTALL_RPATH ""
                 INSTALL_RPATH_USE_LINK_PATH FALSE)
  else ()
    set_target_properties(${lib_name}-shared
      PROPERTIES OUTPUT_NAME ${lib_name}
                 PUBLIC_HEADER "${PUB_HEADERS}")
    set_target_properties(${lib_name}-static
      PROPERTIES OUTPUT_NAME ${lib_name})
  endif ()
  target_link_libraries(${lib_name}-shared ${LINK_LIBS_SHARED})
  target_link_libraries(${lib_name}-static ${LINK_LIBS_STATIC})
  install(TARGETS ${lib_name}-shared
    LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
    PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR})
  install(TARGETS ${lib_name}-static
    ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
    PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR})

  # Keep a list of all libraries being installed
  if (DAGMC_LIBRARIES)
    set(DAGMC_LIBRARIES "${DAGMC_LIBRARIES} ${lib_name}" CACHE INTERNAL "DAGMC_LIBRARIES")
  else ()
    set(DAGMC_LIBRARIES ${lib_name} CACHE INTERNAL "DAGMC_LIBRARIES")
  endif ()
endmacro ()

# Install an executable
macro (dagmc_install_exe exe_name)
  message(STATUS "Building executable: ${exe_name}")

  dagmc_get_link_libs()

  add_executable(${exe_name} ${SRC_FILES})
  if (BUILD_RPATH)
    if (BUILD_STATIC_EXE)
      set_target_properties(${exe_name}
        PROPERTIES INSTALL_RPATH ""
                   INSTALL_RPATH_USE_LINK_PATH FALSE)
      target_link_libraries(${exe_name} ${LINK_LIBS_STATIC})
    else ()
      set_target_properties(${exe_name}
        PROPERTIES INSTALL_RPATH "${INSTALL_RPATH_DIRS}"
                   INSTALL_RPATH_USE_LINK_PATH TRUE)
      target_link_libraries(${exe_name} ${LINK_LIBS_SHARED})
    endif ()
  else ()
    if (BUILD_STATIC_EXE)
      target_link_libraries(${exe_name} ${LINK_LIBS_STATIC})
    else ()
      target_link_libraries(${exe_name} ${LINK_LIBS_SHARED})
    endif ()
  endif ()
  install(TARGETS ${exe_name} DESTINATION ${INSTALL_BIN_DIR})
endmacro ()

# Install a unit test
macro (dagmc_install_test test_name ext)
  message(STATUS "Building unit tests: ${test_name}")

  list(APPEND LINK_LIBS gtest)

  dagmc_get_link_libs()

  add_executable(${test_name} ${test_name}.${ext} ${DRIVERS})
  if (BUILD_RPATH)
    if (BUILD_STATIC_EXE)
      set_target_properties(${test_name}
        PROPERTIES INSTALL_RPATH ""
                   INSTALL_RPATH_USE_LINK_PATH FALSE)
      target_link_libraries(${test_name} ${LINK_LIBS_STATIC})
    else ()
      set_target_properties(${test_name}
        PROPERTIES INSTALL_RPATH "${INSTALL_RPATH_DIRS}"
                   INSTALL_RPATH_USE_LINK_PATH TRUE)
      target_link_libraries(${test_name} ${LINK_LIBS_SHARED})
    endif ()
  else ()
    if (BUILD_STATIC_EXE)
      target_link_libraries(${test_name} ${LINK_LIBS_STATIC})
    else ()
      target_link_libraries(${test_name} ${LINK_LIBS_SHARED})
    endif ()
  endif ()
  install(TARGETS ${test_name} DESTINATION ${INSTALL_TESTS_DIR})
  add_test(NAME ${test_name} COMMAND ${test_name})
endmacro ()

# Install a file needed for unit testing
macro (dagmc_install_test_file filename)
  install(FILES ${filename} DESTINATION ${INSTALL_TESTS_DIR})
  install(FILES ${filename} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endmacro ()
