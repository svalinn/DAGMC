set(DAGMC_LIBRARY_LIST dagmc pyne_dagmc uwuw dagtally makeWatertight dagsolid fludag)

macro (dagmc_set_build_type)
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
  message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endmacro ()

macro (dagmc_setup_build_options)
  option(BUILD_MCNP5  "Build DAG-MCNP5"                       OFF)
  option(BUILD_MCNP6  "Build DAG-MCNP6"                       OFF)
  option(MCNP5_PLOT   "Build DAG-MCNP5 with plotting support" OFF)
  option(MCNP6_PLOT   "Build DAG-MCNP6 with plotting support" OFF)
  option(MPI_BUILD    "Build DAG-MCNP5/6 with MPI support"    OFF)
  option(OPENMP_BUILD "Build DAG-MCNP5/6 with OpenMP support" OFF)

  option(BUILD_GEANT4      "Build DAG-Geant4" OFF)
  option(WITH_GEANT4_UIVIS "Build DAG-Geant4 with visualization support" ${BUILD_GEANT4})

  option(BUILD_FLUKA "Build FluDAG" OFF)

  option(BUILD_TALLY "Build dagtally library" ON)

  option(BUILD_ASTYLE          "Build astyle code formatter" ON)
  option(BUILD_BUILD_OBB       "Build build_obb tool"        ON)
  option(BUILD_MAKE_WATERTIGHT "Build make_watertight tool"  ON)

  option(BUILD_TESTS    "Build unit tests" ON)
  option(BUILD_CI_TESTS "Build everything needed to run the CI tests" ${BUILD_TESTS})

  option(BUILD_STATIC_EXE "Build static executables" OFF)

  if (BUILD_ALL)
    set(BUILD_MCNP5  ON)
    set(BUILD_MCNP6  ON)
    set(BUILD_GEANT4 ON)
    set(BUILD_FLUKA  ON)
  endif ()
endmacro ()

macro (dagmc_setup_flags)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)

  set(CXX_LIBRARY)
  foreach (library IN LISTS CMAKE_CXX_IMPLICIT_LINK_LIBRARIES)
    if (${library} MATCHES "c\\+\\+")
      set(CXX_LIBRARY ${library})
      break()
    endif ()
  endforeach ()
  message(STATUS "CXX_LIBRARY: ${CXX_LIBRARY}")

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

  if (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-intel")
  endif ()

  message(STATUS "CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")
endmacro ()

macro (dagmc_setup_rpath)
  set(INSTALL_RPATH_DIRS ${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR})

  get_filename_component(CXX_COMPILER_ROOT ${CMAKE_CXX_COMPILER} DIRECTORY)
  get_filename_component(CXX_COMPILER_ROOT ${CXX_COMPILER_ROOT} DIRECTORY)
  find_path(CXX_LIBRARY_PATH
    NAMES "lib${CXX_LIBRARY}.so"
    HINTS ${CXX_COMPILER_ROOT}
    PATH_SUFFIXES lib64 lib32 lib
  )
  get_filename_component(Fortran_COMPILER_ROOT ${CMAKE_Fortran_COMPILER} DIRECTORY)
  get_filename_component(Fortran_COMPILER_ROOT ${Fortran_COMPILER_ROOT} DIRECTORY)
  find_path(Fortran_LIBRARY_PATH
    NAMES "libgfortran.so"
    HINTS ${Fortran_COMPILER_ROOT}
    PATH_SUFFIXES lib64 lib32 lib
  )

  if (Fortran_LIBRARY_PATH AND (NOT Fortran_LIBRARY_PATH STREQUAL CXX_LIBRARY_PATH))
    set(INSTALL_RPATH_DIRS "${Fortran_LIBRARY_PATH}:${INSTALL_RPATH_DIRS}")
  endif ()
  if (CXX_LIBRARY_PATH)
    set(INSTALL_RPATH_DIRS "${CXX_LIBRARY_PATH}:${INSTALL_RPATH_DIRS}")
  endif ()

  message(STATUS "INSTALL_RPATH_DIRS: ${INSTALL_RPATH_DIRS}")
endmacro ()

# Figure out what LINK_LIBS_SHARED and LINK_LIBS_STATIC should be based on the
# values of LINK_LIBS and LINK_LIBS_EXTERN_NAMES
macro (dagmc_get_link_libs)
  set(LINK_LIBS_SHARED ${LINK_LIBS})
  set(LINK_LIBS_STATIC)

  foreach (extern_name IN LISTS LINK_LIBS_EXTERN_NAMES)
    list(APPEND LINK_LIBS_SHARED ${${extern_name}_SHARED})
    list(APPEND LINK_LIBS_STATIC ${${extern_name}_STATIC})
  endforeach ()

  foreach (link_lib IN LISTS LINK_LIBS)
    list(FIND DAGMC_LIBRARY_LIST ${link_lib} index)
    if (index STREQUAL "-1")
      list(APPEND LINK_LIBS_STATIC ${link_lib})
    else ()
      list(APPEND LINK_LIBS_STATIC ${link_lib}-static)
    endif ()
  endforeach ()
endmacro ()

# To use the dagmc_install macros, the following lists must be defined:
#   SRC_FILES: source files
#   PUB_HEADERS: public header files
#   LINK_LIBS: e.g. dagmc, pyne_dagmc, uwuw, lapack, gfortran
#   LINK_LIBS_EXTERN_NAMES: e.g. HDF5_LIBRARIES, MOAB_LIBRARIES

# Install a library in both shared and static mode
macro (dagmc_install_library lib_name)
  message(STATUS "Building library: ${lib_name}")

  dagmc_get_link_libs()

  add_library(${lib_name}        SHARED ${SRC_FILES})
  add_library(${lib_name}-static STATIC ${SRC_FILES})
  set_target_properties(${lib_name}
    PROPERTIES OUTPUT_NAME ${lib_name}
               PUBLIC_HEADER "${PUB_HEADERS}"
               INSTALL_RPATH ${INSTALL_RPATH_DIRS}
               INSTALL_RPATH_USE_LINK_PATH TRUE)
  set_target_properties(${lib_name}-static
    PROPERTIES OUTPUT_NAME ${lib_name}
               INSTALL_RPATH ""
               INSTALL_RPATH_USE_LINK_PATH FALSE)
  target_link_libraries(${lib_name}        ${LINK_LIBS_SHARED})
  target_link_libraries(${lib_name}-static ${LINK_LIBS_STATIC})
  install(TARGETS ${lib_name}
    LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
    PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR})
  install(TARGETS ${lib_name}-static
    ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
    PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR})
endmacro ()

# Install an executable
macro (dagmc_install_exe exe_name)
  message(STATUS "Building executable: ${exe_name}")

  dagmc_get_link_libs()

  add_executable(${exe_name} ${SRC_FILES})
  if (BUILD_STATIC_EXE)
    set_target_properties(${exe_name}
      PROPERTIES INSTALL_RPATH ""
                 INSTALL_RPATH_USE_LINK_PATH FALSE)
    target_link_libraries(${exe_name} ${LINK_LIBS_STATIC})
  else ()
    set_target_properties(${exe_name}
      PROPERTIES INSTALL_RPATH ${INSTALL_RPATH_DIRS}
                 INSTALL_RPATH_USE_LINK_PATH TRUE)
    target_link_libraries(${exe_name} ${LINK_LIBS_SHARED})
  endif ()
  install(TARGETS ${exe_name} DESTINATION ${INSTALL_BIN_DIR})
endmacro ()

# Install a unit test
macro (dagmc_install_test test_name ext)
  message(STATUS "Building unit tests: ${test_name}")

  list(APPEND LINK_LIBS gtest)

  dagmc_get_link_libs()

  add_executable(${test_name} ${test_name}.${ext} ${DRIVERS})
  if (BUILD_STATIC_EXE)
    set_target_properties(${test_name}
      PROPERTIES INSTALL_RPATH ""
                 INSTALL_RPATH_USE_LINK_PATH FALSE)
    target_link_libraries(${test_name} ${LINK_LIBS_STATIC})
  else ()
    set_target_properties(${test_name}
      PROPERTIES INSTALL_RPATH ${INSTALL_RPATH_DIRS}
                 INSTALL_RPATH_USE_LINK_PATH TRUE)
    target_link_libraries(${test_name} ${LINK_LIBS_SHARED})
  endif ()
  install(TARGETS ${test_name} DESTINATION ${INSTALL_TESTS_DIR})
endmacro ()
