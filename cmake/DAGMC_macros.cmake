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

  option(BUILD_CI_TESTS "Build everything needed to run the CI tests" ON)

  option(BUILD_STATIC_EXE "Build static executables" OFF)

  if (BUILD_ALL)
    set(BUILD_MCNP5 ON)
    set(BUILD_MCNP6 ON)
    set(BUILD_GEANT4 ON)
    set(BUILD_FLUKA ON)
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
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # Don't include system directories
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  if ("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  endif ()

  message(STATUS "CMAKE_INSTALL_RPATH: ${CMAKE_INSTALL_RPATH}")
endmacro ()

macro (dagmc_setup_test test_name ext driver libs)
  set(drivers ${driver})
  # Convert driver from string to list
  string(STRIP "${drivers}" drivers)
  string(REGEX REPLACE "[ ]+" " " drivers "${drivers}")
  string(REPLACE " " ";" drivers "${drivers}")
  add_executable(${test_name} ${test_name}.${ext} ${drivers})
  target_link_libraries(${test_name} ${libs})
  install(TARGETS ${test_name} DESTINATION tests)
endmacro ()
