project(DAGMC)
cmake_minimum_required(VERSION 3.18)
enable_language(CXX)

# Set DAGMC version
set(DAGMC_MAJOR_VERSION 3)
set(DAGMC_MINOR_VERSION 2)
set(DAGMC_PATCH_VERSION 4)
set(DAGMC_VERSION ${DAGMC_MAJOR_VERSION}.${DAGMC_MINOR_VERSION}.${DAGMC_PATCH_VERSION})

if(MSVC)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  add_definitions( -DH5_BUILT_AS_DYNAMIC_LIB )
endif()


# Set DAGMC Git SHA

# Set git SHA1 hash as a compile definition
find_package(Git)
if(GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE DAGMC_GIT_SHA_SUCCESS
    OUTPUT_VARIABLE DAGMC_GIT_SHA
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT DAGMC_GIT_SHA_SUCCESS EQUAL "0")
    message(WARNING "Could not determine the commit SHA for DAGMC.")
    set(DAGMC_GIT_SHA "")
  endif()

  option(GIT_SUBMODULE "Check submodules during build" ON)
endif()

# Make the scripts in the "cmake" directory available to CMake
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(DAGMC_macros)

dagmc_setup_build()

dagmc_setup_options()

if(BUILD_MCNP5 OR BUILD_MCNP6)
  enable_language(Fortran)
endif()

# Dealing with MOAB
if (PULL_INSTALL_MOAB)
  set(MOAB_VERSION ${PULL_INSTALL_MOAB})  

  # Ensure MOAB_VERSION
  if(NOT ${MOAB_VERSION} MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+")
    MESSAGE("MOAB_VERSION: ${MOAB_VERSION} not recognized. Defaulting to 5.5.1")
    set(MOAB_VERSION "5.5.1")
  endif()
  MESSAGE("MOAB_VERSION to be installed: ${MOAB_VERSION}")
  include(MOAB_PullAndMake)
  moab_pull_make(${MOAB_VERSION})
else()
  find_package(MOAB REQUIRED)
endif()

find_package(OpenMP)

dagmc_setup_flags()

if (BUILD_TESTS)
  enable_testing()
endif ()

add_subdirectory(src)

dagmc_make_configure_files()

message("")
