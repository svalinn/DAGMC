project(DAGMC)
cmake_minimum_required(VERSION 3.1)
enable_language(CXX)

# Set DAGMC version
set(DAGMC_MAJOR_VERSION 3)
set(DAGMC_MINOR_VERSION 2)
set(DAGMC_PATCH_VERSION 0)
set(DAGMC_VERSION ${DAGMC_MAJOR_VERSION}.${DAGMC_MINOR_VERSION}.${DAGMC_PATCH_VERSION})

# Set DAGMC Git SHA

# Set git SHA1 hash as a compile definition
find_package(Git)
if(GIT_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE DAGMC_GIT_SHA_SUCCESS
    OUTPUT_VARIABLE DAGMC_GIT_SHA
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT DAGMC_GIT_SHA_SUCCESS EQUAL "0")
    message(WARNING "Could not determine the commit SHA for DAGMC.")
    set(DAGMC_GIT_SHA "")
  endif()
endif()

# Make the scripts in the "cmake" directory available to CMake
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(DAGMC_macros)

dagmc_setup_build()

dagmc_setup_options()

if(BUILD_MCNP5 OR BUILD_MCNP6)
  enable_language(Fortran)
endif()

find_package(MOAB REQUIRED)
find_package(OpenMP)

dagmc_setup_flags()

if (BUILD_TESTS)
  enable_testing()
endif ()

add_subdirectory(src)

dagmc_make_configure_files()

if (COVERALLS)
  # Create the coveralls target.
  coveralls_setup(
      "${COVERAGE_SRCS}" # The source files.
      ON                 # If we should upload.
      "${PROJECT_SOURCE_DIR}/cmake") # (Optional) Alternate project cmake module path.
endif()


message("")
