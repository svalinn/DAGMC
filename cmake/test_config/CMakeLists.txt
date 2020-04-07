project(Test_DAGMC_config)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

# check for and enable c++11 support
  INCLUDE(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  IF(COMPILER_SUPPORTS_CXX11)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  ELSE()
      MESSAGE(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  ENDIF()

find_package(DAGMC REQUIRED)

link_directories(${DAGMC_LIBRARY_DIRS})
include_directories(${DAGMC_INCLUDE_DIRS})

add_executable(test_DAGMCConfig test_DAGMCConfig.cpp)
target_link_libraries(test_DAGMCConfig ${DAGMC_LIBRARIES})
add_test(NAME DAGMC-Configuration-Test COMMAND test_DAGMCConfig)
enable_testing()
