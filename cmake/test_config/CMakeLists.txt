project(Test_DAGMC_config)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

find_package(DAGMC REQUIRED)

link_directories(${DAGMC_LIBRARY_DIRS})
include_directories(${DAGMC_INCLUDE_DIRS})

add_executable(test_DAGMCConfig test_DAGMCConfig.cpp)
target_link_libraries(test_DAGMCConfig ${DAGMC_LIBRARIES})
add_test(NAME DAGMC-Configuration-Test COMMAND test_DAGMCConfig)
enable_testing()
