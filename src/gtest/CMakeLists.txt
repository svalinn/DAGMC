message(STATUS "Building gtest")

set(GTEST_DIR gtest-1.7.0)
set(GTEST_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/${GTEST_DIR}/include PARENT_SCOPE)

# Source files
set(SRC_FILES ${GTEST_DIR}/src/gtest-all.cc)

# Include directories
include_directories(${GTEST_DIR})
include_directories(${GTEST_DIR}/include)

# CXX compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-as-needed")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-field-initializers")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions")

# Compiler directives
add_definitions(-DGTEST_HAS_PTHREAD=ON)

# Build gtest internal library
add_library(gtest STATIC ${SRC_FILES})
