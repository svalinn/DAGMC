set(GTEST_DIR gtest-1.7.0)
set(GTEST_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/${GTEST_DIR}/include PARENT_SCOPE)

set(SRC_FILES ${GTEST_DIR}/src/gtest-all.cc)

include_directories(${GTEST_DIR})
include_directories(${GTEST_DIR}/include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-as-needed")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-field-initializers")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions")

add_definitions(-DGTEST_HAS_PTHREAD=ON)

message(STATUS "Building library: gtest")
add_library(gtest STATIC ${SRC_FILES})
