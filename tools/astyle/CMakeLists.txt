# astyle needs C++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if (COMPILER_SUPPORTS_CXX11)
  message(STATUS "Building astyle")
else ()
  message(STATUS "The C++ compiler does not support C++11. Not building astyle.")
  return()
endif (COMPILER_SUPPORTS_CXX11)

# Source files
set(SRC_FILES)
list(APPEND SRC_FILES src/ASBeautifier.cpp)
list(APPEND SRC_FILES src/ASEnhancer.cpp)
list(APPEND SRC_FILES src/ASFormatter.cpp)
list(APPEND SRC_FILES src/ASLocalizer.cpp)
list(APPEND SRC_FILES src/ASResource.cpp)
list(APPEND SRC_FILES src/astyle_main.cpp)

# Include directories
include_directories(src)

# C++ compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-rtti")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-exceptions")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Build astyle
add_executable(astyle ${SRC_FILES})
install(TARGETS astyle DESTINATION tools
                       PERMISSIONS ${PERMS})
install(FILES file/google.ini DESTINATION tools)
