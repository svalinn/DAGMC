message(STATUS "Building astyle_google")

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

# Build astyle_google
add_executable(astyle_google ${SRC_FILES})
install(TARGETS astyle_google DESTINATION tools
                              PERMISSIONS ${PERMS})
install(FILES file/google.ini DESTINATION tools)
