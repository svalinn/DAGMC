# This file isn't named "FindHDF5.cmake" because that file already exists in CMake

message("")

if (HDF5_ROOT AND NOT HDF5_DIR)
  set(HDF5_DIR ${HDF5_ROOT})
endif ()

# Use the built-in FindHDF5 script
set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
find_package(HDF5 REQUIRED)
set(HDF5_LIBRARIES_SHARED ${HDF5_LIBRARIES})
# CMake doesn't let you find_package(HDF5) twice so we have to do this instead
string(REPLACE ".so" ".a" HDF5_LIBRARIES_STATIC "${HDF5_LIBRARIES_SHARED}")
set(HDF5_LIBRARIES)

message(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5_LIBRARIES_SHARED: ${HDF5_LIBRARIES_SHARED}")
message(STATUS "HDF5_LIBRARIES_STATIC: ${HDF5_LIBRARIES_STATIC}")

include_directories(${HDF5_INCLUDE_DIRS})

message("")
