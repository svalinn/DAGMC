MACRO (find_set_HDF5)
  # Find HDF5
  set(ENV{PATH} "${HDF5_DIR}:$ENV{PATH}")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX})
  find_package(HDF5 REQUIRED)
  # Remove HDF5 transitive dependencies that are system libraries
  list(FILTER HDF5_LIBRARIES EXCLUDE REGEX ".*lib(pthread|dl|m).*")
  set(HDF5_LIBRARIES_SHARED ${HDF5_LIBRARIES})
  # CMake doesn't let you find_package(HDF5) twice so we have to do this instead
  if (BUILD_STATIC_LIBS)
    string(REPLACE ${CMAKE_SHARED_LIBRARY_SUFFIX} ${CMAKE_STATIC_LIBRARY_SUFFIX}
           HDF5_LIBRARIES_STATIC "${HDF5_LIBRARIES_SHARED}")
  endif ()
  if (NOT BUILD_SHARED_LIBS)
    set(HDF5_LIBRARIES_SHARED)
  endif ()
  set(HDF5_LIBRARIES)

  message(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
  message(STATUS "HDF5_LIBRARIES_SHARED: ${HDF5_LIBRARIES_SHARED}")
  message(STATUS "HDF5_LIBRARIES_STATIC: ${HDF5_LIBRARIES_STATIC}")

  include_directories(${HDF5_INCLUDE_DIRS})

ENDMACRO (find_set_HDF5)
