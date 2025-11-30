# Include the ExternalProject module
# This module handles downloading and building external projects
include(ExternalProject)

############################################################
# Setup HDF5                                               #
############################################################

# If HDF5 is installed in a non-standard location, set HDF5_DIR
# to the cmake directory of the HDF5 installation.
# For example: export HDF5_DIR=/path/to/hdf5/cmake
# If HDF5 is not installed, set DOWNLOAD_HDF5 to ON.
# PyNE will download and build HDF5 if DOWNLOAD_HDF5 is ON.
if(NOT DOWNLOAD_HDF5)
  find_package(HDF5 REQUIRED)
else()
  message(STATUS "HDF5 will be downloaded and installed")

  # Configure HDF5 
  if(NOT HDF5_VERSION)
    set(HDF5_VERSION "1.14.3")
  endif()
  
  # Replace the dots with semicolons to create a list
  string(REPLACE "." ";" HDF5_VERSION_LIST ${HDF5_VERSION})

  # Convert the list into individual variables
  list(GET HDF5_VERSION_LIST 0 HDF5_VERSION_MAJOR)
  list(GET HDF5_VERSION_LIST 1 HDF5_VERSION_MINOR)
  list(GET HDF5_VERSION_LIST 2 HDF5_VERSION_PATCH)
  set(HDF5_ROOT "${CMAKE_BINARY_DIR}/hdf5")
  set(HDF5_INCLUDE_DIRS "${HDF5_ROOT}/include")
  set(HDF5_LIBRARY_DIRS "${HDF5_ROOT}/lib")
  ExternalProject_Add(hdf5-project
    PREFIX ${HDF5_ROOT}
    URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5_VERSION_MAJOR}.${HDF5_VERSION_MINOR}/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
      -DBUILD_SHARED_LIBS:BOOL=ON
      -DBUILD_STATIC_LIBS:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DHDF5_BUILD_TOOLS:BOOL=OFF
      -DHDF5_BUILD_EXAMPLES:BOOL=OFF
    DOWNLOAD_EXTRACT_TIMESTAMP true
    BUILD_BYPRODUCTS "${HDF5_LIBRARY_DIRS}/*${CMAKE_SHARED_LIBRARY_SUFFIX}*"
    )

  # HDF5 is a shared library, so we need to install it
  install_dependent_library(hdf5 ${HDF5_LIBRARY_DIRS}) 
endif()

# Include the HDF5 libraries
include_directories(${HDF5_INCLUDE_DIRS})
message(STATUS "HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")

link_directories(${HDF5_LIBRARY_DIRS})
message(STATUS "HDF5 library dirs: ${HDF5_LIBRARY_DIRS}")

# Set up shared and static versions of the HDF5 libraries
set(HDF5_LIBRARIES_SHARED ${HDF5_LIBRARIES})
if(BUILD_STATIC_LIBS)
  string(REPLACE ${CMAKE_SHARED_LIBRARY_SUFFIX} ${CMAKE_STATIC_LIBRARY_SUFFIX}
          HDF5_LIBRARIES_STATIC "${HDF5_LIBRARIES_SHARED}"
  )
endif()
if(NOT BUILD_SHARED_LIBS)
  set(HDF5_LIBRARIES_SHARED)
endif()

# Print out some more information about the HDF5 setup
message(STATUS "HDF5 libraries shared: ${HDF5_LIBRARIES_SHARED}")
message(STATUS "HDF5 libraries static: ${HDF5_LIBRARIES_STATIC}")


############################################################
# Setup Eigen3                                             #
############################################################

# If Eigen3 is installed in a non-standard location, set EIGEN3_DIR
# to the cmake directory of the Eigen3 installation.
# For example: export EIGEN3_ROOT=/path/to/eigen3/cmake
# If Eigen is not installed, set DOWNLOAD_EIGEN3 to ON.
# PyNE will download Eigen3 if DOWNLOAD_EIGEN3 is ON. (Recommended)
if(NOT DOWNLOAD_EIGEN3)
  find_package(Eigen3 REQUIRED)
else()
  message(STATUS "Eigen3 will be downloaded")

  # Configure Eigen3
  if(NOT EIGEN3_VERSION)
    set(EIGEN3_VERSION "3.4.0")
  endif()
  set(EIGEN3_ROOT "${CMAKE_BINARY_DIR}/eigen3")
  set(EIGEN3_INCLUDE_DIRS "${EIGEN3_ROOT}/include/eigen3")
  ExternalProject_Add(eigen3-project
    PREFIX ${EIGEN3_ROOT}
    URL https://gitlab.com/libeigen/eigen/-/archive/${EIGEN3_VERSION}/eigen-${EIGEN3_VERSION}.tar.bz2
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    DOWNLOAD_EXTRACT_TIMESTAMP true
  )
endif()

# Include the Eigen3 libraries
include_directories("${EIGEN3_INCLUDE_DIRS}")
message(STATUS "Eigen3 include dirs: ${EIGEN3_INCLUDE_DIRS}")


############################################################
# Setup MOAB                                               #
############################################################

# If MOAB is installed in a non-standard location, set MOAB_DIR
# to the cmake directory of the MOAB installation.
# For example: export MOAB_ROOT=/path/to/moab/cmake
# If MOAB is not installed, set DOWNLOAD_MOAB to ON.
# PyNE will download and build MOAB if DOWNLOAD_MOAB is ON. (Recommended)
if(NOT DOWNLOAD_MOAB)
  find_package(MOAB REQUIRED)
else()
  message(STATUS "MOAB will be downloaded and built")

  # Configure MOAB
  if(NOT MOAB_VERSION)
    set(MOAB_VERSION "5.6.0")
  endif()
  set(MOAB_ROOT "${CMAKE_BINARY_DIR}/moab")
  set(MOAB_INCLUDE_DIRS "${MOAB_ROOT}/include")
  set(MOAB_LIBRARY_DIRS "${MOAB_ROOT}/${CMAKE_INSTALL_LIBDIR}")
  ExternalProject_Add(moab-project
    PREFIX ${MOAB_ROOT}
    GIT_REPOSITORY https://bitbucket.org/fathomteam/moab.git
    GIT_TAG ${MOAB_VERSION}
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
      -DBUILD_SHARED_LIBS:BOOL=ON
      -DENABLE_HDF5:BOOL=ON
      -DHDF5_ROOT:PATH=${HDF5_ROOT}
      -DENABLE_BLASLAPACK:BOOL=OFF
      -DENABLE_FORTRAN:BOOL=OFF
      -DCMAKE_MACOSX_RPATH:BOOL=ON
    DOWNLOAD_EXTRACT_TIMESTAMP true
    BUILD_BYPRODUCTS "${MOAB_LIBRARY_DIRS}/*${CMAKE_SHARED_LIBRARY_SUFFIX}*"
  )

  # HDF5 is needed to be build before MOAB
  if(DOWNLOAD_HDF5)
    add_dependencies(moab-project hdf5-project)
  endif()

  # EIGEN3 is needed to be build before MOAB
  if(DOWNLOAD_EIGEN3)
    add_dependencies(moab-project eigen3-project)
  endif()

  # MOAB is a shared library, so we need to install it
  install_dependent_library(MOAB ${MOAB_LIBRARY_DIRS})
endif()
include_directories(${MOAB_INCLUDE_DIRS})
link_directories(${MOAB_LIBRARY_DIRS})
message(STATUS "MOAB include dirs: ${MOAB_INCLUDE_DIRS}")
message(STATUS "MOAB library dirs: ${MOAB_LIBRARY_DIRS}")
