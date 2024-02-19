

# Configure the documentation file

find_package(Doxygen)

set(DAGMC_SOURCE_DIR ${PROJECT_SOURCE_DIR})

if (DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Makefile ${CMAKE_CURRENT_BINARY_DIR}/Makefile COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/index.rst ${CMAKE_CURRENT_BINARY_DIR}/index.rst COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py ${CMAKE_CURRENT_BINARY_DIR}/conf.py COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt ${CMAKE_CURRENT_BINARY_DIR}/requirements.txt COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/make.bat ${CMAKE_CURRENT_BINARY_DIR}/make.bat COPYONLY)
else()
  message(WARNING "Doxygen is required to build the documentation for double-down.")
endif()