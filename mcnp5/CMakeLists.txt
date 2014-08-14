CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
MESSAGE ( STATUS "Including MCNP5 in dagmciface library.")

SET ( CODE_SRCDIR      "mcnp5" )
SET ( CODE_SRC_FILES   "${CODE_SRCDIR}/mcnp_funcs.cpp"  )
SET ( CODE_PUB_HEADERS "${CODE_SRCDIR}/mcnp_funcs.h" )

SET ( DI_PUB_HEADERS ${DI_PUB_HEADERS}
                     ${CODE_PUB_HEADERS}
                     ${CODE_CONFIG}
      PARENT_SCOPE)

SET ( DI_SRC_FILES ${DI_SRC_FILES} 
                   ${CODE_SRC_FILES} 
      PARENT_SCOPE)

# Copy MCNP source to build dir
#add_custom_command(OUTPUT mcnp5_out
  # make a directory for mcnp source
  # COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/Source 
  # copy mcnp source
  #  COMMAND ${CMAKE_COMMAND} -E copy_directory ${MCNP5_DIR}/Source ${CMAKE_CURRENT_SOURCE_DIR}/Source
  # move to the source directory
  # COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_SOURCE_DIR}/Source
  # why
  # execute_process (COMMAND patch "-p1 < ${CMAKE_CURRENT_SOURCE_DIR}/patch/mcnp5.1.60.patch")
#  COMMAND ${CMAKE_COMMAND} -E copy_directory ${MCNP5_DIR}/Source ${CMAKE_CURRENT_SOURCE_DIR}/Source
#  COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_SOURCE_DIR} patch -p1 < ${CMAKE_CURRENT_SOURCE_DIR}/patch/dagmc.patch.5.1.60
  
#)

# build MCNP
# add_custom_target(mcnp5_exec DEPENDS main.F90)
# add_custom_target(mcnp5 
#  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/build/CMakeLists.txt
# )
include(${CMAKE_CURRENT_SOURCE_DIR}/build/CMakeLists.txt)
# include(build/CmakeLists.txt)

