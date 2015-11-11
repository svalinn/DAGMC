# CMakeLists.txt for astyle
#
# get astyle

add_executable(astyle_google 
  ASBeautifier.cpp
  ASEnhancer.cpp
  ASFormatter.cpp
  ASLocalizer.cpp
  ASLocalizer.h
  ASResource.cpp
  CMakeLists.txt
  README.md
  astyle.h
  astyle_main.cpp
  astyle_main.h
)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS astyle_google DESTINATION tools)