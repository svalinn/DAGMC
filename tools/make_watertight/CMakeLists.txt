message(STATUS "Building make_watertight library")

set(SRC_FILES
  Gen.cpp
  Arc.cpp
  Zip.cpp
  Cleanup.cpp
  MakeWatertight.cpp
  CheckWatertight.cpp
)

set(PUB_HEADERS
  Gen.hpp
  Arc.hpp
  Zip.hpp
  Cleanup.hpp
  CheckWatertight.hpp
  MakeWatertight.hpp
)

set(LINK_LIBS_SHARED dagmc)
set(LINK_LIBS_STATIC dagmc-static)

dagmc_install_library(makeWatertight)

add_subdirectory(build)
add_subdirectory(tests)
