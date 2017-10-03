option(BUILD_BUILD_OBB "Build build_obb" ON)
option(BUILD_MAKE_WATERTIGHT "Build make_watertight" ON)

if (BUILD_BUILD_OBB)
  add_subdirectory(build_obb)
endif (BUILD_BUILD_OBB)
if (BUILD_MAKE_WATERTIGHT)
  add_subdirectory(make_watertight)
endif (BUILD_MAKE_WATERTIGHT)
