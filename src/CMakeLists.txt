# Required libraries
add_subdirectory(dagmc)
add_subdirectory(pyne)
add_subdirectory(uwuw)

# Non-required libraries
if (BUILD_TALLY)
  add_subdirectory(tally)
endif ()

# Physics code interfaces
if (BUILD_MCNP5 OR BUILD_MCNP6 OR BUILD_CI_TESTS)
  add_subdirectory(mcnp)
endif ()
if (BUILD_GEANT4)
  add_subdirectory(geant4)
endif ()
if (BUILD_FLUKA)
  add_subdirectory(fluka)
endif ()

# Tools
if (BUILD_ASTYLE)
  add_subdirectory(astyle)
endif ()
if (BUILD_BUILD_OBB)
  add_subdirectory(build_obb)
endif ()
if (BUILD_MAKE_WATERTIGHT)
  add_subdirectory(make_watertight)
endif ()
