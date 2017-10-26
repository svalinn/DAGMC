message(STATUS "Building DAGMC tools")

macro (build_tool name ext libs)
  add_executable(${name} ${name}.${ext})
  target_link_libraries(${name} ${libs})
  install(TARGETS ${name} DESTINATION tools)
endmacro (build_tool)

# Libraries
set(LINK_LIBS_SHARED dagmc)
set(LINK_LIBS_STATIC dagmc-static)

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_EXE ${LINK_LIBS_STATIC})
else ()
  set(LINK_LIBS_EXE ${LINK_LIBS_SHARED})
endif ()

build_tool(pt_vol_test   cpp "${LINK_LIBS_EXE}")
build_tool(ray_fire_test cpp "${LINK_LIBS_EXE}")
build_tool(test_geom     cpp "${LINK_LIBS_EXE}")
