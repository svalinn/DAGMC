message(STATUS "Building DAGMC tools")

macro (build_tool name ext libs)
  add_executable(${name} ${name}.${ext})
  target_link_libraries(${name} ${libs})
  install(TARGETS ${name} DESTINATION tools)
endmacro (build_tool)

# Libraries
set(LINK_LIBS)
if (BUILD_STATIC)
  list(APPEND LINK_LIBS dagmc-static)
else ()
  list(APPEND LINK_LIBS dagmc)
endif (BUILD_STATIC)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

build_tool(pt_vol_test   cpp "${LINK_LIBS}")
build_tool(ray_fire_test cpp "${LINK_LIBS}")
build_tool(test_geom     cpp "${LINK_LIBS}")
