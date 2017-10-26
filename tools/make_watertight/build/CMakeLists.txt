message(STATUS "Building make_watertight and check_watertight")

set(LINK_LIBS_SHARED makeWatertight)
set(LINK_LIBS_STATIC makeWatertight-static)

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_EXE ${LINK_LIBS_STATIC})
else ()
  set(LINK_LIBS_EXE ${LINK_LIBS_SHARED})
endif ()

include_directories(${CMAKE_SOURCE_DIR}/tools/make_watertight)

add_executable(make_watertight  make_watertight.cpp)
add_executable(check_watertight check_watertight.cpp)
target_link_libraries(make_watertight  ${LINK_LIBS_EXE})
target_link_libraries(check_watertight ${LINK_LIBS_EXE})
install(TARGETS make_watertight  DESTINATION bin)
install(TARGETS check_watertight DESTINATION bin)
