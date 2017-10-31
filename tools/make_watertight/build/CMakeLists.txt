set(LINK_LIBS_SHARED makeWatertight)
set(LINK_LIBS_STATIC makeWatertight-static)

if (BUILD_STATIC_EXE)
  set(LINK_LIBS_EXE ${LINK_LIBS_STATIC})
else ()
  set(LINK_LIBS_EXE ${LINK_LIBS_SHARED})
endif ()

include_directories(${CMAKE_SOURCE_DIR}/tools/make_watertight)

set(SRC_FILES make_watertight.cpp)
dagmc_install_exe(make_watertight)
set(SRC_FILES check_watertight.cpp)
dagmc_install_exe(check_watertight)
