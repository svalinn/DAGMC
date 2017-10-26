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

set(LINK_LIBS_SHARED ${MOAB_LIBRARIES})
set(LINK_LIBS_STATIC ${MOAB_LIBRARIES})

add_library(makeWatertight        SHARED ${SRC_FILES})
add_library(makeWatertight-static STATIC ${SRC_FILES})
set_target_properties(makeWatertight        PROPERTIES PUBLIC_HEADER "${PUB_HEADERS}")
set_target_properties(makeWatertight-static PROPERTIES OUTPUT_NAME makeWatertight)
target_link_libraries(makeWatertight        ${LINK_LIBS_SHARED})
target_link_libraries(makeWatertight-static ${LINK_LIBS_STATIC})
install(TARGETS makeWatertight        LIBRARY       DESTINATION lib
                                      PUBLIC_HEADER DESTINATION include/make_watertight)
install(TARGETS makeWatertight-static ARCHIVE       DESTINATION lib)

add_subdirectory(build)
add_subdirectory(tests)
