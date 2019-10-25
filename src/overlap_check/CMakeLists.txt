message("")

set(LINK_LIBS dagmc)
set(LINK_LIBS_EXTERN_NAMES)



if(OpenMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(app)

if (BUILD_TESTS)
  add_subdirectory(test)
endif()
