message(STATUS "Including UWUW in dagmciface library")

set(CODE_SRC_DIR     "uwuw")
set(CODE_SRC_FILES   "uwuw.cpp")
set(CODE_PUB_HEADERS "uwuw.hpp")

set(DI_PUB_HEADERS ${DI_PUB_HEADERS} ${CMAKE_CURRENT_LIST_DIR}/${CODE_PUB_HEADERS} PARENT_SCOPE)
set(DI_SRC_FILES   ${DI_SRC_FILES}   ${CMAKE_CURRENT_LIST_DIR}/${CODE_SRC_FILES}   PARENT_SCOPE)

add_subdirectory(tests)
