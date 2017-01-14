cmake_minimum_required(VERSION 2.6)

project(FDTD)

if(NOT WIN32)
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra \
                        -fsanitize=undefined")
endif()

include_directories(${CMAKE_SOURCE_DIR}/include)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/Source)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/Error)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/Update)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/Process)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/PML)

add_executable(FDTD ${CMAKE_SOURCE_DIR}/src/main.cpp)
target_link_libraries(FDTD source error update process pml)
