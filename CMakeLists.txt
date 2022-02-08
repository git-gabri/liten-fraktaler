cmake_minimum_required(VERSION 3.0.0)
project(liten-fraktaler)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -pedantic -g -O2 -pg")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -pedantic -s -O2 -msse4.1 -msse4.2 -mavx -mavx2")

file (GLOB_RECURSE LF_SOURCES CONFIGURE_DEPENDS "include/*.cpp")
file (GLOB_RECURSE LF_HEADERS CONFIGURE_DEPENDS "include/*.hpp")

set (LF_INCLUDE_DIRS "")
foreach (_headerFile ${LF_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND LF_INCLUDE_DIRS ${_dir})
endforeach()
list (REMOVE_DUPLICATES LF_INCLUDE_DIRS)

add_executable(liten-fraktaler main.cpp ${LF_SOURCES})
target_include_directories(liten-fraktaler PRIVATE ${LF_INCLUDE_DIRS})
target_link_libraries(liten-fraktaler png pthread)

target_compile_definitions(liten-fraktaler PUBLIC VERSION="0.8.15")
target_compile_definitions(liten-fraktaler PUBLIC FIRST_ITERATION_IS_0=true)