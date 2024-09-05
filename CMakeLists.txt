cmake_minimum_required(VERSION 3.29)
project(AIS4104-assignments)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(test1)
add_subdirectory(test2)
add_subdirectory(assignment1)
add_subdirectory(assignment2)
add_subdirectory(math)

add_dependencies(assignment2 math)