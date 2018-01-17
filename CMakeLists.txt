cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(find-peaks)

include_directories(.)

add_executable (example example.cpp Util.cpp)
target_link_libraries (example)

