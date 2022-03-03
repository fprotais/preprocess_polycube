cmake_minimum_required(VERSION 2.8)
project(preprocess)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_subdirectory(ultimaille)
include_directories(ultimaille ultimaille/ext)


add_executable(preprocess main.cpp)
target_link_libraries(preprocess ultimaille)




