cmake_minimum_required(VERSION 3.10)
project(PETPVCv2)

set(CMAKE_CXX_STANDARD 14)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(PETPVCv2 main.cpp)
target_link_libraries(PETPVCv2 ${ITK_LIBRARIES})