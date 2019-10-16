cmake_minimum_required(VERSION 3.8)
project(irrpoly CXX)

set(CMAKE_CXX_STANDARD 17)

file(GLOB SOURCES *.cpp)
add_executable("${PROJECT_NAME}" "${SOURCES}")

# Boost 1.71.0
include_directories(include)
