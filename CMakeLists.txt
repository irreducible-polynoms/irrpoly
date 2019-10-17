cmake_minimum_required(VERSION 3.8)
project(irrpoly CXX)

set(CMAKE_CXX_STANDARD 17)

file(GLOB SOURCES *.cpp *.hpp)
add_executable("${PROJECT_NAME}" "${SOURCES}" polynomial.hpp)
