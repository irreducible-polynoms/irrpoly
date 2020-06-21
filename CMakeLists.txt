cmake_minimum_required(VERSION 3.8)
project(irrpoly VERSION 2.3.0 LANGUAGES C CXX) # remember to bump VERSION here

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" "${CMAKE_MODULE_PATH}")

add_library("${CMAKE_PROJECT_NAME}" INTERFACE)
target_include_directories("${CMAKE_PROJECT_NAME}"
    INTERFACE "${CMAKE_SOURCE_DIR}/include")

add_subdirectory(examples)
add_subdirectory(tests)
