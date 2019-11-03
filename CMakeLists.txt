cmake_minimum_required(VERSION 3.8)
project(irrpoly)

set(CMAKE_CXX_STANDARD 17)

file(GLOB SOURCES *.cpp *.hpp)
add_executable("${PROJECT_NAME}" "${SOURCES}")

find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_NAME} ${CMAKE_THREAD_LIBS_INIT})
