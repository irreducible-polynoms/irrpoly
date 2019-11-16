# версия 3.8 требуется для совместимости с Visual Studio 2017
cmake_minimum_required(VERSION 3.8)
project(irrpoly VERSION 1.0.0 LANGUAGES C CXX)

# используется версия C++ 17
set(CMAKE_CXX_STANDARD 17)

# указываем, где искать вспомогательные модули Cmake
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# исходники библиотеки
include_directories("${CMAKE_SOURCE_DIR}/include")

# примеры использования библиотеки
add_subdirectory(examples)

# тесты
add_subdirectory(tests)
