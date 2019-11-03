# версия 3.8 требуется для совместимости с Visual Studio 2017
cmake_minimum_required(VERSION 3.8)
project(irrpoly)

# используется версия C++ 17
set(CMAKE_CXX_STANDARD 17)

# все файлы с расширениями cpp и hpp требуется передать компилятору
file(GLOB SOURCES *.cpp *.hpp)
# будем компилировать исполняемый файл, а не библиотеку
add_executable("${PROJECT_NAME}" "${SOURCES}")

# для работы программы используется библиотека pthread
# если она не найдена - используется std::thread
find_package(Threads)
if(CMAKE_USE_PTHREADS_INIT)
    # если pthread найден - препроцессор устанавливает #define PTHREAD
    add_definitions(-DPTHREAD)
    # выполняем линковку с библиотекой pthread
    target_link_libraries(${PROJECT_NAME} Threads::Threads)
endif(CMAKE_USE_PTHREADS_INIT)
