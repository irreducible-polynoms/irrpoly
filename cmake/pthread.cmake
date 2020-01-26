# для работы программы используется библиотека pthread
# если она не найдена - используется std::thread
find_package(Threads)
if(CMAKE_USE_PTHREADS_INIT)
    # если pthread найден - препроцессор устанавливает #define PTHREAD
    add_definitions(-DPTHREAD)
    # выполняем линковку с библиотекой pthread
    target_link_libraries(${PROJECT_NAME} Threads::Threads)
endif(CMAKE_USE_PTHREADS_INIT)
