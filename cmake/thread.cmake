# предпочитаем использовать POSIX threads если есть такая возможность
set(THREADS_PREFER_PTHREAD_FLAG ON)
# ищем библиотеку, предоставляющую возможность работать с потоками
find_package(Threads REQUIRED)
# линкуемся с найденной библиотекой
target_link_libraries(${PROJECT_NAME} Threads::Threads)
