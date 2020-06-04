set(THREADS_PREFER_PTHREAD_FLAG ON) # prefer POSIX threads if possible
find_package(Threads REQUIRED)
target_link_libraries("${PROJECT_NAME}" Threads::Threads)
