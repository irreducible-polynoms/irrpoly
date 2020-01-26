# имена всех файлов с расширениями .c, .h, .cpp, .hpp помещаем в переменную SOURCES
file(GLOB SOURCES *.cpp *.hpp)
# устанавливаем папку, в которую будет помещён исполняемый файл
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${PROJECT_SOURCE_DIR}/bin/debug")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PROJECT_SOURCE_DIR}/bin/release")
# создаём исполняемый файл
add_executable(${PROJECT_NAME} "${SOURCES}")
# задаём имя исполняемого файла
get_filename_component(CURRDIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${CURRDIR}")
