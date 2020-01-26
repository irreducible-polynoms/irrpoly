# созраняем пути ко всем файлам в текущей папке в переменную CHILDREN
file(GLOB CHILDREN "${CMAKE_CURRENT_SOURCE_DIR}/*")
# проверяем все полученные имена файлов
foreach(CHILD ${CHILDREN})
    # если папка
    if(IS_DIRECTORY "${CHILD}")
        # получаем из полного пути только имя
        get_filename_component(CHILD ${CHILD} NAME)
        # добавляем подпроект
        add_subdirectory("${CHILD}")
    endif()
endforeach()
