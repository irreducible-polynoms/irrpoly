# версия 3.8 требуется для совместимости с Visual Studio 2017
cmake_minimum_required(VERSION 3.8)

# подключаем модуль currdir.cmake, сохраняющий имя текущей папки в переменную CURRDIR
include(currdir)
# устанавливаем имя проекта
project(${PROJECT_NAME}_${CURRDIR})

# подключаем модуль sources.cmake, находящий все исходники в текущей папке
# и добавляющий исполняемый файл с именем текущего проекта
include(sources)
