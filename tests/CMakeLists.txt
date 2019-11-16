# версия 3.8 требуется для совместимости с Visual Studio 2017
cmake_minimum_required(VERSION 3.8)

# подключаем модуль currdir.cmake, сохраняющий имя текущей папки в переменную CURRDIR
include(currdir)
# устанавливаем имя проекта
project(${PROJECT_NAME}_${CURRDIR})

# подключаем модуль subdirs.cmake, который добавляет в проект все папки в текущем каталоге
include(subdirs)
