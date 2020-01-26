# сохраняем имя текущей папки в переменную CURRENT_DIR
# (текущей для CMakeLists.txt, который include-ит данный файл)
get_filename_component(CURRDIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
# все пробелы заменяем нижними подчёркиваниями
string(REPLACE " " "_" CURRDIR ${CURRDIR})
