# saves current directory name into CURRENT_DIR variable
# (current for CMakeLists.txt that includes this file)
get_filename_component(CURRDIR "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
# replace all spaces with underlines
string(REPLACE " " "_" CURRDIR "${CURRDIR}")
