# collect all file names ending with .c, .h, .cpp, .hpp into SOURCES variable
file(GLOB SOURCES *.cpp *.hpp)
# set directories where binaries will be placed
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${PROJECT_SOURCE_DIR}/bin/debug")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PROJECT_SOURCE_DIR}/bin/release")
# add build target
add_executable("${PROJECT_NAME}" "${SOURCES}")
target_link_libraries("${PROJECT_NAME}" "${CMAKE_PROJECT_NAME}")
# set executable name
get_filename_component(CURRDIR "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${CURRDIR}")
