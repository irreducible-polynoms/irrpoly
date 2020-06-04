# collect paths of all files in current directory into CHILDREN variable
file(GLOB CHILDREN "${CMAKE_CURRENT_SOURCE_DIR}/*")
foreach(CHILD ${CHILDREN})
    if(IS_DIRECTORY "${CHILD}")
        get_filename_component(CHILD ${CHILD} NAME) # get filename from full path
        add_subdirectory("${CHILD}")
    endif()
endforeach()
