# добавляет запуск PVS-Studio
include(pvs-studio)
pvs_studio_add_target(TARGET "${PROJECT_NAME}_analyze" ALL
        OUTPUT FORMAT errorfile RECURSIVE ANALYZE ${PROJECT_NAME})
