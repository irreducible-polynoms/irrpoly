# add PVS-Studio analyze target
include(pvs-studio)
pvs_studio_add_target(TARGET "${PROJECT_NAME}.analyze" ALL
        OUTPUT FORMAT errorfile RECURSIVE ANALYZE "${PROJECT_NAME}")
