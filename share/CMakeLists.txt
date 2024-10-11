# (C) British Crown Copyright 2024 Met Office

list( APPEND shared_files
    odb_col_chl_query_filter_null.yaml
)

set( destination share/orca-jedi/odb-query )

install(
    FILES       ${shared_files}
    DESTINATION ${destination}
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)

file( MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${destination} )

foreach( file ${shared_files} )
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "copy_if_different"
        "${CMAKE_CURRENT_SOURCE_DIR}/${file}"
        "${CMAKE_BINARY_DIR}/${destination}/${file}")
endforeach()
