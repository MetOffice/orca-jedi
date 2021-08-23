
list( APPEND orcajedi_test_input
  basic.yaml
)

foreach(FILENAME ${orcajedi_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

execute_process(COMMAND ncgen -4 -o ${CMAKE_CURRENT_BINARY_DIR}/simple_nemo.nc ${CMAKE_CURRENT_SOURCE_DIR}/simple_nemo.cdl)

execute_process(COMMAND ncgen -4 -o ${CMAKE_CURRENT_BINARY_DIR}/orca2_t_nemo.nc ${CMAKE_CURRENT_SOURCE_DIR}/orca2_t_nemo.cdl)

execute_process(COMMAND ncgen -4 -o ${CMAKE_CURRENT_BINARY_DIR}/ostia_seaice_obs.nc ${CMAKE_CURRENT_SOURCE_DIR}/ostia_seaice_obs.cdl)
