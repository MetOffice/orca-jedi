
list( APPEND orcajedi_test_output
  test_ostia_seaice_obs_hofx3d.ref
)

foreach(FILENAME ${orcajedi_test_output})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)
