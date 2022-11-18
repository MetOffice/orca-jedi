
list( APPEND orcajedi_test_input
  basic.yaml
  ostia_seaice_obs_hofx3d.yaml
  hofx_odb_ice.yaml
  hofx_nc_ice.yaml
  hofx3d_nc_sst.yaml
  hofx3d_nc_potm.yaml
  hofx3d_nc_prof_2vars.yaml
  odb_ice_query.yaml
  odb_ice_mapping.yaml
  test_name_map.yaml
)

foreach(FILENAME ${orcajedi_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)
