# (C) British Crown Copyright 2024 Met Office

list( APPEND orcajedi_test_input
  hofx_odb_sic.yaml
  hofx_nc_sic.yaml
  hofx_nc_sst.yaml
  hofx_nc_ssh.yaml
  hofx_nc_ssh_checkerboard.yaml
  hofx_nc_ssh_eorca025.yaml
  hofx_nc_ssh_amm1.yaml
  hofx_nc_potm_amm1r.yaml
  hofx3d_nc_sic.yaml 
  hofx3d_nc_sst.yaml
  hofx3d_nc_potm.yaml
  hofx3d_nc_prof_2vars.yaml
  hofx3d_nc_ssh_amm1.yaml
  hofx3d_nc_ssh_amm1r.yaml
  odb_sic_query.yaml
  odb_sic_mapping.yaml
  test_name_map.yaml
  dirac.yaml
  3dvar_sic.yaml
)

foreach(FILENAME ${orcajedi_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)
