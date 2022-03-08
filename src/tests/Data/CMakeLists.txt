
list( APPEND orcajedi_test_data
ostia_seaice_obs.odb
hofx_sst_obs.nc
orca2_t_bkg_var.nc
orca2_t_coords.nc
orca2_t_nemo.nc
ostia_seaice_obs.nc
simple_nemo.nc
create_orca2_t_bkg.nco
create_orca2_t_bkg_var.nco
hofx_sst_obs.cdl
orca2_t_coords.cdl
ostia_seaice_obs.cdl
ostia_seaice_obs.csv
simple_nemo.cdl
hofx_two_vars_obs.nc
)

foreach(FILENAME ${orcajedi_test_data})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)
